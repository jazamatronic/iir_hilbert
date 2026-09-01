# Bug: `iir_hilbert_even.m` produces a non-functional Hilbert transformer for roughly half of its valid input range

## Status
Diagnosed and fix verified in a standalone scratch script (not yet applied to the repo).
One open caveat at very low filter order (see "Known gap" below) — verify before relying on `N=4`.

## Affected file
`iir_hilbert_even.m`, function `iir_hilbert_even(wc, d_phi)`, pole-splitting loop at
[iir_hilbert_even.m:75-88](iir_hilbert_even.m#L75-L88):

```matlab
alpha0 = [];
alpha1 = [];
for i = 1:N;
  if (rem(i, 2))
    alpha0 = [alpha0 alpha_hb_all(i)];
  else
    alpha1 = [alpha1 alpha_hb_all(i)];
  end
endfor

A = poly(-alpha0);
B = fliplr(A);
Ai = poly(-alpha1);
Bi = fliplr(Ai);
```
where `alpha_hb_all = imag(p)` and `p` are the poles returned by `ellip(N, rp_db, rs_db, wc)` at [iir_hilbert_even.m:68-71](iir_hilbert_even.m#L68-L71).

## Symptom
The function is supposed to return two IIR all-pass branches `(B,A)` and `(Bi,Ai)` whose outputs differ in phase by ~90° (± `d_phi`) across the whole band except right at DC/Nyquist. For roughly half of achievable filter orders, the two branches instead come out **nearly in phase** (phase difference ≈ 0 rad instead of ≈ ±π/2) — i.e. the returned filters are not a Hilbert transformer at all. Both filters remain individually stable, so nothing errors; the failure is silent and only visible by actually measuring the phase response.

## Reproduction
```
octave-cli --eval "
pkg load signal;
[B, A, Bi, Ai] = iir_hilbert_even(0.47, 0.02);
[h0, w] = freqz(B, A, 400);
[h1, w] = freqz(Bi, Ai, 400);
diffp = mod(arg(h1) - arg(h0), 2*pi);
diffp(diffp > pi) -= 2*pi;
printf('phase diff at mid-band = %.4f (should be near +-pi/2 = %.4f)\n', diffp(200), pi/2);
"
```
Output: phase diff at mid-band ≈ `-0.009` rad — essentially zero, not `±1.5708`.

`wc = 0.47` is well within the documented valid domain (`0 < wc < 0.499*pi`; the sibling function `iir_hilbert_odd.m` even lists `wc = 0.44` as its doc example). The repo's own shipped default (`wc = 0.495, d_phi = 0.012235`, used in `chorus.m` / `stereo_phaser.m` / `test_sos.m`) happens to avoid the bug — see "Root cause" for why.

## Root cause
`ellip()` returns its pole vector `p` (length `N`) as two mirrored halves:
- `p(1 : N/2)`  = ascending-magnitude poles with **positive** imaginary part
- `p(N/2+1 : N)` = the same magnitudes, **negative** imaginary part, same ascending order

So `alpha_hb_all = imag(p)` is `[+m1, +m2, ..., +mL, -m1, -m2, ..., -mL]` where `L = N/2` and `m1 < m2 < ... < mL`.

The branch-assignment loop alternates by **raw array index parity** (`rem(i,2)`) over this full length-`N` array. Follow one matched pair `(+m_k, -m_k)`:
- `+m_k` is at index `k`.
- `-m_k` is at index `L + k`.

Whether these two land in the *same* branch or *opposite* branches depends entirely on the parity of `L`:
- **`L` odd**: `L + k` always has opposite parity from `k` → `+m_k` and `-m_k` always end up in different branches (one branch gets `+m_k`, the other gets `-m_k`). This is the correct pole-interlacing behavior the docstring's Fig 13.8 / Fig 8.44 references describe.
- **`L` even**: `L + k` always has the *same* parity as `k` → `+m_k` and `-m_k` land in the **same** branch. That branch now contains a self-conjugate pair for that magnitude level, which cancels exactly the odd/even asymmetry between the two branches that produces the 90° phase relationship — collapsing the phase difference toward 0.

In short: the loop's correctness is an accident of `ellip()`'s return-value layout combined with `N/2` happening to be odd. It is not a property of the filter design math. Since the repo's default parameters happen to produce `N=14` (`L=7`, odd), the bug has stayed hidden in normal use.

Verified failure/success pattern by holding `wc`/`rp_db`/`rs_db` fixed and only varying `N`:

| N  | L=N/2 | L parity | mid-band phase diff (rad) | phase quadrature? |
|----|-------|----------|---------------------------|--------------------|
| 4  | 2     | even     | ≈ -0.008                  | ✗ broken |
| 6  | 3     | odd      | ≈ -1.590                  | ✓ ok |
| 8  | 4     | even     | ≈ -0.009                  | ✗ broken |
| 10 | 5     | odd      | ≈ -1.593                  | ✓ ok |
| 12 | 6     | even     | ≈ -0.009                  | ✗ broken |
| 14 | 7     | odd      | ≈ -1.591 (repo default)   | ✓ ok |
| 16 | 8     | even     | ≈ -0.009                  | ✗ broken |

## Fix
Don't alternate index parity over the raw, doubled-length `imag(p)` array. Alternate **sign** over the `L` unique magnitude levels instead — this is parity-of-`L`-independent by construction.

Replace [iir_hilbert_even.m:75-88](iir_hilbert_even.m#L75-L88) with something equivalent to:

```matlab
L = N / 2;
m = alpha_hb_all(1:L);   % ascending positive magnitudes, one per pole pair

alpha0 = zeros(1, L);
alpha1 = zeros(1, L);
for i = 1:L
  if (rem(i, 2))
    alpha0(i) =  m(i);
    alpha1(i) = -m(i);
  else
    alpha0(i) = -m(i);
    alpha1(i) =  m(i);
  end
endfor

A = poly(-alpha0);
B = fliplr(A);
Ai = poly(-alpha1);
Bi = fliplr(Ai);
```

The key change: iterate `i = 1:L` (one iteration per magnitude level / pole pair) and explicitly assign `+m(i)` to one branch and `-m(i)` to the other, alternating which branch gets which sign as `i` increases. This reproduces the correct interlacing for every pair regardless of whether `L` is odd or even, because it no longer depends on comparing the parity of two different index ranges (`k` vs `L+k`) — there is only one index range now.

This relies on the empirically-observed fact that `p(1:N/2)` from `ellip()` is the ascending positive-imaginary half; if that ordering assumption is ever in doubt, it's safer to compute `m = sort(abs(imag(p(imag(p) > 0))))` explicitly rather than slicing `alpha_hb_all(1:L)`.

## Verification of the fix (standalone, not yet applied to the repo)
Re-running the same `N` sweep with the corrected assignment logic:

| N  | L | maxdev from ±π/2 across passband (old code) | maxdev (fixed code) |
|----|---|----------------------------------------------|----------------------|
| 4  | 2 | 2.45 (broken) | 0.545 — **still elevated, see gap below** |
| 6  | 3 | 0.034 (was fine) | 0.034 (unchanged) |
| 8  | 4 | 2.83 (broken) | 0.025 (fixed) |
| 10 | 5 | 0.028 (was fine) | 0.028 (unchanged) |
| 12 | 6 | 2.82 (broken) | 0.035 (fixed) |
| 14 | 7 | 0.030 (was fine, matches repo default) | 0.030 (unchanged) |
| 16 | 8 | 2.82 (broken) | 0.037 (fixed) |

All previously-broken even-`L` cases (`L=4,6,8`) now match the error band of the cases that already worked. The repo's default (`L=7`) is unaffected by the change, as expected.

## Known gap — needs follow-up before trusting low order
`N=4` (`L=2`) still shows a larger deviation (0.545 rad) than the rest even after the fix. Not yet root-caused. Two candidate explanations, neither confirmed:
1. A genuine remaining defect in the pairing logic specific to `L=2` (only one alternation, no "settling" — worth checking against the actual Mitra/Gray-Markel reference figures for the `L=2` case specifically).
2. A measurement artifact: low-order elliptic filters have a proportionally wider transition band, so the fixed-width edge-exclusion window used when computing `maxdev` (excluding the outer 10% of the band near DC/Nyquist) may still be catching genuine transition-region phase droop rather than a real algorithm failure.

Recommended next step: re-measure `N=4` with a wider edge-exclusion margin (e.g. 25–30%) and/or compare directly against a known-good `L=2` reference (e.g. hand-derived from Gray & Markel [3] or Mitra Ch.13 Fig 13.8) before shipping the fix for low-order use cases. Mid-to-high order (`N ≥ 6`, matching realistic `d_phi`/`wc` combinations for audio use) is unambiguously fixed by the change above.

## Suggested follow-up hardening (optional, not required to fix the bug)
Add a self-check at the end of `iir_hilbert_even.m` that computes the phase difference between the two returned branches at a handful of frequencies and asserts it's within `d_phi` of `±π/2` — this class of bug (silent, stable, wrong) is exactly what a runtime assertion would have caught immediately instead of requiring numerical excavation.

## TODO — narrow validity region for `wc`, independent of the pole-splitting bug above

Investigated during the Python port (`iir_hilbert_even.py`). Not the same bug as above — this shows up even with the pole-splitting fix applied. **Not yet root-caused; no fix implemented.** Recorded here so it isn't rediscovered from scratch.

### Symptom
Even with the corrected pole-splitting logic, the design fails to meet its own `d_phi` spec (`maxdev >= d_phi`, measured the same way as above) for a substantial middle range of `wc`, despite `ellip()` returning a filter without error. The repo default (`wc=0.495`) and the sibling `iir_hilbert_odd.m` docstring's example (`wc=0.44`) sit on opposite sides of this boundary — `wc=0.44` mostly fails.

### Tests run (exploratory, not committed as code — python3 + scipy 1.11.4, hand-rolled order/pole logic implemented ad hoc to match the `.m` file; corrected pole-splitting per the "Fix" section above)

1. **Order-estimation sweep, hand-rolled formula (`.m` lines 53-66) vs `scipy.signal.ellipord`.** 45 points (`wc` in [0.30, 0.499] × `d_phi` in [~0.004π, 0.03π]). `ellipord` never picked a lower order than the hand-rolled formula (0/45), matched it in 32/45, picked higher (usually +2, sometimes +4) in 13/45. Where they diverged, the hand-rolled formula's lower order sometimes met spec while `ellipord`'s higher order did not (e.g. `wc=0.46, d_phi=0.0393`: hand-rolled N=6 → maxdev=0.0386 ✓, ellipord N=8 → maxdev=0.0396 ✗). **Conclusion: keep the hand-rolled order formula; `ellipord` optimizes magnitude ripple, not this construction's phase-quadrature property, and is not a safe drop-in.**

2. **Full sweep, hand-rolled formula only**, same 45-point grid, tabulating pass/fail (`maxdev < d_phi`) by `wc`:

   | wc | pass rate (of 5 `d_phi` values) |
   |---|---|
   | 0.30, 0.35 | 0/5 |
   | 0.40 | 2/5 |
   | 0.44 | 1/5 |
   | 0.46, 0.48, 0.49, 0.495, 0.499 | 5/5 |

   Pass/fail is **not** a function of `wc` alone — at `wc=0.44` the one passing `d_phi` value is surrounded by failing values at both tighter and looser tolerances, so a simple `wc`-only input-range guard can't cleanly separate valid from invalid designs. (This is why the earlier "Suggested follow-up hardening" self-check — measuring actual `maxdev` post-design — is the right guard, not a static `wc` floor.)

3. **Does iteratively increasing `N` recover a failing design?** Tested at `wc=0.44` for all 4 failing `d_phi` values from (2), sweeping `N` from the hand-rolled value up to +24 (step 2), `rp_db`/`rs_db`/`wc` held fixed. **No — it makes things worse.** `maxdev` roughly triples at `N+2` then oscillates at an elevated plateau for 20+ more poles, never trending back toward spec. Ruled out as a runtime remediation strategy (e.g. "catch the self-check failure, retry with larger N").

4. **Root-cause lead: discarded pole real part.** The `.m` code takes `alpha_hb_all = imag(p)`, discarding the real part of every pole outright (comment: "close but don't quite get there — this cleans it up"). Measured `max(|real(p)|) / max(|imag(p)|)` at `d_phi=0.012235` across `wc`:

   | wc | N | re/im ratio | passes spec? |
   |---|---|---|---|
   | 0.30 | 4 | 0.092 | ✗ |
   | 0.35 | 6 | 0.178 | ✗ |
   | 0.40 | 6 | 0.007 | ✓ |
   | 0.44 | 8 | 0.056 | ✗ |
   | 0.46 | 8 | 0.008 | ✓ |
   | 0.48 | 10 | 0.008 | ✓ |
   | 0.495 | 14 | 0.008 | ✓ |

   Ratio tracks pass/fail closely (elevated ~0.05-0.18 in failing cases vs. flat ~0.007-0.009 in passing cases). This is a genuine lead, not yet a confirmed mechanism.

5. **Ruled out as the cause of (4):** (a) suboptimal integer `N` choice — swept `N` both above and below the hand-rolled value at `wc=0.44`; the hand-rolled `N` was already the local minimum of both `maxdev` and the re/im ratio, every neighbor was worse. (b) `rs_db = ceil(mag2db(1/rs))` integer-dB rounding — reran with exact (unrounded) `rs_db`; changed `maxdev` by <10%, and in the wrong direction to explain the effect. Neither is a tunable knob that fixes this within the current construction.

### Working hypothesis
The elevated pole real-part is a structural property of what `scipy.signal.ellip()` (and by extension Octave's `ellip()`, same underlying elliptic prototype math) returns for these particular `(wc, d_phi, N)` combinations — not a rounding or parameter-selection artifact fixable by retuning within the current "generic `ellip()`, then discard `real(p)`" approach.

### Possible real fix (substantial follow-up work, not a quick patch)
Replace the "call generic `ellip()`, reverse-engineer allpass coefficients by discarding pole real parts" approach with a proper closed-form derivation that computes the allpass section coefficients directly from `wp`/`ws`/`rp`/`rs` via Jacobi elliptic functions — the actual Gray-Markel / Ansari / Renfors-Neuvo halfband recursive allpass design method that references [1]-[3] (and this design's own docstring) are gesturing at, rather than this file's apparent shortcut of extracting coefficients from a generic elliptic lowpass's poles.

- **Building blocks already available, no need to implement from scratch:** `scipy.special.ellipj` (Jacobi `sn`/`cn`/`dn`), `scipy.special.ellipk`/`ellipkm1`/`ellipkinc`, `scipy.special.ellipe`/`ellipeinc`. Note scipy's convention is parameter `m = k²`, not modulus `k` — most DSP literature states formulas in `k`/`k'`, so formulas need `m = k**2` conversion when transcribing.
- **Not available packaged anywhere found:** the halfband-specific allpass coefficient formula itself (as opposed to the elliptic-function primitives it's built from). Would need to be transcribed from the literature.
- **Literature pointers:** Mitra & Kaiser, *Handbook for Digital Signal Processing*, Ch. 13 (already cited as [1] in this file) and Mitra, *Digital Signal Processing: A Computer-Based Approach* (already cited as [2]) are the references already in hand and the first place to look for the closed-form allpass coefficient equations (as opposed to the order-estimation series already ported from Ch. 13/[2] 13.110). R. Ansari and separately M. Renfors / Y. Neuvo are names associated with recursive/allpass-based halfband elliptic filter design in the DSP literature and worth a literature search for the specific closed-form derivation — exact paper titles/years not verified in this session, don't take them as confirmed citations.
- Scope: this is a different algorithm from what's in `iir_hilbert_even.m`, not a bugfix to it. Treat as a separate task, not part of the current port + pole-splitting-bug fix.

### Practical takeaway for the current port
Don't try to widen the valid `wc` range as part of this port. Ship the pole-splitting fix and the `maxdev` self-check (see "Suggested follow-up hardening" above); the self-check will correctly reject `(wc, d_phi)` combinations that fall in this narrower-than-documented region regardless of root cause. The repo's actual usage (`wc=0.495` in `chorus.m`/`test_sos.m`) is well inside the reliable region. Revisit only if a real use case needs `wc` in roughly the 0.40-0.46 band.
