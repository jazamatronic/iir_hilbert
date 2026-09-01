"""
pytest tests for iir_hilbert_even.py

Covers:
- the pole-splitting bug documented in HILBERT_EVEN_BUG.md (regression, both isolated and end-to-end)
- the general phase-quadrature spec the design is supposed to satisfy
- the narrow-wc-validity-region limitation from the TODO section of HILBERT_EVEN_BUG.md (documented, not silently broken)
- structural/stability invariants of the returned filter coefficients

Written with assistance from Claude (Sonnet 5, model id claude-sonnet-5).
"""

import numpy as np
import pytest
from scipy.signal import freqz

from iir_hilbert_even import iir_hilbert_even


def in_band_maxdev(B, A, Bi, Ai, worN=2000):
    """Max deviation of the branch phase difference from pi/2, excluding the
    outer 10% of the band at each end (near DC/Nyquist) -- matches the
    convention used to verify the pole-splitting fix in HILBERT_EVEN_BUG.md."""
    w, h0 = freqz(B, A, worN=worN)
    _, h1 = freqz(Bi, Ai, worN=worN)
    diffp = np.angle(h1) - np.angle(h0)
    diffp = np.mod(diffp, 2 * np.pi)
    diffp[diffp > np.pi] -= 2 * np.pi
    n = len(w)
    lo, hi = int(0.10 * n), int(0.90 * n)
    return np.max(np.abs(np.abs(diffp[lo:hi]) - np.pi / 2))


# --- Bug regression: pole-splitting, isolated from the rest of the design pipeline ---

@pytest.mark.parametrize("L", [2, 3, 4, 5, 6, 7, 8])  # covers both parities of L
def test_split_poles_interlaces_regardless_of_l_parity(L):
    m = np.sort(np.random.default_rng(0).uniform(0.1, 1.0, L))
    p = np.concatenate([1j * m, -1j * m])
    alpha0, alpha1 = iir_hilbert_even.split_poles(p)
    assert not np.all(np.sign(alpha0) == np.sign(alpha0[0]))
    assert np.allclose(alpha0, -alpha1)


# --- Bug regression: end-to-end, anchored to the bug doc's own reproduction case ---

def test_bug_doc_repro_case_no_longer_collapses():
    # HILBERT_EVEN_BUG.md's exact reproduction: was ~-0.009 rad (collapsed), should be ~pi/2
    h = iir_hilbert_even(wc=0.47, d_phi=0.02)
    w, h0 = freqz(h.B, h.A, worN=400)
    _, h1 = freqz(h.Bi, h.Ai, worN=400)
    diffp = np.mod(np.angle(h1) - np.angle(h0), 2 * np.pi)
    diffp = np.where(diffp > np.pi, diffp - 2 * np.pi, diffp)
    assert abs(diffp[200]) > 1.0  # nowhere near the old ~0 collapse


# --- Spec correctness: phase quadrature over a sweep of known-good (wc, d_phi) ---

@pytest.mark.parametrize("wc, d_phi", [
    (0.495, 0.012235),   # repo default, used by chorus.m / test_sos.m
    (0.46, 0.0393),
    (0.48, 0.0628),
])
def test_meets_phase_quadrature_spec(wc, d_phi):
    h = iir_hilbert_even(wc, d_phi)
    assert in_band_maxdev(h.B, h.A, h.Bi, h.Ai) < d_phi


# --- Known limitation: narrow wc validity region (see HILBERT_EVEN_BUG.md TODO section) ---

@pytest.mark.xfail(reason="narrow wc validity region, see HILBERT_EVEN_BUG.md TODO section")
def test_wc_044_currently_fails_spec():
    wc, d_phi = 0.44, 0.012235
    h = iir_hilbert_even(wc, d_phi)
    assert in_band_maxdev(h.B, h.A, h.Bi, h.Ai) < d_phi


# --- Structural/regression anchors ---

@pytest.mark.parametrize("wc, d_phi", [
    (0.495, 0.012235),
    (0.44, 0.0393),
    (0.30, 0.0628),
])
def test_ellip_ord_always_even(wc, d_phi):
    assert iir_hilbert_even.ellip_ord(wc, d_phi) % 2 == 0


def test_repo_default_order_is_14():
    assert iir_hilbert_even.ellip_ord(0.495, 0.012235) == 14


def test_b_is_reverse_of_a():
    h = iir_hilbert_even(0.495, 0.012235)
    assert np.allclose(h.B, h.A[::-1])
    assert np.allclose(h.Bi, h.Ai[::-1])


def test_filters_are_stable():
    h = iir_hilbert_even(0.495, 0.012235)
    assert np.all(np.abs(np.roots(h.A)) < 1)
    assert np.all(np.abs(np.roots(h.Ai)) < 1)
