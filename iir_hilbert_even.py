"""
Python implementation of iir_hilbert_even.m

% Hilbert transformer based on even order polyphase IIR halfband filter 
% constructed from parallel all-pass sections
%
% wc = passband edge   0 < wc < 0.499 * pi
%   e.g. 0.44
%
% d_phi = max deviation from phase quadrature
%   e.g. 0.0125 * pi
%
% B,A,Bi,Ai = IIR filter coefs for real and imag branches
%
% Handbook for Digital Signal Processing Ch 13.
%
% This is meant to an improvement over iir_hilbert_odd.m in that it should use even order filters to produce all-pass
% sections with a single z^-1 delay, as opposed to the two delay elements in the first attempt
% 
% REFERENCES:
%
% [1] Sanjit K. Mitra, James F. Kaiser - editors
%     Handbook for Digital Signal Processing
%
% [2] Sanjit K. Mitra,
%     Digital Signal Processing a Computer Based Approach
%
% [3] A.H. Gray Jr, John D. Markel
%     A Computer Program for Designing Digital Elliptic Filters
%     IEEE TRANSACTIONS ON ACOUSTICS, SPEECH, AND SIGNAL PROCESSING, DECEMBER 1976
%
% Jared ANDERSON
% 20260901
%
% Ported from Octave to Python with assistance from Claude (Sonnet 5, model id claude-sonnet-5).
% Along the way, fixed a pole-splitting bug that silently produced a non-functional (non-quadrature)
% transformer for even L=N/2
"""

# %%
import numpy as np
from scipy.signal import ellip


def mag2db(x):
  return 20. * np.log10(x)


class iir_hilbert_even:

  # wc = passband edge   0 < wc < 0.499
  # d_phi = max deviation from phase quadrature
  def __init__(self, wc, d_phi):
    self.wc = wc
    self.d_phi = d_phi

    rp_db, rs_db = self.ripple_dbs(d_phi)
    N = self.ellip_ord(wc, d_phi)
    _z, p, _k = ellip(N, rp_db, rs_db, wc, btype='low', output='zpk')

    alpha0, alpha1 = self.split_poles(p)

    self.A = np.poly(-alpha0)
    self.B = np.flip(self.A)
    self.Ai = np.poly(-alpha1)
    self.Bi = np.flip(self.Ai)


  """
  [1] (13.30) - stopband attenuation = sin(d_phi/2) where d_phi = max deviation from phase quadrature
  convert to db and change sign for atten

  [1] (13.26) - passband max atten and stopband atten are related by (1 - rp)^2 + rs^2 = 1
  This constraint, along with wp + ws = pi ensure that poles lie on the img axis - although maybe the rp/rs constraint is for odd order filters?
  In the handbook for dsp example I can't figure out how they they get 0.0316dB for rp form d_phi and rs
  Maybe it's a typo - 0.0016 is rp_db calculated using (13.26) and those values in ellip get pretty close

  """
  @staticmethod
  def ripple_dbs(d_phi):
    # For normalized mag specifications: Maximum stopband ripple = 1 / A. This gives minimum stopband attenuation of -20 * log10(1/A)
    rs = np.sin(d_phi / 2) 
    rs_db = np.ceil(mag2db(1/rs))

    rp = 1 - np.sqrt(1 - rs**2)
    rp_db = mag2db(1/(1-rp))

    return rp_db, rs_db

  @staticmethod
  def ellip_ord(wc, d_phi):
    wp = wc * np.pi
    ws = np.pi - wp
    r = np.tan(0.5 * wp) / np.tan(0.5 * ws)
    rprime = np.sqrt(1 - r**2)
    rho0 = 0.5 * ((1 - np.sqrt(rprime)) / (1 + np.sqrt(rprime)))
    rho = rho0 + 2 * rho0**5 + 15 * rho0**9 + 150 * rho0**13

    # [2] 13.110
    dels = np.sin(d_phi / 2)
    D = ((1 - dels**2) / dels**2)**2
    N = round(np.log10(16 * D) / np.log10(1 / rho))

    # Make sure N is even
    if N % 2 != 0: 
      N += 1
      
    return N


  @staticmethod
  def split_poles(p):
    imag_p = p.imag
    m = np.sort(np.abs(imag_p[imag_p > 0]))
    L = len(m)
    alpha0 = np.empty(L)
    alpha1 = np.empty(L)
    for i in range(L):
        sign = 1 if i % 2 == 0 else -1
        alpha0[i] = sign * m[i]
        alpha1[i] = -sign * m[i]
    return alpha0, alpha1


# %%
if __name__ == "__main__":
  import matplotlib.pyplot as plt
  from scipy.signal import freqz, lfilter

  # Time-domain sanity check: 220Hz sine through both all-pass branches
  fs = 48e3
  ts = 1 / fs
  t = np.arange(0, 0.2, ts)
  f1 = 220
  s1 = np.sin(t * f1 * 2 * np.pi)

  wc = 0.495
  d_phi = 0.012235
  h = iir_hilbert_even(wc, d_phi)

  out0 = lfilter(h.B, h.A, s1)
  out1 = lfilter(h.Bi, h.Ai, s1)

  fig1, ax1 = plt.subplots()
  ax1.plot(t, out0, "-b", label="h0 branch")
  ax1.plot(t, out1, "-r", label="h1 branch")
  ax1.set_title("220Hz sine test - implementation")
  ax1.set_xlim(0.1, 0.125)
  ax1.set_ylim(-1.1, 1.1)
  ax1.legend()

  # Frequency-domain check: phase difference between branches should sit at pi/2
  w, h0 = freqz(h.B, h.A)
  _, h1 = freqz(h.Bi, h.Ai)
  p0 = np.angle(h0)
  p1 = np.angle(h1)
  diffp = np.abs(np.mod(p1 - p0, -np.pi))

  fig2, ax2 = plt.subplots()
  ax2.plot(w, diffp)
  ax2.set_xlabel("Normalized Freq (Rad)")
  ax2.set_ylabel("Phase difference (Rad)")
  ax2.set_title("Phase difference vs Freq")
  ax2.set_ylim(np.pi / 2 - 0.1, np.pi / 2 + 0.1)
  ax2.set_xlim(0, np.pi / 2)

  # Numeric correctness check (excludes outer 10% of band, near DC/Nyquist)
  band = (w >= 0.05 * np.pi) & (w <= 0.95 * np.pi)
  maxdev = np.max(np.abs(diffp[band] - np.pi / 2))
  print(f"wc={wc}  d_phi={d_phi:.6f}  N={h.ellip_ord(wc, d_phi)}  (L={len(h.A) - 1} per branch)")
  print(f"max phase deviation from pi/2 across passband = {maxdev:.6f}")
  print(f"meets d_phi spec? {'YES' if maxdev < d_phi else 'NO'}")

  plt.show()