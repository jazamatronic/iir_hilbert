function [B, A, Bi, Ai] = iir_hilbert_even(wc, d_phi)
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
% 20210603
pkg load signal

% wc = passband edge   0 < wc < 0.499 * pi
%   e.g. 0.44
wp = wc * pi;
ws = pi - wp;

% [1] (13.30) - stopband attenuation = sin(d_phi/2) where d_phi = max deviation from phase quadrature
% convert to db and change sign for atten
rs = sin(d_phi / 2); % For normalized mag specifications: Maximum stopband ripple = 1 / A. This gives minimum stopband attenuation of -20 * log10(1/A)
rs_db = ceil(mag2db(1/rs));

% [1] (13.26) - passband max atten and stopband atten are related by (1 - rp)^2 + rs^2 = 1
% This constraint, along with wp + ws = pi ensure that poles lie on the img axis - although maybe the rp/rs constraint is for odd order filters?
% In the handbook for dsp example I can't figure out how they they get 0.0316dB for rp form d_phi and rs
% Maybe it's a typo - 0.0016 is rp_db calculated using (13.26) and those values in ellip get pretty close
rp = 1 - sqrt(1 - rs^2);

% rp_db = 0.0316 % typo in [1]?
rp_db = mag2db(1/(1-rp));

r = tan(0.5 * wp) / tan(0.5 * ws);
rprime = sqrt(1 - r^2);
rho0 = 0.5 * ((1 - sqrt(rprime)) / (1 + sqrt(rprime)));
rho = rho0 + 2 * rho0^5 + 15 * rho0^9 + 150 * rho0^13;

% [2] 13.110
dels = rs;
D = ((1 - dels^2) / dels^2)^2;
N = round(log10(16 * D) / log10(1 / rho)); % 

% Make sure N is even
if (mod(N, 2))
  N++;
endif

[z,p,g] = ellip(N, rp_db, rs_db, wc); 

% Poles are meant to lie on the imag axis - These are close but don't quite get there - this cleans it up
alpha_hb_all = imag(p);

% Sort the coefficients into top and bottom branches
% See [1] (Fig 13.8) or [2] (Fig 8.44) for pole interlacing
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

% Test
%  fs = 48e3;
%  ts = 1/fs;
%  t = [0:ts:0.2-ts];
%  f1 = 220;
%  s1 = sin(t * f1 * 2 * pi);
%
%  wc = 0.495;
%  d_phi = 0.012235;
%  [B, A, Bi, Ai] = iir_hilbert_even(wc, d_phi);
%  
%  out0 = filter(B, A, s1);
%  out1 = filter(Bi, Ai, s1);
%  figure(1)
%  plot(t, out0, "-b;h0 branch;", t, out1, "-r;h1 branch;")
%  title("220Hz sine test - implementation")
%  xlim([0.1 0.125])
%  ylim([-1.1 1.1])
%  [h0, w0] = freqz(B, A);
%  [h1, w1] = freqz(Bi, Ai);
%  p0 = arg(h0);
%  p1 = arg(h1);
%  figure(2)
%  plot(w0 / 2, abs(mod((p1 - p0), -pi)))
%  xlabel("Normalized Freq (Rad)")
%  ylabel("Phase difference (Rad)")
%  title("Phase difference vs Freq")
%  ylim([pi / 2 - 0.1, pi / 2 + 0.1])
%  xlim([0, pi/2])

