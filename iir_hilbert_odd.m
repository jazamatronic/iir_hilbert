function [B, A, Bi, Ai] = iir_hilbert_odd(wc, dp, ds)
% Hilbert transformer based on odd order polyphase IIR halfband filter 
% constructed from parallel all-pass sections
%
% wc = passband edge   0 < wc < 0.499 * pi
%   e.g. 0.44
%
% dp = passband mag deviation
%   e.g. 0.0001
% ds = stopband mag deviation
%   e.g. 0.1
%
% B,A,Bi,Ai = IIR filter coefs for real and imag branches
%
% Rough design approach - create an elliptic LP IIR filter with constraints that put all the poles on the IMG axis
% Then, multiply by -j to put the poles on the REAL axis 
% These poles are then divided into two branches - the h0 branch which is real, the h1 branch which is imaginary
% The phase difference between these branches, from roughly DC up to fs/2 is 90Degrees - magic
%
% Elliptic IIR design from
% 
% [1] Rashid Ansari
%   Elliptic Filter Design for a Class of Generalized Halfband Filters
%   IEEE TRANSACTIONS ON ACOUSTICS, SPEECH, AND SIGNAL PROCESSING, VOL. ASSP-33, NO. 4. OCTOBER 1985
%
% LP to AP phase difference network transform from:
%
% [2] Dan Harris, Edgar Berdahl, and Jonathan S. Abel
% An Infinite Impulse Response (IIR) Hilbert Transformer Filter Design Technique for Audio
% AES 129th Convention Proceedings 
%
%
% https://dsp.stackexchange.com/questions/37411/iir-hilbert-transformer
% https://github.com/robwasab/HalfBand
% https://dsp.stackexchange.com/questions/8692/hilbert-transform-filter-for-audio-applications-using-iir-half-band-parallel-al
% Octave signal ellip.m, ellipord.m and ncauer.m
%
%
% Jared ANDERSON
% 20210603

pkg load signal

wp = wc * pi;

% constraint from [1] is wp = pi - ws
ws = pi - wp;

As  = -20 * log10(ds); %(19) [2] - unused
d   = min(ds, sqrt(2 * dp - dp^2));
del = (1 - d^2) / d^2;

% order estimation
k  = tan(wp / 2)^2; 
kp = sqrt(1 - k^2); 
p  = 0.5 * ((1 - sqrt(kp)) / (1 + sqrt(kp))); % equiv to q0
q  = p + 2 * p^5 + 15 * p^9 + 150 * p^13;
N  = ceil((2 * log(4 * del)) / -log(q));

% Make sure N is odd
if (mod(N + 1, 2))
  N++;
endif

L = (N - 1) / 2;

% from ncauer
wi=zeros(1, L);
for i = 1:L
  soma1 = 0;
  for m = 0:30
    soma1 = soma1 + 2 * q^(1/4) * ((-1)^m * q^(m * (m + 1)) * sin(((2 * m + 1) * pi * i) / N));
  endfor
  soma2=0;
  for m=1 : 30
    soma2 = soma2 + 2 * ((-1)^m * q^(m^2) * cos((2 * m * pi * i) / N));
  endfor
  wi(i)=(soma1/(1+soma2));
endfor

vi = sqrt((1 - (k .* (wi.^2))) .* (1 - (wi.^2) / k));

% from (28,29) [2]
ci = (2 * vi) ./ (1 .+ wi.^2);
alpha = (2 - ci) ./ (2 + ci);

% Sort the coefficients into top and bottom branches
alpha0 = [];
alpha1 = [];
for i = 1:L
  if (rem(i, 2))
    alpha0 = [alpha0 alpha(i)];
  else
    alpha1 = [alpha1 alpha(i)];
  end
endfor

% generate num/den poly coeffs
% insert 0s to make it all powers of z^-2
den0 = poly(alpha0);
A = zeros(1, length(alpha0) * 2);
A(1:2:length(den0) * 2) = den0;
B = fliplr(-A);

den1 = poly(alpha1);
Ai = zeros(1, length(alpha1) * 2);
Ai(1:2:length(den1) * 2) = den1;
% add unit delay to bottom branch
Ai(length(Ai) + 1) = 0;
Bi = fliplr(-Ai);

% Test
%  fs = 48e3;
%  ts = 1/fs;
%  t = [0:ts:0.2-ts];
%  f1 = 220;
%  s1 = sin(t * f1 * 2 * pi);
%  
%  wc = 0.495;
%  dp  = 0.0001;
%  ds  = 0.1;
%  [B, A, Bi, Ai] = iir_hilbert_odd(wc, dp, ds)
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
