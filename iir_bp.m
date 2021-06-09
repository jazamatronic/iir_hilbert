function [B, A] = iir_bp(fs, fc, bw)
% 2nd order IIR filter 
%
% fs = sampling frequency
%
% fc = center frequency in Hz
%
% bw = 3db bandwidth in Hz
%
% B,A = IIR filter coefs that should be able to go straight into a biquad
%
% Handbook for Digital Signal Processing Ch 9.
%
% [1] Sanjit K. Mitra,
%     Digital Signal Processing a Computer Based Approach
%


% From fs and fc/bw derive normalized frequencies
wc = 2 * pi * fc / fs;
bl = fc - (bw / 2);
bh = fc + (bw / 2);
wl = 2 * pi * bl / fs;
wh = 2 * pi * bh / fs;
ww = wh - wl;

% [1] (9.30)
alpha = (1 - tan(ww / 2)) / (1 + tan(ww / 2));
beta  = cos(wc);

% [1] (9.28)
g = (1 - alpha) / 2;
B = g .* [1, 0, -1];
A = [1, -beta * (1 + alpha), alpha];


% Test
%  ts = 1/fs;
%  t = [0:ts:0.2-ts];
%  f1 = 220;
%  s1 = sin(t * f1 * 2 * pi);
%  
%  figure(1) 
%  freqz(B, A);
%  out = filter(B, A, s1);
%  figure(2)
%  plot(t, s1, "-b;in;", t, out, "-r;out;")
%  title("220Hz sine test - implementation")
%  xlim([0.1 0.125])
%  ylim([-1.1 1.1])
