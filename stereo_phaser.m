% Stero Phaser Effect
%
% Example of using a Hilbert Transformer as part of Single Sideband Modulator (SSB)
%
% Dan Harris, Edgar Berdahl, and Jonathan S. Abel
% An Infinite Impulse Response (IIR) Hilbert Transformer Filter Design Technique for Audio
% AES 129th Convention Proceedings 
clear all

% Enter some audio
%[y, fs] = audioread("<mono_wave_file_or_similar>");

ts = 1/fs;

wc = 0.495;

% dp,ds pairing is roughly equivalent to the given d_phi 
dp = 0.0001;
ds = 0.1;
d_phi = 0.012235;

%Pick your favourite hilbert transformer
[B, A, Bi, Ai] = iir_hilbert_even(wc, d_phi);
%[B, A, Bi, Ai] = iir_hilbert_odd(wc, dp, ds);
out0 = filter(B, A, y);
out1 = filter(Bi, Ai, y);

t = [0:ts:ts * (length(y) -1)];
f1 = 0.2;
s1 = sin(t * f1 * 2 * pi)';
c1 = cos(t * f1 * 2 * pi)';
y0 = out0 .* c1;
y1 = out1 .* s1;
y0_out = y0 .+ y1;
y1_out = -y1 .+ y0;
phased(:,1) = (y .+ y0_out) / sqrt(2);
phased(:,2) = (y .+ y1_out) / sqrt(2);
plot(t, y, "-r;in;", t, phased, "-b;out;")
% Playback the original follwed by the stereo phase panning version
player = audioplayer(y, fs);
playblocking(player); 
player = audioplayer(phased, fs);
playblocking(player);
