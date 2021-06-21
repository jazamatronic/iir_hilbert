%
% Parametric chorus simulator
% From Fig 14-25 of Hal Chamberlin's Musical Applications of Microprocessors - Second Edition
%
% Jared ANDERSON
% 20210609
clear all
pkg load miscellaneous

% Enter some audio
[y, fs] = audioread("<your_fave_wave>");

ts = 1/fs;

% Band creation, N = number of bands, i_bw is the initial band BW
N = 8; 
i_bw = 40; 

% Lots of filtered random noise used
% simple FIR used here, n_noise is its order
fc_noise = 1;
wc_noise = fc_noise / fs; % argument W to fir1 is already normalized to 2pi
n_noise = 128;

% po_max is max phase offset
% mod_frac is fraction of bw for the band i.e. 1/12 = a semitone
po_max = 10 * ts;
mod_frac = 1 / 6;

% variable delay min/max in seconds
do_delay = 1;
d_min = 0.1e-3;
d_max = 8e-3;
%d_min = 40e-3;
%d_max = 60e-3;
d_max_samps = ceil(d_max * fs);

% filter the output before recombination
final_filter_n = 16;
do_final_filter = 1;

% Hilbert xformer params
% dp,ds pairing is roughly equivalent to the given d_phi 
wc = 0.495;
dp = 0.0001;
ds = 0.1;
d_phi = 0.012235;

%Pick your favourite hilbert transformer
[B, A, Bi, Ai] = iir_hilbert_even(wc, d_phi);
%[B, A, Bi, Ai] = iir_hilbert_odd(wc, dp, ds);

bw = i_bw;
t = [0:ts:ts * (length(y) -1)];
recon = zeros(length(y), 1);

for i = (1:N)
  % create bands
  fc = bw + bw / 2;
  [B, A] = iir_bp(fs, fc, bw);
  bp = filter(B, A, y);

  % each band gets its own modulation amount proportional to the band bw
  mod_max = bw * mod_frac; 

  % generate the filtered random noise, normalize it
  r = randn(length(y),1); 
  b = fir1(n_noise, wc_noise);
  r_filt = filter(b, 1, r);
  r_filt = r_filt / max(abs(r_filt));

  % Introduce some random phase offset
  tbar = t' .+ (po_max * r_filt);

  % spectrum shifter
  y0(i, :) = filter(B, A, bp)   .* sin(tbar * mod_max * 2 * pi);
  y1(i, :) = filter(Bi, Ai, bp) .* cos(tbar * mod_max * 2 * pi);


  if (do_delay)
    band_delay = zeros(length(y), 1);
    % generate some filtered random noise, normalize it
    % translate it to a range of [0, 1] so there's no -ve delay
    r = randn(length(y),1); 
    b = fir1(n_noise, wc_noise);
    r_filt = 0.5 + filter(b, 1, r);
    r_filt = r_filt / max(abs(r_filt));
    r_filt = d_min .+ (r_filt .* (d_max - d_min));

    % This takes some time so show some signs of life
    text_waitbar(i / N, "Delay Processing");

    % initialise to avoid -ve indexes
    band_delay(1:d_max_samps) = (y0(i, 1:d_max_samps) + y1(i, 1:d_max_samps));
    
    for k = (d_max_samps:length(y))
      del = floor(r_filt(k));
      band_delay(k) = (y0(i, k - del) + y1(i, k - del));
    endfor

    % player = audioplayer(band_delay / sqrt(2), fs);
    % playblocking(player);
    recon = recon .+ band_delay / sqrt(2);
  else
    recon = recon .+ ((y0(i, :) + y1(i, :)) / sqrt(2));
  endif

  bw = bw * 2;
endfor

if (do_final_filter)
  fc_recon = fs / 2;
  wc_recon = fc_recon / fs;
  b = fir1(final_filter_n, wc_recon);
  recon = filter(b, 1, recon);
endif

wd = 0.5;
mix = (wd .* y) + ((1 - wd) .* recon);
player = audioplayer(y, fs);
playblocking(player);
player = audioplayer(mix, fs);
play(player);
%plot(y, "-r;in;", mix, "-b;mix;")
%plot(ybar, "-r;in;", recon, "-b;recon;")
