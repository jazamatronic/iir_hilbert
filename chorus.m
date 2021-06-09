%
% Parametric chorus simulator
% From Fig 14-25 of Hal Chamberlin's Musical Applications of Microprocessors - Second Edition
%
% Jared ANDERSON
% 20210609
clear all

% Enter some audio
[y, fs] = audioread("<your_fave_wave>");

ts = 1/fs;

% Band creation, N = number of bands, i_bw is the initial band BW
N = 8; 
i_bw = 40; 

% Lots of filtered random noise used
% simple FIR used here, n_noise is its order
fc_noise = 5;
wc_noise = fc_noise / fs; % argument W to fir1 is already normalized to 2pi
n_noise = 128;

% po_max is max phase offset
% mod_frac is fraction of bw for the band i.e. 1/12 = a semitone
po_max = 10 * ts;
mod_frac = 1 / 12;

% variable delay min/max in seconds, d_per is delay update period
do_delay = 1;
d_min = 0.1e-3;
d_max = 8e-3;
d_per = d_max * 8;

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
band_delay = zeros(1, length(y) + ceil(d_max / ts));
recon = zeros(1, length(y) + ceil(d_max / ts));

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
    % need to find a better way to do this.
    % generate the filtered random noise, normalize it
    % translate it to a range of [0, 1] so there's no -ve delay
    d_block = ceil(d_per / ts);
    r = randn(ceil(length(y) / d_block),1); 
    b = fir1(n_noise, wc_noise);
    r_filt = 0.5 + filter(b, 1, r);
    r_filt = r_filt / max(abs(r_filt));
    r_filt = d_min .+ (r_filt .* (d_max - d_min));
    for k = (1:(length(y) / d_block))
      del = floor(r_filt(k) / ts);
      for h = (1:d_block)
        idx = (k - 1) * d_block + h; 
        band_delay(idx + del) = (y0(i, idx) + y1(i, idx)) / sqrt(2);
      endfor
    endfor

    recon = recon .+ band_delay;
  else
    recon = recon .+ [((y0(i, :) + y1(i, :)) / sqrt(2)) zeros(1, ceil(d_max / ts))];
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
ybar = [y' zeros(1, ceil(d_max / ts))];
mix = (wd .* ybar) + ((1 - wd) .* recon);
player = audioplayer(y, fs);
playblocking(player);
player = audioplayer(mix, fs);
play(player);
%plot(y, "-r;in;", mix, "-b;mix;")
