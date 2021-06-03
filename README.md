# iir_hilbert
IIR Hilbert Transformers

This is where a quest to undestand the "Phase Difference Network" (PDN) shown in Fig 14-25 of Hal Chamberlin's Musical Applications of Microprocessors 2nd Edition has led me.
The PDN makes up part of a proposed Parametric Chorus Simulator.
I still don't fully understand how the values used in his version were defined, but it set me to trying to understand how it worked.

Michael Gerzon also wrote about PDNs in some of his work on Synthetic Stero Reverb.

Research led me to using the Hilbert Transform.
Octave has hilbert.m - but that's FFT based.
Csound has a hilbert transform opcode that is implemented as an IIR Hilbert transformer.
A particularly intriguing comment states:

    /* pole values taken from Bernie Hutchins, "Musical Engineer's Handbook" */

The magic pole numbers are from the MEH which can not be obtained for love nor money these days. 

So I set about following the crumbs, reading papers and scouring old second hand DSP textbooks.
This is what I've ended up with, a couple of functions that can generate the polynomials needed for IIR implementation of PDNs.

These .m functions were developed using GNU Octave and its signal package.
I've tried to provide references for all the steps.
