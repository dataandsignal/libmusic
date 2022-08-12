# Welcome to libmusic
### a frequency detection library

![2tone_psd](https://user-images.githubusercontent.com/40000574/190140071-8672a878-4146-49a2-8a79-20bb2e777f86.jpg)



## Overview

libmusic is a frequency detection library, that uses [MUSIC](https://en.wikipedia.org/wiki/MUSIC_(algorithm)) algorithm. It belongs to parametric methods of spectral estimation which working principle is to decompose signal space into pure signal and noise subspaces, using SVD of autocorrelation matrix. 

This allows for superresolution detection, where [Heisenberg uncertainty principle](https://en.wikipedia.org/wiki/Uncertainty_principle) no longer holds and one is able to compute pseudospectrum at every point up to half sampling rate just from a (very) limited amount of samples, and with better accuracy than parametric methods of spectrum estimation do (e.g. periodogram).

For applications with real time requirements (playout buffers, jitter, etc.) this provides significantly reduced latency over other detection methods and opens up possibilities to detect from minimum amount of data. 

In telephony for example, by using this method instead of popular Goertzel algorithm, 
one can detect DTMF tones with samples requirement reduced from 110 to 8 and frequency resolution (accuracy) increased from 72.73 Hz to 10^-2 Hz in the same time.


## Performance in noise

In a noisy environment, MUSIC performs well as long as SNR is above 66.53 dB.


## Applications

In telephony, DTMF symbols must be removed from stream, due too sensitive data protection. Often though, fractions of those DTMFs are left in a stream and must be removed. This cannot be done with Goertzel algorithm as it needs more than 110 samples to achieve right frequency resolution, to be able to detect DTMF frequencies. An example of such DTMF fraction is shown on the picture. This one is 14 samples in length (1.75 ms at a sampling rate of 8000 Hz).

![dtmf_test_vector_valid](https://user-images.githubusercontent.com/40000574/190151206-2e7b78a0-0d79-459f-bf8f-cf422fd9da72.jpg)

With MUSIC, samples requirement is reduced from 110 to 8 and frequency resolution (accuracy) increased from 72.73 Hz to 10^-2 Hz in the same time.
This picture presents correctness as a percentage of detected fractions of dual tone signals (DTMFs), by input vector length **N** (8,9,10,11,12,14), autocorrelation order **M** (4-8) and fraction length **L** (8-28 samples).

![dtmf_test_valid_freq_2](https://user-images.githubusercontent.com/40000574/190215259-e8a2c69e-921d-4c7b-99a7-69d4fb2ece7a.jpg)



For example, using a block of N=12 samples, all fractions of length L=10 and above can be detected (with autocorrelation order M={6,7}). N=8 detects all fractions longer than 8 samples (1 ms) with M=4.


## Discussion

If you have any questions, or would like to share your thoughts, 
request new features, etc. - please post them to [Discussions](https://github.com/dataandsignal/libmusic/discussions).


## Issues

If you would like to report an issue, please do it on [issues](https://github.com/dataandsignal/libmusic/issues) page.


## MATLAB

A MATLAB sandbox for libmusic: [libmusic_m](https://github.com/dataandsignal/libmusic_m)


## Repository 

URL: https://github.com/dataandsignal/libmusic


## Further reading, references

A good references about spectral analysis and space decomposition methods are:

- Hayes M. H., Statistical Digital Signal Processing And Modeling, Georgia Institute of Technology, Wiley, 2008
- Lawrence Marple S. Jr., Digital Spectral Analysis, Dover Publications, 2019
- Schmidt R. O., Multiple Emitter Location and Signal Parameter Estimation, IEEE Transactions on Antennas and Propagation, Vol. AP-34, No. 3, 1986

These references are missing though (or skipping intentionally) a crucial result about autocorrelation and sinusoids embedded in a vector space whose elements are shifted samples of that same sinusoid (with all the phases). This is a fundamental finding space decomposition methods are built on.

This is explained in more detail in:

- Penny W. D., Signal Processing Course, University College London, 2000

- and also in my [engineering thesis](https://drive.google.com/file/d/1e2LjLYKVGIdSypj2sSbbauz2-SU5a1as/view?usp=sharing) (written in polish, will probably be translated to english) 




## Copyright 

LIBMUSIC

Copyright (C) 2018-2022, Piotr Gregor piotr@dataandsignal.com

August 2022
