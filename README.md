Code for frequency domain analysis of cortical population recordings.

The code exemplifies the main analyses in our paper entitled "Distinct structure of cortical population activity on fast and infraslow timescales"
(by M. Okun, N. Steinmetz, A. Lak, M. Dervinis and K. Harris, bioRxiv link: https://doi.org/10.1101/395251 ). 
It demonstrates estimation of power spectrum of spike trains, and of (rate adjusted) coherence and phase between spike trains of individual neurons 
and the population rate.

Data of an example recording (populationSpikeData.zip) is also provided. The main script is in the exampleAnalysis.m file.

The code relies on the following two packages:
* Circular Statistics Toolbox ( https://github.com/circstat/circstat-matlab )
* Chronux (2.12) to which several corrections and GPU support were added ( https://github.com/m-okun/chronuxGPU )
  Note that GPU can be used for speed gain, but is not required.
