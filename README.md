# MDT_GPU_F-engine

Source for the GPU Polyphase Filter (PPF) kernel in CUDA and the data used for
the test.

The data was acquired using Airspy R2 SDR. The files are in binary
format with samples arranged in FLOAT32 IQ, with file names RL4.25....dat. The
centre frequency is featured in the file names. Each data file contains 2^20 IQ
samples (1 MiS).

A code snippet showing how to extract the data is given in the file
sample_extract.c

Source code for the GPU kernels are found in file
AirSpy_FLOAT32IQ_GPUPPF.cu. It contains 2 kernels. The first one computes the
coefficients for the windowed-Sinc filter and the second kernel computes the
polyphase structures for a batch of spectra. Polyphase filtering is
completed with cuFFT kernels following the PPF kernel.

Finally, source code for the MDT PFB kernel is found in MDT_PFB_kernel.cu. It
contains only the PFB kernel which should be called after the coefficients
have been initialised using the Windowed-Sinc kernel. The PFB kernel can
process multiple data streams in parallel, i.e. the different signals coming
from each station.
