# MDT_GPU_F-engine

Source for the GPU Polyphase Filter (PPF) kernel in CUDA and the data used for the test.
The data was acquired using Airspy R2 SDR and the files are in binary format with samples arranged in FLOAT32 IQ.
Each data file contains 2^20 IQ samples (1 MiS).

A code snippet to extract the data is given in the file sample_extract.c
