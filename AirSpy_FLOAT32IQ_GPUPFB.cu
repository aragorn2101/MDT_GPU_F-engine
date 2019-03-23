/*
 *  GPU kernels to compute the coefficients of the windowed-Sinc filter and
 *  kernel for computing the polyphase structure. The PPF is completed with
 *  an FFT kernel by calling the CUDA library cuFFT.
 *
 *
 *  Copyright (C) 2018 Nitish Ragoomundun
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */


#include <cuda.h>

#define Fs 1e7  // Sampling frequency (MSps) for AirSpy in mode FLOAT32 IQ
#define NTAPS 8  // Number of taps for computing PPF
#define MAXTHREADPERBLOCK 1024  // Maximum number of threads per block for GTX Titan X

/*  Constants for GPU  */
__device__ __constant__ float d_PI = 3.141592653589793238462643383279502884f;



/*
 *  Kernel to compute FIR Filter.
 *
 *  Grid dimensions:
 *  NumThreadx = MAXTHREADPERBLOCK,
 *  NumThready = 1,
 *  NumThreadz = 1,
 *  NumBlockx = (Ntaps * Nchannels) / NumThreadx,
 *  NumBlocky = 1,
 *  NumBlockz = 1
 *
 *  Ntaps: number of taps,
 *  Nchannels: number of channels in output,
 *  Filter: array of size Ntaps x Nchannels which will hold coefficients
 *          of FIR filter.
 */
__global__ void WindowedSinc(int Ntaps, long Nchannels, float *Filter)
{
  long idx = threadIdx.x + blockIdx.x*blockDim.x;

  /*  Temporary variables to prevent redundant computation  */
  float tmp1, tmp2, tmp_filter, tmp_window;

  if (idx < Ntaps*Nchannels)
  {

    /*
     *  Filter: Sinc
     *
     *  sinc( ( channel - (Ntaps x Nchannels)/2 ) / Nchannels )
     *
     */
    tmp1 = (idx - 0.5f*Ntaps*Nchannels) / Nchannels;

    if ( tmp1 == 0.0f )  /*  To prevent division by 0  */
      tmp_filter = 1.0f ;
    else
      tmp_filter = sinpif( tmp1 ) / ( d_PI * tmp1 );


    /*
     *  Window: Exact Blackman
     *
     *  a0 - a1cos( 2 x PI x i / (Ntaps x Nchannels) ) + a2cos( 4 x PI x i / (Ntaps x Nchannels) )
     *
     */
    tmp2 = 2.0f*idx / (Ntaps*Nchannels);
    tmp_window = 0.42659071f - 0.49656062f*cospif(tmp2) + 0.07684867f*cospif(2.0f*tmp2);


    /*  Write Windowed Sinc to global memory array  */
    Filter[idx] = tmp_filter * tmp_window;
  }
}



/*
 *  Kernel to compute polyphase structure.
 *
 *  Grid dimensions:
 *  NumThreadx = Number of channels computed per block,
 *  NumThready = 1,
 *  NumThreadz = 1,
 *  NumBlockx = Nchannels / NumThreadx,
 *  NumBlocky = Nspectra,
 *  NumBlockz = Nelements
 *
 *  Ntaps: number of taps,
 *  Nchannels: number of channels in output,
 *  Filter: array of size Ntaps x Nchannels containing FIR filter,
 *  InSignal: array of size 2 x Ntaps x Nchannels containing
 *            interleaved IQ samples of the input signal in an
 *            an array of float2 vectors (I:x, Q:y),
 *  PolyStruct: array of size Nchannels which will hold output.
 *
 */
__global__ void PPFBatch(int Ntaps, long Nchannels, float *Filter, float2 *InSignal, cufftComplex *PolyStruct)
{
  int i;
  int channelIdx = threadIdx.x + blockIdx.x*blockDim.x;
  long stride_element, stride_spect, stride_taps;
  float2 tmp_input;
  float tmp_filter;
  cufftComplex tmp_product;

  /*  Stride wrt element index  */
  stride_element = blockIdx.z * (gridDim.y - 1 + Ntaps) * Nchannels;

  /*  Stride wrt number of previous spectra  */
  stride_spect = blockIdx.y * Nchannels;

  tmp_product.x = 0.0f;
  tmp_product.y = 0.0f;

  for ( i=0 ; i<Ntaps ; i++ )
  {
    /*  Stride in spectrum wrt previous number of taps  */
    stride_taps = i*Nchannels;

    /*  Read input signal data and filter coefficient  */
    tmp_input = InSignal[stride_element + stride_spect + stride_taps + channelIdx];
    tmp_filter = Filter[stride_taps + channelIdx];

    /*  Accumulate FIR  */
    tmp_product.x = fmaf(tmp_filter, tmp_input.x, tmp_product.x);  // I
    tmp_product.y = fmaf(tmp_filter, tmp_input.y, tmp_product.y);  // Q
  }

  /*  Write output to array in global memory  */
  PolyStruct[(blockIdx.z * gridDim.y * Nchannels) + stride_spect + channelIdx] = tmp_product;
}
