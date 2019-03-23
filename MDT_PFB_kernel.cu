/*
 *  Kernel to compute the polyphase structure for the MDT. Computation of
 *  signal from several different stations is done in parallel.
 *
 *  Copyright (C) 2019 Nitish Ragoomundun
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


/*
 *  GPU kernel to compute polyphase structure.
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
__global__ void PPFBatch(int Ntaps, long Nchannels, float *Filter, cufftComplex *InSignal, cufftComplex *PolyStruct)
{
  int i;
  int channelIdx = threadIdx.x + blockIdx.x*blockDim.x;
  long stride_elem, stride_spec;
  float tmp_filter;
  cufftComplex tmp_input, tmp_product;

  /*  Stride wrt element index  */
  stride_elem = blockIdx.z * (gridDim.y - 1 + Ntaps) * Nchannels;

  /*  Stride wrt number of previous spectra  */
  stride_spec = blockIdx.y * Nchannels;

  tmp_product.x = 0.0f;
  tmp_product.y = 0.0f;

  for ( i=0 ; i<Ntaps ; i++ )
  {
    /*  Read input signal data and filter coefficient  */
    tmp_input  = InSignal[stride_elem + stride_spec + i*Nchannels + channelIdx];
    tmp_filter = Filter[i*Nchannels + channelIdx];

    /*  Accumulate FIR  */
    tmp_product.x = fmaf(tmp_filter, tmp_input.x, tmp_product.x);  // I
    tmp_product.y = fmaf(tmp_filter, tmp_input.y, tmp_product.y);  // Q
  }

  /*  Write output to array in global memory  */
  PolyStruct[blockIdx.z*gridDim.y*Nchannels + stride_spec + channelIdx] = tmp_product;
}
