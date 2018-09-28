/*
 *  sample_extract.c
 *
 *  Code snippet showing how extract digital samples from radio signal data
 *  obtained from Airspy R2 SDR stored in a binary file.
 *
 *  Copyright (C) 2018 Nitish Ragoomundun, Mauritius
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

#include <stdio.h>


  /***  BEGIN Read raw data from input file  ***/
  /*
   *  Data in order: I1 Q1 I2 Q2 I3 ...
   *  Each I and Q is a 32-bit (4B) integer sample.
   *  If we want N samples of whole data, we need to
   *  read N (I,Q) pairs, therefore 2N objects should
   *  be read from binary file.
   *
   *  f1 : file pointer,
   *  FILENAME : name of file with data,
   *  RawIQ : array to hold raw samples,
   *  Nchannels : number of channels in spectra,
   *  Ntaps : number of taps implemented in PPF,
   *  retval : integer return value.
   *
   */
  printf("Reading raw binary data into array ...\n");
  i = fread((void *)RawIQ, 4, 2*Nchannels*Ntaps, f1);
  if (i < 2*Nchannels*NTAPS)
  {
    printf("Could not read %ld data samples from file %s!\n", 2*Nchannels*Ntaps, FILENAME);
    retval = ferror(f1);
    fclose(f1);
    exit(retval);
  }

  fclose(f1);

  /***  END Read raw data from input file  ***/
