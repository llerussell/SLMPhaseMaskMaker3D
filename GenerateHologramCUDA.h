/*
    Header file for GenerateHologramCUDA
    Hologram generating algorithms for CUDA Devices
    Copyright 2009, 2010, 2011, 2012 Martin Persson
    martin.persson@physics.gu.se
    Small edits by Lloyd Russell 2016

    This file is part of GenerateHologramCUDA.
    GenerateHologramCUDA is free software: you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    GenerateHologramCUDA is distributed in the hope that it will be
    useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.
    You should have received a copy of the GNU Lesser General Public License
    along with GenerateHologramCUDA.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdint.h>

#ifndef GENERATEHOLOGRAMCUDA_H
#define GENERATEHOLOGRAMCUDA_H
    
#ifdef __cplusplus
extern "C" {  // only need to export C interface if
              // used by C++ source code
#endif


// Allocate GPU memory
__declspec( dllexport ) int startCUDA(float *h_pSLMstart, int deviceId);

// Generate a hologram
__declspec( dllexport ) int GenerateHologram(float *h_checkData, uint16_t *h_pSLM_uc, float *x_spots, float *y_spots, float *z_spots, float *I_spots, int N_spots, int N_iterations, float *h_Iobtained, int method);

// Free GPU memory
__declspec( dllexport ) int stopCUDA();


#ifdef __cplusplus
}
#endif

#endif