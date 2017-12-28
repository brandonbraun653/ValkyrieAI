#ifndef CUDA_MODEL_TF_H_
#define CUDA_MODEL_TF_H_

/*------------------------------------------------------
* INCLUDES
*------------------------------------------------------*/
/* Cuda Dependencies */
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda_memory.h"

/* C/C++ Dependencies */
#include <stdlib.h>
#include <iostream>

/* Testing Functions */
extern void TF_CPU_testSquare();
extern void TF_GPU_testSquare(int blocks, int threads, float *d_x, float *d_y);

/*------------------------------------------------------
* STUB CALLS
* Note: These calls are only meant to be used with the 
* Transfer Function class model
*------------------------------------------------------*/
extern void TF_GPU_evaluateModel(int blocks, int threads);





#endif /* !CUDA_MODEL_TF_H_*/