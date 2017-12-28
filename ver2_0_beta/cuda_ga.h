#ifndef CUDA_GA_H_
#define CUDA_GA_H_

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



/*------------------------------------------------------
* STUB CALLS
*------------------------------------------------------*/
extern void GPU_DisplayHeader();
extern void GPU_crossover(int blocks, int threads);
extern void GPU_generateMutationRates(int blocks, int threads);

#endif /* !CUDA_GA_H_ */