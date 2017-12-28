#ifndef CUDA_MEMORY_H
#define CUDA_MEMORY_H

/*------------------------------------------------------
* INCLUDES
*------------------------------------------------------*/
/* Cuda Dependencies */
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

/* C/C++ Dependencies */
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

/*------------------------------------------------------
* TYPES
*------------------------------------------------------*/
enum memcpy_direction
{
	HOST_TO_DEVICE,
	DEVICE_TO_HOST,
	DEVICE_TO_DEVICE,
	DEFAULT
};

/*------------------------------------------------------
* FUNCTIONS
*------------------------------------------------------*/
/* Allocates memory on the GPU */
template <typename data_type>
data_type* cuda_malloc(size_t mem_size)
{
	data_type *ptr;
	if (cudaMalloc(&ptr, mem_size * sizeof(data_type)) != cudaSuccess)
		std::cout << "CUDA MEMORY ALLOCATION ERROR!!" << std::endl;
	return ptr;
}

/* Frees up memory allocated on the GPU */
template <typename data_type>
cudaError_t cuda_free(data_type *cuda_mem)
{
	cudaError_t error = cudaFree(cuda_mem);
	if (error != cudaSuccess)
		std::cout << "cudaFree error: " << error << std::endl;

	return error;
}

/* Copies data between devices, specified by memcpy_direction */
template <typename data_type>
cudaError_t cuda_memcpy(data_type* dst, data_type* src, size_t count, memcpy_direction direction)
{
	cudaMemcpyKind copyKind;
	switch (direction)
	{
	case HOST_TO_DEVICE:
		copyKind = cudaMemcpyHostToDevice;
		break;

	case DEVICE_TO_HOST:
		copyKind = cudaMemcpyDeviceToHost;
		break;

	case DEVICE_TO_DEVICE:
		copyKind = cudaMemcpyDeviceToDevice;
		break;

	case DEFAULT:
		copyKind = cudaMemcpyDefault;
		break;

	default:
		copyKind = cudaMemcpyDefault;
		break;
	}

	cudaError_t error = cudaMemcpy(dst, src, count * sizeof(data_type), copyKind);

	if (error != cudaSuccess)
	{
		char *outputCode;

		switch (error)
		{
		case cudaErrorInvalidValue:
			outputCode = "Invalid Value";
			break;

		default:
			outputCode = "Unknown Error";
			break;
		}

		std::cout << "cudaMemcpy error: " << outputCode << " - " << error << std::endl;
	}

	return error;
}

#endif /* !CUDA_MEMORY_H_ */
