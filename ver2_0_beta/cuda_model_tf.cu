#include "cuda_model_tf.h"

/*------------------------------------------------------
* Kernel Forward Declarations 
*------------------------------------------------------*/
__global__ void square(float *x, float *y);

/*------------------------------------------------------
* STUB CALLS
*------------------------------------------------------*/
void TF_CPU_testSquare()
{
	const int N = 10000;
	float x[N], y[N];


	//Populate with some dummy data
	for (int i = 0; i < N; i++)
	{
		x[i] = 3.34;
		y[i] = 0.0;
	}
		

	//Figure out the size and then square away!
	for (int i = 0; i < N; i++)
		y[i] = x[i] * x[i];
}

void TF_GPU_testSquare(int blocks, int threads, float *d_x, float *d_y)
{
	square<<<blocks, threads>>>(d_x, d_y);
}


/*------------------------------------------------------
* GLOBAL KERNELS:
* Caller	-->	CPU
* Execute	-->	GPU
*------------------------------------------------------*/
__global__ void square(float *x, float *y)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	y[i] = x[i] * x[i];
}

/*------------------------------------------------------
* LOCAL KERNELS:
* Caller	-->	GPU
* Execute	-->	GPU
*------------------------------------------------------*/