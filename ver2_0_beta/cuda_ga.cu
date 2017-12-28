#include "cuda_ga.h"


//Eventually I'll add some more things here.
void GPU_DisplayHeader()
{
	const int kb = 1024;
	const int mb = kb * kb;
	std::cout << "NBody.GPU" << std::endl << "=========" << std::endl << std::endl;

	std::wcout << "CUDA Version:   v" << CUDART_VERSION << std::endl;
	//std::wcout << "Thrust version: v" << THRUST_MAJOR_VERSION << "." << THRUST_MINOR_VERSION << std::endl << std::endl;

	int devCount;
	cudaGetDeviceCount(&devCount);
	std::wcout << "CUDA Devices: " << devCount << std::endl << std::endl;

	for (int i = 0; i < devCount; ++i)
	{
		cudaDeviceProp props;
		cudaGetDeviceProperties(&props, i);
		std::wcout << i << ": " << props.name << ": " << props.major << "." << props.minor << std::endl;
		std::wcout << "  Global memory:   " << props.totalGlobalMem / mb << "mb" << std::endl;
		std::wcout << "  Shared memory:   " << props.sharedMemPerBlock / kb << "kb" << std::endl;
		std::wcout << "  Constant memory: " << props.totalConstMem / kb << "kb" << std::endl;
		std::wcout << "  Block registers: " << props.regsPerBlock << std::endl << std::endl;

		std::wcout << "  Warp size:         " << props.warpSize << std::endl;
		std::wcout << "  Threads per block: " << props.maxThreadsPerBlock << std::endl;
		std::wcout << "  Max block dimensions: [ " << props.maxThreadsDim[0] << ", " << props.maxThreadsDim[1] << ", " << props.maxThreadsDim[2] << " ]" << std::endl;
		std::wcout << "  Max grid dimensions:  [ " << props.maxGridSize[0] << ", " << props.maxGridSize[1] << ", " << props.maxGridSize[2] << " ]" << std::endl;
		std::wcout << std::endl;
	}
}
