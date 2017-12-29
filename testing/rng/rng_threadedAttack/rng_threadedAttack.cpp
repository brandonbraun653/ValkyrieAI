#include <stdint.h>
#include <stdlib.h>

#include "rng.hpp"


int main(void)
{

	RNGManager_UniformInt<boost::mt19937, 0, 15> testRNG;



	return 0;
}