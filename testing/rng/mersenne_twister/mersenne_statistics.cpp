#include <stdlib.h>
#include <stdint.h>
#include <iostream>
#include <fstream>
#include <string.h>

#include <boost/container/vector.hpp>
#include <boost/random.hpp>

#define RNG_LIMIT 500000

std::string filename_mt19937_continuous = "C:/git/Github/ValkyrieAI/testing/rng/mersenne_twister/mt19937_continuous.csv";
std::string filename_mt19937_reinit = "C:/git/Github/ValkyrieAI/testing/rng/mersenne_twister/mt19937_reinit.csv";
std::string filename_mt19937_continuous_reseed = "C:/git/Github/ValkyrieAI/testing/rng/mersenne_twister/mt19937_continuous_reseed.csv";
std::string filename_mt19937_reinit_reseed = "C:/git/Github/ValkyrieAI/testing/rng/mersenne_twister/mt19937_reinit_reseed.csv";
std::ofstream myFile;

int main(void)
{
	std::cout << "Starting batch test..." << std::endl;
	
	/*------------------------------------------
	 * Test RNG generation of single mt19937 over a range of -100.0 <-> 100.0 for 
	 * accurate uniform distribution. This is not generally how the RNG is 
	 * used in the previous generation of AI software. 
	 *-----------------------------------------*/
	{ //Force scoping of local vars
		std::cout << "Beginning continuous mt19937 test:" << std::endl;
		myFile.open(filename_mt19937_continuous);
		boost::mt19937 rng;
		boost::uniform_real<double> rng_range(-100.0, 100.0);
		boost::variate_generator<boost::mt19937, boost::uniform_real<double>>
			rng_engine(rng, rng_range);

		unsigned int printThreshold = 25000;
		for (unsigned int i = 0; i < RNG_LIMIT; i++)
		{
			myFile << rng_engine() << ",";
			if (i == printThreshold)
			{
				std::cout << i << std::endl;
				printThreshold += 25000;
			}

		}
		myFile.close();
		std::cout << "Finished!" << std::endl;
	}

	/*------------------------------------------
	* Test RNG generation of single mt19937 over a range of -100.0 <-> 100.0 for
	* accurate uniform distribution trying out a reseeded/reset version. I don't 
	* expect much to happen here.
	*-----------------------------------------*/
	{
		std::cout << "Beginning continuous_reseed mt19937 test:" << std::endl;
		myFile.open(filename_mt19937_continuous_reseed);
		boost::mt19937 rng;
		boost::uniform_real<double> rng_range(-100.0, 100.0);
		boost::variate_generator<boost::mt19937, boost::uniform_real<double>>
			rng_engine(rng, rng_range);

		rng_engine.engine().seed();				/* Randomize the seed */
		rng_engine.distribution().reset();		/* Reset so the next iteration can't depend on prev */

		unsigned int printThreshold = 25000;
		for (unsigned int i = 0; i < RNG_LIMIT; i++)
		{
			myFile << rng_engine() << ",";
			if (i == printThreshold)
			{
				std::cout << i << std::endl;
				printThreshold += 25000;
			}

		}
		myFile.close();
		std::cout << "Finished!" << std::endl;
	}


	/*------------------------------------------
	* Test RNG generation of multiple mt19937 over a range of -100.0 <-> 100.0 for
	* accurate uniform distribution. This is more like how it was done in the 
	* previous generation of AI software. Generally an engine would have only a 
	* few opportunities for generation before being destroyed and reinitialized.
	*-----------------------------------------*/
	{
		std::cout << "Beginning continuous mt19937 test:" << std::endl;
		myFile.open(filename_mt19937_reinit);

		unsigned int printThreshold = 25000;
		unsigned int destroyThreshold = 10;
		for (unsigned int i = 0; i < floor(RNG_LIMIT/destroyThreshold); i++)
		{
			/* Recreate the RNG engine after a few iterations */
			boost::mt19937 rng;
			boost::uniform_real<double> rng_range(-100.0, 100.0);
			boost::variate_generator<boost::mt19937, boost::uniform_real<double>>
				rng_engine(rng, rng_range);

			/* Generate a small number of numbers */
			for (unsigned int j = 0; j < destroyThreshold; j++)
			{
				myFile << rng_engine() << ",";
			}

			/* All engines destroyed here */
		}
		myFile.close();
		std::cout << "Finished!" << std::endl;
	}

	/*------------------------------------------
	* Test RNG generation of multiple mt19937 over a range of -100.0 <-> 100.0 for
	* accurate uniform distribution. See if reset/reseed improves operation.
	*-----------------------------------------*/
	{
		std::cout << "Beginning continuous mt19937 test:" << std::endl;
		myFile.open(filename_mt19937_reinit_reseed);

		unsigned int printThreshold = 25000;
		unsigned int destroyThreshold = 10;
		for (unsigned int i = 0; i < floor(RNG_LIMIT / destroyThreshold); i++)
		{
			/* Recreate the RNG engine after a few iterations */
			boost::mt19937 rng;
			boost::uniform_real<double> rng_range(-100.0, 100.0);
			boost::variate_generator<boost::mt19937, boost::uniform_real<double>>
				rng_engine(rng, rng_range);

			rng_engine.engine().seed();				/* Randomize the seed */
			rng_engine.distribution().reset();		/* Reset so the next iteration can't depend on prev */

			/* Generate a small number of numbers */
			for (unsigned int j = 0; j < destroyThreshold; j++)
			{
				myFile << rng_engine() << ",";
			}

			/* All engines destroyed here */
		}
		myFile.close();
		std::cout << "Finished!" << std::endl;
	}
}