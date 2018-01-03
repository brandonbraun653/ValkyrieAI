#include <stdint.h>
#include <stdlib.h>

#include "rng.hpp"
#include <boost/bind.hpp>
#include <boost/thread.hpp>

using namespace boost::random;

void rng_attack_thread();
uniform_real_distribution<> test(0.0, 1000.0);
RNGManager<mt19937, uniform_real_distribution<>> rng2(test);

int main(void)
{
	boost::thread_group threadGroup;
	
	auto start = boost::chrono::high_resolution_clock::now();
	for (int i = 0; i < 8; i++)
	{
		threadGroup.create_thread(rng_attack_thread);
	}
	threadGroup.join_all();
	auto end = boost::chrono::high_resolution_clock::now();
	auto dt = boost::chrono::duration_cast<boost::chrono::milliseconds>(end - start);

	std::cout << "Done in " << dt.count() << " mS" << std::endl;
	return 0;
}

void rng_attack_thread()
{
	while (rng2.acquireEngine() == false)
	{
		boost::this_thread::yield();
	}

	/* Hurray! We achieved a lock. Grab some delicious numbers. */
	std::cout << "Thread: " << boost::this_thread::get_id() << " got a lock!" << std::endl;
	
	volatile int result = 0;
	for (int i=0; i<100; i++)
	{
		result = rng2.get<int>();
	}
	

	rng2.releaseEngine();
}