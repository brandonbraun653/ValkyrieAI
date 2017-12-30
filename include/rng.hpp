#pragma once
#ifndef RNG_HPP_H
#define RNG_HPP_H

/* C/C++ Standard Includes */
#include <stdint.h>

/* Boost Includes */
#include <boost/thread/mutex.hpp>
#include <boost/random.hpp>

template<typename Engine, typename Distribution>
class RNGManager
{
public:
	template<typename Scalar>
	Scalar get()
	{
		return (isLocked) ? ((Scalar)dist(eng)) : (Scalar)0;
	}

	bool acquireEngine()
	{
		return (isLocked = rng_mutex.try_lock());
	}
	
	void releaseEngine()
	{
		rng_mutex.unlock();
		isLocked = false;
	}

	RNGManager(Distribution& distribution)
	{
		dist = distribution;
		isLocked = false;

		/* Reseed and reset so the new object is unique */
		eng.seed(static_cast<std::uint32_t>(std::time(0)));
		dist.reset();

		/* Warm up the engine a bit */
		warmup();

		
	};
	
private:
	bool isLocked;
	boost::mutex rng_mutex;
	Engine eng;
	Distribution dist;

	void warmup()
	{
		for (int i = 0; i < 100; i++)
			dist(eng);
	}
};


#endif

