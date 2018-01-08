#pragma once
#ifndef RNG_HPP_H
#define RNG_HPP_H

/* C/C++ Standard Includes */
#include <stdint.h>

/* Boost Includes */
#include <boost/thread/mutex.hpp>
#include <boost/random.hpp>
#include <boost/shared_ptr.hpp>

class RNGManager
{
public:
	virtual int getInt() = 0;
	virtual double getDouble() = 0;

	virtual bool acquireEngine() = 0;
	virtual void releaseEngine() = 0;

	
private:
	virtual void warmup() = 0;
};
typedef boost::shared_ptr<RNGManager> RNGManager_sPtr;


template<typename Engine, typename Distribution>
class RNGInstance : RNGManager
{
public:
	double getInt() override
	{
		return (isLocked) ? ((int)dist(eng)) : (int)0;
	}

	double getDouble() override
	{
		return (isLocked) ? ((double)dist(eng)) : (double)0;
	}

	bool acquireEngine() override
	{
		return (isLocked = rng_mutex.try_lock());
	}
	
	void releaseEngine() override
	{
		rng_mutex.unlock();
		isLocked = false;
	}

	RNGInstance(Distribution& distribution)
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

	void warmup() override
	{
		for (int i = 0; i < 100; i++)
			dist(eng);
	}
};


#endif

