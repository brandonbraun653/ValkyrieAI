#pragma once
#ifndef GA_STEPS_SORTING_H_
#define GA_STEPS_SORTING_H_

/* C/C++ Includes */
#include <stdlib.h>
#include <stdint.h>
#include <random>
#include <math.h>

/* Boost Includes */
#include <boost/thread.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/container/vector.hpp>

/* Local Includes */
#include "config.h"
#include "ga_helper.h"
#include "types.h"

struct GA_SortingInput
{
	//boost::container::vector<double> parentChildFitScores;

	boost::container::vector<PID_FitnessScores> parentChildFitScores;
};

struct GA_SortingOutput
{
	boost::container::vector<int> sortedPopulation;
};

class GA_SortBase
{
public:
	virtual void sort(const GA_SortingInput, GA_SortingOutput&) = 0;
	virtual ~GA_SortBase() = default;
private:

};
typedef boost::shared_ptr<GA_SortBase> GA_SortBase_sPtr;

///////////////////////////////////////////////////
/* CLASS:  FastNonDominatedSort */
///////////////////////////////////////////////////
class FastNonDominatedSort : public GA_SortBase
{
public:
	void sort(const GA_SortingInput input, GA_SortingOutput& output) override;

	FastNonDominatedSort() = default;
	~FastNonDominatedSort() = default;
private:

	boost::container::vector<int> sortObjectiveFunc(boost::container::vector<int> memberSet,
		boost::container::vector<double> objFuncScores);

	boost::container::vector<int> sortCrowdingDistance(boost::container::vector<int> memberSet,
		boost::container::vector<double> crowdingDistances);
};


#endif