#pragma once
#ifndef FILTER_POPULATION_H_
#define FILTER_POPULATION_H_

/* C/C++ Includes */
#include <stdlib.h>
#include <stdint.h>
#include <iostream>
#include <random>
#include <math.h>

/* Boost Includes */
#include <boost/thread.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/random.hpp>

/* Local Includes */
#include "config.h"
#include "ga_helper.h"
#include "types.h"
#include "rng.hpp"

struct GA_PopulationFilterDataInput
{
	boost::container::vector<double> currentGlobalFitScores;		/* Latest fitness scores */

	double static_performanceThreshold = 0.0;		/* All members who perform below this threshold will be replaced. Range (0.0,1.0)
													Only used in the static filter. Leave empty if unused. */

	

	int forcedReplacementQuota = 0;					/* TODO: For a later filtering method */
	double naturalDisasterPercentage = 0.0;			/* TODO: For a later filtering method */
};

struct GA_PopulationFilterDataOutput
{
	boost::container::vector<int> replacedMemberIndexes;							/* Index of the population member that is to be killed */
	boost::container::vector<GA_PIDChromosome<double>> replacementPIDValues;		/* For a given member, this is the new genetic material */
};

class GA_PopulationFilterBase
{
public:
	virtual void filter(const GA_PopulationFilterDataInput, GA_PopulationFilterDataOutput&) = 0;

	virtual ~GA_PopulationFilterBase() = default;
private:

};
typedef boost::shared_ptr<GA_PopulationFilterBase> GA_PopulationFilterBase_sPtr;


///////////////////////////////////////////////////
/* CLASS:  StaticFilter */
///////////////////////////////////////////////////
class StaticFilter : public GA_PopulationFilterBase
{
public:
	void filter(const GA_PopulationFilterDataInput input, GA_PopulationFilterDataOutput& output) override;

	StaticFilter(const double KpMax, const double KpMin, const double KiMax, const double KiMin, const double KdMax, const double KdMin);
	~StaticFilter();

private:
	RNGManager_sPtr rngKp, rngKi, rngKd;

};


///////////////////////////////////////////////////
/* CLASS:  DynamicFilter */
///////////////////////////////////////////////////
class DynamicFilter : public GA_PopulationFilterBase
{
public:
	void filter(const GA_PopulationFilterDataInput input, GA_PopulationFilterDataOutput& output) override;

	DynamicFilter();
	~DynamicFilter();

private:
};

#endif