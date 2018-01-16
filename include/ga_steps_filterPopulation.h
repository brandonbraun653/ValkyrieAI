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

/* Local Includes */
#include "config.h"
#include "ga_helper.h"
#include "types.h"

struct GA_PopulationFilterDataInput
{
	//Add whatever the heck is needed here
};

struct GA_PopulationFilterDataOutput
{
	//Add whatever the heck is needed here
};

class GA_PopulationFilterBase
{
public:
	virtual void filter(const GA_PopulationFilterDataInput, GA_PopulationFilterDataOutput&) = 0;

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

	StaticFilter();
	~StaticFilter();

private:
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