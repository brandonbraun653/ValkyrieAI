#pragma once
#ifndef EVAL_FITNESS_H_
#define EVAL_FITNESS_H_

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
#include "ga_config.h"
#include "ga_helper.h"
#include "types.h"

struct GA_EvaluateFitnessDataInput
{
	//Add whatever the heck is needed here
};

struct GA_EvaluateFitnessDataOutput
{
	//Add whatever the heck is needed here
};

class GA_EvaluateFitnessBase
{
public:
	virtual void evaluateFitness(const GA_EvaluateFitnessDataInput, GA_EvaluateFitnessDataOutput&) = 0;

private:
	//?
};
typedef boost::shared_ptr<GA_EvaluateFitnessBase> GA_EvaluateFitnessBase_sPtr;


///////////////////////////////////////////////////
/* CLASS:  WeightedSum */
///////////////////////////////////////////////////
class WeightedSum : public GA_EvaluateFitnessBase
{
public:
	void evaluateFitness(const GA_EvaluateFitnessDataInput input, GA_EvaluateFitnessDataOutput& output) override;

	WeightedSum();
	~WeightedSum();
private:
};

///////////////////////////////////////////////////
/* CLASS:  NonDominatedSort */
///////////////////////////////////////////////////
class NonDominatedSort : public GA_EvaluateFitnessBase
{
public:
	void evaluateFitness(const GA_EvaluateFitnessDataInput input, GA_EvaluateFitnessDataOutput& output) override;

	NonDominatedSort();
	~NonDominatedSort();
private:
};

#endif