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
#include "config.h"
#include "ga_helper.h"
#include "types.h"

struct GA_EvaluateFitnessDataInput
{
	PID_PerformanceGoals goals;				/* User defined reference goals for the optimizer to achieve */
	PID_PerformanceTolerance tolerance;		/* User defined acceptable tolerance bounds for results */

	double POS;								/* Percent Overshoot */
	double SSER;							/* Steady State Error */
	double TR;								/* Time Rise */
	double TS;								/* Time Settle */
	double FV;								/* Final Value */

	Eigen::MatrixXd simulationData;			/* Raw simulation data */
};

struct GA_EvaluateFitnessDataOutput
{
	PID_FitnessScores fit;	/* Calculated fitness metrics for the given input data */
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