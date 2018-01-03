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

/* Local Includes */
#include "ga_config.h"
#include "ga_helper.h"
#include "data.h"

///////////////////////////////////////////////////
/* CLASS:  WeightedSum */
///////////////////////////////////////////////////
class WeightedSum
{
public:
	void calculateFitness(StepPerformance_Vec input_data, PID_ControlGoals_sPtr input_goals, PIDFitness_Vec* output_fitness);

	WeightedSum(GA_RunMode execution_type);
	~WeightedSum();
private:
	GA_RunMode executionType;

	StepPerformance_Vec ws_data;
	PID_ControlGoals_sPtr ws_goals;
	PIDFitness_Vec* ws_fitness;
	boost::mutex ws_fitness_mutex;

	PID_Fitness calculateMemberFit(int memberNum, bool dataValid, double POS, double TS, double TR, double SSERR);

	void calculate_cpu_single_threaded();
	void calculate_cpu_multi_threaded();
	void calculate_gpu_single_threaded();
	void calculate_gpu_multi_threaded();
};

///////////////////////////////////////////////////
/* CLASS:  NonDominatedSort */
///////////////////////////////////////////////////
class NonDominatedSort
{
public:

private:
};

#endif