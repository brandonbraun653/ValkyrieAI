#pragma once
#ifndef SELECT_PARENTS_H_
#define SELECT_PARENTS_H_

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
/* CLASS:  RankedSelection */
///////////////////////////////////////////////////
class RankedSelection
{
public:
	GA_RunMode executionType;

	RankedSelection(GA_RunMode execution_type);
	~RankedSelection();
private:
};

///////////////////////////////////////////////////
/* CLASS:  RandomSelection */
///////////////////////////////////////////////////
class RandomSelection
{
public:
	void selectParents(iVec* out_selections);

	RandomSelection(GA_RunMode execution_type);
	~RandomSelection();
private:
	GA_RunMode executionType;
	iVec* rs_selections;

	void calculate_cpu_single_threaded();
	void calculate_gpu_single_threaded();
};

///////////////////////////////////////////////////
/* CLASS:  RouletteSelection */
///////////////////////////////////////////////////
class RouletteSelection
{
public:

private:
};

///////////////////////////////////////////////////
/* CLASS:  TournamentSelection */
///////////////////////////////////////////////////
class TournamentSelection
{
public:
	void selectParents(PIDFitness_Vec in_fitnessValues, iVec* out_selections);

	TournamentSelection(GA_RunMode execution_type);
	~TournamentSelection();

private:
	GA_RunMode executionType;
	PIDFitness_Vec ts_data;
	iVec* ts_selections;

	void calculate_cpu_single_threaded();
	void calculate_cpu_multi_threaded();
	void calculate_gpu_single_threaded();
	void calculate_gpu_multi_threaded();
};

///////////////////////////////////////////////////
/* CLASS:  StochasticSamplingSelection */
///////////////////////////////////////////////////
class StochasticSamplingSelection
{
public:

private:
};

///////////////////////////////////////////////////
/* CLASS:  Elitist */
///////////////////////////////////////////////////

#endif