#pragma once
#ifndef EVALUATE_MODEL_H_
#define EVALUATE_MODEL_H_

/* C/C++ Includes */
#include <stdlib.h>
#include <stdint.h>
#include <iostream>
#include <random>
#include <math.h>

/* Boost Includes */
#include <boost/chrono.hpp>
#include <boost/thread.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/container/vector.hpp>

/* Local Includes */

#include "ga_helper.h"


/* ValkyrieAI Includes */
#include "config.h"
#include "types.h"
#include "model.h"
#include "model_simulation.h"
#include "signal_analysis.h"

struct StateSpaceModelInput
{
	double dt;								/* Simulation dt step for each iteration of the solver (s) */
	double startTime;						/* Simulation start time (s) */
	double endTime;							/* Simulation end time (s) */
	
	GA_PIDChromosome<double> pid;			/* Specific PID values to use in the simulation */

	SS_ModelBase_sPtr model;				/* Generic State Space Model */
	ModelSimulationType simulationType;	/* Instructs the simulator what kind of simulation to execute */
};

struct StateSpaceModelOutput
{
	int errorCode;									/* Any possible error codes from the simulation */
	boost::chrono::microseconds executionTime;		/* Physical time spent solving for results */
	ModelSimulationType simulationType;			/* Kind of simulation that was run */
	StepPerformance_sPtr stepPerformance;			/* Calculated performance metrics given a step input */
};

struct NeuralNetworkModelInput
{
	double dt;								/* Simulation dt step for each iteration of the solver (s) */
	double startTime;						/* Simulation start time (s) */
	double endTime;							/* Simulation end time (s) */

	GA_PIDChromosome<double> pid;			/* Specific PID values to use in the simulation */

	NN_ModelBase_sPtr model;
	ModelSimulationType simulationType;

	float step_magnitude = 0.0;
};

struct NeuralNetworkModelOutput
{
	int errorCode;									/* Any possible error codes from the simulation */
	boost::chrono::microseconds executionTime;		/* Physical time spent solving for results */
	ModelSimulationType simulationType;				/* Kind of simulation that was run */
	StepPerformance_sPtr stepPerformance;			/* Calculated performance metrics given a step input */
};


class GA_EvaluateModelBase
{
public:
	virtual void evaluate(const StateSpaceModelInput, StateSpaceModelOutput&) {};
	virtual void evaluate(const NeuralNetworkModelInput, NeuralNetworkModelOutput&) {};

	virtual ~GA_EvaluateModelBase() = default;
private:

};
typedef boost::shared_ptr<GA_EvaluateModelBase> GA_EvaluateModelBase_sPtr;


///////////////////////////////////////////////////
/* CLASS:  StateSpaceModel */
///////////////////////////////////////////////////
class StateSpaceEvaluator : public GA_EvaluateModelBase
{
public:

	/**
	* \brief Evaluates model according to the given simulation type
	* The function is intended to be used for only a single set of PID values and thus 
	* only returns a single set of simulation data. This should allow for easy multi-threading
	* should the user want to do this. 
	*/
	void evaluate(const StateSpaceModelInput input, StateSpaceModelOutput& output) override;


	StateSpaceEvaluator() = default;
	~StateSpaceEvaluator() = default;

private:
	StateSpaceSimulator simulator;
	StepResponseAnalyzer stepAnalyzer;

};
typedef boost::shared_ptr<StateSpaceEvaluator> SSModel_sPtr;


///////////////////////////////////////////////////
/* CLASS:  NeuralNetworkModel */
///////////////////////////////////////////////////
class NeuralNetworkEvaluator : public GA_EvaluateModelBase
{
public:
	void evaluate(const NeuralNetworkModelInput input, NeuralNetworkModelOutput& output) override;

	NeuralNetworkEvaluator() = default;
	~NeuralNetworkEvaluator() = default;

private:

};
typedef boost::shared_ptr<NeuralNetworkEvaluator> NNModel_sPtr;

#endif 