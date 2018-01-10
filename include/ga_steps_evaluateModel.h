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

/* Local Includes */

#include "ga_helper.h"


/* ValkyrieAI Includes */
#include "ga_config.h"
#include "types.h"
#include "model_simulation.h"


struct StateSpaceModelInput
{
	double dt;				/* Simulation dt step for each iteration of the solver (s) */
	double startTime;		/* Simulation start time (s) */
	double endTime;			/* Simulation end time (s) */
	
	PID_Values pid;

	SS_NLTIVModel model;		
	StateSpaceSimulation simulationType;
};

struct StateSpaceModelOutput
{
	Eigen::MatrixXd data;	/* Raw simulation data */

	int errorCode;			/* Any possible error codes from the simulation */
};

struct NeuralNetworkModelInput
{

};

struct NeuralNetworkModelOutput
{

};


class GA_EvaluateModelBase
{
public:
	virtual void evaluate(const StateSpaceModelInput, StateSpaceModelOutput&) {};
	virtual void evaluate(const NeuralNetworkModelInput, NeuralNetworkModelOutput&) {};
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


	StateSpaceEvaluator();
	~StateSpaceEvaluator();

private:
	StateSpaceSimulator_sPtr simulator;

};
typedef boost::shared_ptr<StateSpaceEvaluator> SSModel_sPtr;


///////////////////////////////////////////////////
/* CLASS:  NeuralNetworkModel */
///////////////////////////////////////////////////
class NeuralNetworkEvaluator : public GA_EvaluateModelBase
{
public:
	void evaluate(const NeuralNetworkModelInput input, NeuralNetworkModelOutput& output) override;

	NeuralNetworkEvaluator();
	~NeuralNetworkEvaluator();

private:

};
typedef boost::shared_ptr<NeuralNetworkEvaluator> NNModel_sPtr;

#endif 