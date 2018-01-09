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
#include <boost/thread.hpp>
#include <boost/shared_ptr.hpp>

/* Local Includes */
#include "ga_config.h"
#include "ga_helper.h"
#include "types.h"

struct GA_EvaluateModelDataInput
{
	//Add whatever the heck is needed here
};

struct GA_EvaluateModelDataOutput
{
	//Add whatever the heck is needed here
};

class GA_EvaluateModelBase
{
public:
	virtual void evaluate(const GA_EvaluateModelDataInput, GA_EvaluateModelDataOutput&) = 0;

private:

};
typedef boost::shared_ptr<GA_EvaluateModelBase> GA_EvaluateModelBase_sPtr;


///////////////////////////////////////////////////
/* CLASS:  StateSpaceModel */
///////////////////////////////////////////////////
class StateSpaceModel : public GA_EvaluateModelBase
{
public:
	virtual void evaluate(const GA_EvaluateModelDataInput input, GA_EvaluateModelDataOutput& output) override;

	StateSpaceModel();
	~StateSpaceModel();

private:

};
typedef boost::shared_ptr<StateSpaceModel> SSModel_sPtr;

///////////////////////////////////////////////////
/* CLASS:  NeuralNetworkModel */
///////////////////////////////////////////////////
class NeuralNetworkModel : public GA_EvaluateModelBase
{
public:
	virtual void evaluate(const GA_EvaluateModelDataInput input, GA_EvaluateModelDataOutput& output) override;

	NeuralNetworkModel();
	~NeuralNetworkModel();

private:

};
typedef boost::shared_ptr<NeuralNetworkModel> NNModel_sPtr;

#endif 