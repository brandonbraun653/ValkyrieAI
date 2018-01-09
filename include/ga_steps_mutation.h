#pragma once
#ifndef MUTATION_H_
#define MUTATION_H_

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


struct GA_MutateDataInput
{
	//Add whatever the heck is needed here
};

struct GA_MutateDataOutput
{
	//Add whatever the heck is needed here
};


class GA_MutateBase
{
public:
	virtual void mutate(const GA_MutateDataInput, GA_MutateDataOutput&) = 0;

private:
	virtual void mutateKp() = 0;
	virtual void mutateKi() = 0;
	virtual void mutateKd() = 0;

};
typedef boost::shared_ptr<GA_MutateBase> GA_MutateBase_sPtr;

///////////////////////////////////////////////////
/* CLASS:  MutateProbGenerator */
///////////////////////////////////////////////////
// class MutateProbGenerator
// {
// public:
// 	double get();
// 
// 	MutateProbGenerator(GA_METHOD_MutateProbability prob_type);
// 	~MutateProbGenerator();
// 
// private:
// 	GA_METHOD_MutateProbability probType;
// };

///////////////////////////////////////////////////
/* CLASS:  BitFlipMutator */
///////////////////////////////////////////////////
class BitFlipMutator : public GA_MutateBase
{
public:
	void mutate(const GA_MutateDataInput input, GA_MutateDataOutput& output) override;

	BitFlipMutator();
	~BitFlipMutator();
private:

	void mutateKp() override;
	void mutateKi() override;
	void mutateKd() override;
};

///////////////////////////////////////////////////
/* CLASS:  AddSubMutator */
///////////////////////////////////////////////////
class AddSubMutator : public GA_MutateBase
{
public:
	void mutate(const GA_MutateDataInput input, GA_MutateDataOutput& output) override;

	AddSubMutator();
	~AddSubMutator();
private:
	
	void mutateKp() override;
	void mutateKi() override;
	void mutateKd() override;
};
#endif