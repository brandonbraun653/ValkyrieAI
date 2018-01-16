#pragma once
#ifndef BREEDING_H_
#define BREEDING_H_

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
 
struct GA_BreedingDataInput
{
	//Add whatever the heck is needed here
};

struct GA_BreedingDataOutput
{
	//Add whatever the heck is needed here
};

class GA_BreedBase
{
public:
	virtual void breed(const GA_BreedingDataInput, GA_BreedingDataOutput&) = 0;

private:
	virtual void breedKp() = 0;
	virtual void breedKi() = 0;
	virtual void breedKd() = 0;
};
typedef boost::shared_ptr<GA_BreedBase> GA_BreedBase_sPtr;


///////////////////////////////////////////////////
/* CLASS:  SimpleCrossover */
///////////////////////////////////////////////////
class SimpleCrossover : public GA_BreedBase
{
public:

	void breed(const GA_BreedingDataInput input, GA_BreedingDataOutput& output) override;

	SimpleCrossover();
	~SimpleCrossover();

private:
	void breedKp() override;
	void breedKi() override;
	void breedKd() override;

};

///////////////////////////////////////////////////
/* CLASS:  DynamicCrossover */
///////////////////////////////////////////////////
class DynamicCrossover : public GA_BreedBase
{
public:

	void breed(const GA_BreedingDataInput input, GA_BreedingDataOutput& output) override;

	DynamicCrossover();
	~DynamicCrossover();

private:
	void breedKp() override;
	void breedKi() override;
	void breedKd() override;

};

///////////////////////////////////////////////////
/* CLASS:  FixedRatioCrossover */
///////////////////////////////////////////////////
class FixedRatioCrossover : public GA_BreedBase
{
public:

	void breed(const GA_BreedingDataInput input, GA_BreedingDataOutput& output) override;

	FixedRatioCrossover();
	~FixedRatioCrossover();

private:
	void breedKp() override;
	void breedKi() override;
	void breedKd() override;

};

///////////////////////////////////////////////////
/* CLASS:  SimulatedBinaryCrossover */
///////////////////////////////////////////////////
class SimulatedBinaryCrossover : public GA_BreedBase
{
public:

	void breed(const GA_BreedingDataInput input, GA_BreedingDataOutput& output) override;

	SimulatedBinaryCrossover();
	~SimulatedBinaryCrossover();

private:
	void breedKp() override;
	void breedKi() override;
	void breedKd() override;

};

#endif /* BREEDING_H_ */