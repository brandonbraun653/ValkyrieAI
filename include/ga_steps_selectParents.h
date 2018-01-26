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
#include <boost/shared_ptr.hpp>

/* Local Includes */
#include "config.h"
#include "ga_helper.h"
#include "types.h"
#include "rng.hpp"

struct GA_SelectParentDataInput
{
	size_t populationSize = 0;
	boost::container::vector<double> popGlobalFitScores;	/* Latest fitness scores for all population members */
	//Do I need a mutex?
	//Add more as needed 
};

struct GA_SelectParentDataOutput
{
	/* Ordered list of which parents mate together. Should be 2x population size
	so that parent 1 mates with parent 2, 3 with 4, 5 with 6, etc */
	boost::container::vector<int> parentSelections;	
};

class GA_SelectParentBase
{
public:
	virtual void selectParent(const GA_SelectParentDataInput, GA_SelectParentDataOutput&) = 0;

	virtual ~GA_SelectParentBase() = default;
private:
	virtual void selectParentKp() = 0;
	virtual void selectParentKi() = 0;
	virtual void selectParentKd() = 0;
};
typedef boost::shared_ptr<GA_SelectParentBase> GA_SelectParentBase_sPtr;


///////////////////////////////////////////////////
/* CLASS:  RankedSelection */
///////////////////////////////////////////////////
class RankedSelection : public GA_SelectParentBase
{
public:
	void selectParent(const GA_SelectParentDataInput input, GA_SelectParentDataOutput& output) override;
	
	RankedSelection();
	~RankedSelection();
private:
	void selectParentKp() override;
	void selectParentKi() override;
	void selectParentKd() override;
};

///////////////////////////////////////////////////
/* CLASS:  RandomSelection */
///////////////////////////////////////////////////
class RandomSelection : public GA_SelectParentBase
{
public:
	void selectParent(const GA_SelectParentDataInput input, GA_SelectParentDataOutput& output) override;

	RandomSelection(const int populationSize);
	~RandomSelection();
private:
	void selectParentKp() override; //Not actually used in this class 
	void selectParentKi() override; //Not actually used in this class
	void selectParentKd() override; //Not actually used in this class

	RNGManager_sPtr rng_engine;
};

///////////////////////////////////////////////////
/* CLASS:  RouletteSelection */
///////////////////////////////////////////////////
class RouletteSelection : public GA_SelectParentBase
{
public:
	void selectParent(const GA_SelectParentDataInput input, GA_SelectParentDataOutput& output) override;

	RouletteSelection();
	~RouletteSelection();
private:
	void selectParentKp() override;
	void selectParentKi() override;
	void selectParentKd() override;
};

///////////////////////////////////////////////////
/* CLASS:  StochasticSelection */
///////////////////////////////////////////////////
class StochasticSelection : public GA_SelectParentBase
{
public:
	void selectParent(const GA_SelectParentDataInput input, GA_SelectParentDataOutput& output) override;

	StochasticSelection();
	~StochasticSelection();
private:
	void selectParentKp() override;
	void selectParentKi() override;
	void selectParentKd() override;
};

///////////////////////////////////////////////////
/* CLASS:  TournamentSelection */
///////////////////////////////////////////////////
class TournamentSelection : public GA_SelectParentBase
{
public:
	void selectParent(const GA_SelectParentDataInput input, GA_SelectParentDataOutput& output) override;

	TournamentSelection(const int populationSize);
	~TournamentSelection();
private:
	void selectParentKp() override;
	void selectParentKi() override;
	void selectParentKd() override;

	RNGManager_sPtr tourneySizeSelectorRNG;
	RNGManager_sPtr tourneyCompetitorSelectorRNG;
};

///////////////////////////////////////////////////
/* CLASS:  ElitistSelection */
///////////////////////////////////////////////////
class ElitistSelection : public GA_SelectParentBase
{
public:
	void selectParent(const GA_SelectParentDataInput input, GA_SelectParentDataOutput& output) override;

	ElitistSelection();
	~ElitistSelection();
private:
	void selectParentKp() override;
	void selectParentKi() override;
	void selectParentKd() override;
};

#endif