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
#include "ga_config.h"
#include "ga_helper.h"
#include "types.h"

struct GA_SelectParentDataInput
{
	//Add whatever the heck is needed here
};

struct GA_SelectParentDataOutput
{
	//Add whatever the heck is needed here
};

class GA_SelectParentBase
{
public:
	virtual void selectParent(const GA_SelectParentDataInput, GA_SelectParentDataOutput&) = 0;

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

	RandomSelection();
	~RandomSelection();
private:
	void selectParentKp() override;
	void selectParentKi() override;
	void selectParentKd() override;
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

	TournamentSelection();
	~TournamentSelection();
private:
	void selectParentKp() override;
	void selectParentKi() override;
	void selectParentKd() override;
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