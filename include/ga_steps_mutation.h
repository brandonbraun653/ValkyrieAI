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
#include "config.h"
#include "ga_helper.h"
#include "types.h"
#include "rng.hpp"


struct GA_MutateDataInput
{
	double mutationProbabilityThreshold = 0.5;							/* Range from 0.0 - 1.0 that determines how likely it is for a mutation to occur */
	GA_METHOD_MutateProbability mutateProbType;							/* Defines probability distribution curve to use when deciding to mutate */


	GA_ChromosomeMappingType optimizerChromType;						/* Mutation is last step, so it needs to known to convert back to optimizer data type */
	GA_ChromosomeMappingType chromType;									/* Informs the user what type of chromosome was inputted */

	boost::container::vector<GA_PIDChromosome<double>> d_chrom;			/* Real valued representation of chromosomes*/
	boost::container::vector<GA_PIDChromosome<uint16_t>> u16_chrom;		/* Bit field mapped chromosomes */

	FCSOptimizer_MappingCoeff* mapCoeff_Kp;								/* Precalculated mapping coefficients to convert between chrom types */
	FCSOptimizer_MappingCoeff* mapCoeff_Ki;
	FCSOptimizer_MappingCoeff* mapCoeff_Kd;
};

struct GA_MutateDataOutput
{
	GA_ChromosomeMappingType chromType;									/* Informs the user as to how the chromosomes are being represented */

	boost::container::vector<GA_PIDChromosome<double>> d_chrom;			/* Real valued representation of chromosomes */
	boost::container::vector<GA_PIDChromosome<uint16_t>> u16_chrom;		/* Bit field mapped chromosomes */
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
class MutateProbGenerator
{
public:
	double get(GA_METHOD_MutateProbability probType);

	MutateProbGenerator();
	~MutateProbGenerator();

private:
	RNGManager_sPtr rngUniform;
};

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
	const int maxSeverity = 15;	//Set to 15 due to using uint16_t as chrom type
	RNGManager_sPtr severityRNG;
	boost::shared_ptr<MutateProbGenerator> mutateRNG;

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