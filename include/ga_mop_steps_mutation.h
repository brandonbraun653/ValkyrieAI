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

/* Local Includes */
#include "ga_config.h"
#include "ga_helper.h"
#include "data.h"

///////////////////////////////////////////////////
/* CLASS:  MutateProbGenerator */
///////////////////////////////////////////////////
class MutateProbGenerator
{
public:
	double get();

	MutateProbGenerator(GA_MutateProbabilityMethod prob_type);
	~MutateProbGenerator();

private:
	GA_MutateProbabilityMethod probType;
};

///////////////////////////////////////////////////
/* CLASS:  BitFlipMutator */
///////////////////////////////////////////////////
class BitFlipMutator
{
public:
	void mutate(hPID_Chromosomes* in_bredChrom, GA_ConverganceCriteria_sPtr in_convgCriteria, PID_ControlGoals_sPtr in_config,
		hPID_DataVector* out_pidVals, mapCoeff_t* in_kp, mapCoeff_t* in_ki, mapCoeff_t* in_kd);

	BitFlipMutator(GA_RunMode execution_type, GA_MutateProbabilityMethod mutationProb_type);
	~BitFlipMutator();
private:
	GA_RunMode executionType;
	GA_MutateProbabilityMethod mutationProbType;

	hPID_DataVector* cm_data;
	hPID_Chromosomes* cm_chrom;
	GA_ConverganceCriteria_sPtr cm_convergence;
	PID_ControlGoals_sPtr cm_config;

	mapCoeff_t* mapCM_Kp;
	mapCoeff_t* mapCM_Ki;
	mapCoeff_t* mapCM_Kd;

	void calculate_cpu_single_threaded();
	void calculate_cpu_multi_threaded();
	void calculate_gpu_single_threaded();
	void calculate_gpu_multi_threaded();

	void mutateKp();
	void mutateKi();
	void mutateKd();
};

///////////////////////////////////////////////////
/* CLASS:  AddSubMutator */
///////////////////////////////////////////////////
class AddSubMutator
{
public:
	void mutate(hPID_Chromosomes* in_bredChrom, GA_ConverganceCriteria_sPtr in_convgCriteria, PID_ControlGoals_sPtr in_config,
		hPID_DataVector* out_pidVals, mapCoeff_t* in_kp, mapCoeff_t* in_ki, mapCoeff_t* in_kd);

	AddSubMutator(GA_RunMode execution_type, GA_MutateProbabilityMethod mutationProb_type, GA_Resolution res_type);
	~AddSubMutator();
private:
	GA_RunMode executionType;
	GA_MutateProbabilityMethod mutationProbType;
	GA_Resolution resolutionType;

	hPID_DataVector* as_data;
	hPID_Chromosomes* as_chrom;
	GA_ConverganceCriteria_sPtr as_convergence;
	PID_ControlGoals_sPtr as_config;

	mapCoeff_t* mapAS_Kp;
	mapCoeff_t* mapAS_Ki;
	mapCoeff_t* mapAS_Kd;

	void calculate_cpu_single_threaded();
	void calculate_cpu_multi_threaded();
	void calculate_gpu_single_threaded();
	void calculate_gpu_multi_threaded();

	void mutateKp();
	void mutateKi();
	void mutateKd();
};
#endif