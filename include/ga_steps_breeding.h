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
#include "ga_config.h"
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
	virtual void breed(GA_BreedingDataInput, GA_BreedingDataOutput&) = 0;

private:
	virtual void breedKp() = 0;
	virtual void breedKi() = 0;
	virtual void breedKd() = 0;
};
typedef boost::shared_ptr<GA_BreedBase> GA_BreedBase_sPtr;


#endif /* BREEDING_H_ */
// ///////////////////////////////////////////////////
// /* CLASS:  SimpleCrossover */
// ///////////////////////////////////////////////////
// class SimpleCrossover
// {
// public:
// 	void breed(iVec in_parents, hPID_DataVector* in_pidVals, hPID_Chromosomes* out_breed, FCSOptimizer_MappingCoeff* in_kp, FCSOptimizer_MappingCoeff* in_ki, FCSOptimizer_MappingCoeff* in_kd);
// 
// 	SimpleCrossover();
// 	~SimpleCrossover();
// private:
// 
// 	iVec sc_parents;
// 	hPID_DataVector* sc_data;
// 	hPID_Chromosomes* sc_chrom;
// 
// 	FCSOptimizer_MappingCoeff* mapCF_Kp;
// 	FCSOptimizer_MappingCoeff* mapCF_Ki;
// 	FCSOptimizer_MappingCoeff* mapCF_Kd;
// 
// 	void calculate_cpu_single_threaded();
// 	void calculate_cpu_multi_threaded();
// 	void calculate_gpu_single_threaded();
// 	void calculate_gpu_multi_threaded();
// 
// 	void breedKp();
// 	void breedKi();
// 	void breedKd();
// };
// 
// 
// ///////////////////////////////////////////////////
// /* CLASS:  DynamicCrossover */
// ///////////////////////////////////////////////////
// class DynamicCrossover
// {
// public:
// 	void breed(iVec in_parents, hPID_DataVector* in_pidVals, hPID_Chromosomes* out_breed, FCSOptimizer_MappingCoeff* in_kp, FCSOptimizer_MappingCoeff* in_ki, FCSOptimizer_MappingCoeff* in_kd);
// 
// 	DynamicCrossover();
// 	~DynamicCrossover();
// 
// private:
// 
// 	iVec dc_parents;
// 	hPID_DataVector* dc_data;
// 	hPID_Chromosomes* dc_chrom;
// 
// 	FCSOptimizer_MappingCoeff* mapCF_Kp;
// 	FCSOptimizer_MappingCoeff* mapCF_Ki;
// 	FCSOptimizer_MappingCoeff* mapCF_Kd;
// 
// 	void calculate_cpu_single_threaded();
// 	void calculate_cpu_multi_threaded();
// 	void calculate_gpu_single_threaded();
// 	void calculate_gpu_multi_threaded();
// 
// 	void breedKp();
// 	void breedKi();
// 	void breedKd();
// 
// };
// 
// ///////////////////////////////////////////////////
// /* CLASS:  FixedRatioCrossover */
// ///////////////////////////////////////////////////
// class FixedRatioCrossover
// {
// public:
// 	void breed(iVec in_parents, hPID_DataVector* in_pidVals, hPID_Chromosomes* out_breed, FCSOptimizer_MappingCoeff* in_kp, FCSOptimizer_MappingCoeff* in_ki, FCSOptimizer_MappingCoeff* in_kd);
// 
// 	FixedRatioCrossover(double ratio);
// 	~FixedRatioCrossover();
// 
// private:
// 	double fr_ratio;
// 	iVec fr_parents;
// 	hPID_DataVector* fr_data;
// 	hPID_Chromosomes* fr_chrom;
// 
// 	FCSOptimizer_MappingCoeff* mapCF_Kp;
// 	FCSOptimizer_MappingCoeff* mapCF_Ki;
// 	FCSOptimizer_MappingCoeff* mapCF_Kd;
// 
// 	void calculate_cpu_single_threaded();
// 	void calculate_cpu_multi_threaded();
// 	void calculate_gpu_single_threaded();
// 	void calculate_gpu_multi_threaded();
// 
// 	void breedKp();
// 	void breedKi();
// 	void breedKd();
// };
// 
// ///////////////////////////////////////////////////
// /* CLASS:  SimulatedBinaryCrossover */
// ///////////////////////////////////////////////////
// class SimulatedBinaryCrossover
// {
// public:
// 
// 	SimulatedBinaryCrossover();
// 	~SimulatedBinaryCrossover();
// private:
// 	void calculate_cpu_single_threaded();
// 	void calculate_cpu_multi_threaded();
// 	void calculate_gpu_single_threaded();
// 	void calculate_gpu_multi_threaded();
// };
// 
// #endif