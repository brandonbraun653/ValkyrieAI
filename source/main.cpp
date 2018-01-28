#pragma comment(linker, "/STACK:5000000")
#pragma comment(linker, "/HEAP:5000000")

#include <stdlib.h>
#include <iostream>
#include <thread>
#include <chrono>
#include <boost/chrono.hpp>
#include "valkyrie_engine.h"



int main()
{
	ValkyrieEngine DroneTuner;

	/*-----------------------------
	* Setup all parameters for a single tuner
	*----------------------------*/
	FCSOptimizer_Init_t initStruct;
	initStruct.optimizerName = "ROLL";
	initStruct.messageQueueName = "rollCMD";
	
	/*-----------------------------
	* STEP CHOOSER
	*----------------------------*/
	initStruct.solverParam.rngEngine				= GA_MERSENNE_TWISTER;
	initStruct.solverParam.rngDistribution			= GA_DISTRIBUTION_UNIFORM_REAL;
	initStruct.solverParam.breedType				= GA_BREED_SIMPLE_CROSSOVER;
	initStruct.solverParam.fitnessType				= GA_FITNESS_WEIGHTED_SUM;
	initStruct.solverParam.selectType				= GA_SELECT_RANDOM;
	initStruct.solverParam.filterType				= GA_POPULATION_STATIC_FILTER;
	initStruct.solverParam.resolutionType			= GA_RESOLUTION_2DP;
	initStruct.solverParam.mutateType				= GA_MUTATE_BIT_FLIP;
	initStruct.solverParam.mutateProbabilityType	= GA_MUTATE_PROBABILITY_EXPONENTIAL;

	/*-----------------------------
	* DYNAMIC RECONFIGURATION OPTIONS
	*----------------------------*/
	initStruct.solverParamOptions.allowedBreedTypes = (
		GA_BREED_SIMPLE_CROSSOVER_MSK |
		GA_BREED_SIMULATED_BINARY_CROSSOVER_MSK |
		GA_BREED_DYNAMIC_CROSSOVER_MSK |
		GA_BREED_FIXED_RATIO_CROSSOVER_MSK);

	initStruct.solverParamOptions.allowedFitnessTypes = (
		GA_FITNESS_WEIGHTED_SUM_MSK |
		GA_FITNESS_MEAN_SQUARE_ERROR_MSK);

	/*-----------------------------
	* PID SETTINGS
	*----------------------------*/
	initStruct.pidControlSettings.tuningLimits.Kp = { 0.0, 100.0 };
	initStruct.pidControlSettings.tuningLimits.Ki = { 0.0, 100.0 };
	initStruct.pidControlSettings.tuningLimits.Kd = { 0.0, 1.0 };

	initStruct.pidControlSettings.performanceGoals.percentOvershoot_goal = 0.10;
	initStruct.pidControlSettings.performanceGoals.riseTime_goal = 0.50;
	initStruct.pidControlSettings.performanceGoals.settlingTime_goal = 1.0;
	initStruct.pidControlSettings.performanceGoals.steadyStateError_goal = 0.05;

	initStruct.pidControlSettings.performanceTolerance.percentOvershoot_pcntTol = 0.1;
	initStruct.pidControlSettings.performanceTolerance.riseTime_pcntTol = 0.1;
	initStruct.pidControlSettings.performanceTolerance.settlingTime_pcntTol = 0.1;
	initStruct.pidControlSettings.performanceTolerance.steadyStateError_pcntTol = 0.1;

	/*-----------------------------
	* ADVANCED PARAMETER SETTINGS
	*----------------------------*/
	initStruct.advConvergenceParam.populationSize = 10;
	initStruct.advConvergenceParam.generationLimit = 250;
	initStruct.advConvergenceParam.limitingBehavior = FCS_LIMITER_REGENERATE_CHROMOSOME;

	/*-----------------------------
	* SIMULATION MODEL
	*----------------------------*/
	SS_NLTIVModel_sPtr model = boost::make_shared<SS_NLTIVModel>(1, 1, 2);
	model->A << -8.202, -2.029, -0.149, -3.25;
	model->B << 1.14, -1.23;
	model->C << 1.0, 0.0;
	model->D << 0.0;
	model->X0 << 0.0, 0.0;

	//Testing: Purposefully cast back to the base class for abstraction into the model simulation method
	initStruct.stateSpaceModel = boost::dynamic_pointer_cast<SS_ModelBase, SS_NLTIVModel>(model);
	initStruct.solverParam.modelType = GA_MODEL_STATE_SPACE;

	/*-----------------------------
	* RUN
	*----------------------------*/
	FCSOptimizer_Handle hRollTuner = DroneTuner.newOptimizer(initStruct);

	DroneTuner.initialize(hRollTuner);
	DroneTuner.start(hRollTuner);

	DroneTuner.waitForCompletion();
	return 0;
}