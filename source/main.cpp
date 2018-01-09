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
	
	initStruct.solverParam.rngEngine = GA_MERSENNE_TWISTER;
	initStruct.solverParam.rngDistribution = GA_DISTRIBUTION_UNIFORM_REAL;

	initStruct.pidControlSettings.tuningLimits.Kp = { 0.0, 100.0 };
	initStruct.pidControlSettings.tuningLimits.Ki = { 0.0, 100.0 };
	initStruct.pidControlSettings.tuningLimits.Kd = { 0.0, 100.0 };


	/* Setup the starting operation parameters */
	initStruct.solverParam.breedType = GA_BREED_SIMPLE_CROSSOVER;

	/* Choose which parameters the reconfiguration engine is allowed to use */
	initStruct.solverParamOptions.allowedBreedTypes = (
		GA_BREED_SIMPLE_CROSSOVER_MSK | 
		GA_BREED_SIMULATED_BINARY_CROSSOVER_MSK | 
		GA_BREED_DYNAMIC_CROSSOVER_MSK |
		GA_BREED_FIXED_RATIO_CROSSOVER_MSK );

	FCSOptimizer_Handle hRollTuner = DroneTuner.newOptimizer(initStruct);

	DroneTuner.initialize(hRollTuner);
	DroneTuner.start(hRollTuner);

	std::this_thread::sleep_for(std::chrono::seconds(2));

	DroneTuner.stop(hRollTuner);

	DroneTuner.waitForCompletion();
	return 0;
}