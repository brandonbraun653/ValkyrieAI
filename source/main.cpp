#include <stdlib.h>
#include <iostream>
#include <thread>
#include <chrono>
#include "valkyrie_engine.h"



int main()
{
	ValkyrieEngine DroneTuner;

	FCSOptimizer_Init_t initStruct;

	initStruct.optimizerName = "ROLL";
	initStruct.messageQueueName = "rollCMD";
	
	initStruct.solverParam.rngEngine = GA_MERSENNE_TWISTER;
	initStruct.solverParam.rngDistribution = GA_DISTRIBUTION_UNIFORM_REAL;

	FCSOptimizer_Handle hRollTuner = DroneTuner.newOptimizer(initStruct);

	DroneTuner.initialize(hRollTuner);
	DroneTuner.start(hRollTuner);

	std::this_thread::sleep_for(std::chrono::seconds(2));

	DroneTuner.stop(hRollTuner);

	DroneTuner.waitForCompletion();
	return 0;
}