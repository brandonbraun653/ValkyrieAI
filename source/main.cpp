#include <stdlib.h>
#include <iostream>
#include "valkyrie_engine.h"



int main()
{
	ValkyrieEngine DroneTuner;

	FCSOptimizer_Init_t initStruct;

	FCSOptimizer_Handle hRollTuner = DroneTuner.newOptimizer(initStruct);

	DroneTuner.start(hRollTuner);


	DroneTuner.waitForCompletion();
	return 0;
}