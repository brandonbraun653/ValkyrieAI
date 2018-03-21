#include <stdlib.h>
#include <iostream>
#include <thread>
#include <chrono>
#include <boost/chrono.hpp>

#include <nlohmann/json.hpp>

#include "valkyrie_engine.h"
#include "user_optimizer_config.h"

#include "logger.h"


int main()
{
	/*-----------------------------
	* CREATE TUNING MANAGEMENT ENGINE
	*----------------------------*/
	ValkyrieEngine DroneTuner;

	/*-----------------------------
	* CREATE OPTIMIZERS
	*----------------------------*/
	DroneTuner.newOptimizer(dynamicInit());
	//DroneTuner.newOptimizer(staticInit1());
	//DroneTuner.newOptimizer(staticInit2());
	//DroneTuner.newOptimizer(staticInit3());


	/*-----------------------------
	* RUN
	*----------------------------*/
	DroneTuner.initAll();
	DroneTuner.startAll();

	DroneTuner.waitForCompletion();
	return 0;
}