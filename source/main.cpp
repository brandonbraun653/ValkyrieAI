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
	FCSOptimizer_Handle hRollTuner = DroneTuner.newOptimizer(rollTunerVer1_Init());
	FCSOptimizer_Handle hPitchTuner = DroneTuner.newOptimizer(pitchTunerVer1_Init());
	FCSOptimizer_Handle hYawTuner = DroneTuner.newOptimizer(yawTunerVer1_Init());


	/*-----------------------------
	* RUN
	*----------------------------*/
	DroneTuner.initAll();
	DroneTuner.startAll();

	DroneTuner.waitForCompletion();
	return 0;
}