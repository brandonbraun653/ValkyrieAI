#include "../include/user_optimizer_config.h"


FCSOptimizer_Init_t rollTunerVer1_Init()
{
	/*-----------------------------
	* Setup all parameters for a single tuner
	*----------------------------*/
	FCSOptimizer_Init_t initStruct;
	initStruct.optimizerName = "ROLL";
	initStruct.messageQueueName = "rollCMD";
	initStruct.logPath = "rollVer1Log.txt";

	/*-----------------------------
	* STEP CHOOSER
	*----------------------------*/
	initStruct.solverParam.rngEngine = GA_MERSENNE_TWISTER;
	initStruct.solverParam.rngDistribution = GA_DISTRIBUTION_UNIFORM_REAL;
	initStruct.solverParam.breedType = GA_BREED_SIMPLE_CROSSOVER;
	initStruct.solverParam.fitnessType = GA_FITNESS_WEIGHTED_SUM;
	initStruct.solverParam.selectType = GA_SELECT_TOURNAMENT;
	initStruct.solverParam.filterType = GA_POPULATION_STATIC_FILTER;
	initStruct.solverParam.resolutionType = GA_RESOLUTION_2DP;
	initStruct.solverParam.mutateType = GA_MUTATE_BIT_FLIP;
	initStruct.solverParam.mutateProbabilityType = GA_MUTATE_PROBABILITY_EXPONENTIAL;

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
	initStruct.pidControlSettings.tuningLimits.Kd = { 0.0, 100.0 };

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
	initStruct.advConvergenceParam.iterations_before_refresh = 5;
	initStruct.advConvergenceParam.populationSize = 10;
	initStruct.advConvergenceParam.generationLimit = 15;
	initStruct.advConvergenceParam.limitingBehavior = FCS_LIMITER_REGENERATE_CHROMOSOME;

	/*-----------------------------
	* SIMULATION MODEL
	*----------------------------*/
	MatlabModel_sPtr model = boost::make_shared<MatlabModel>();
	model->setModelRoot("C:\\git\\GitHub\\ValkyrieAI\\Matlab");
	model->setInitFunction("quad_initialization");
	model->setModelFunction("quad_simulation");

	model->simAxis = "roll";
	model->endTime = 10.0;
	model->stepMagnitude = 10.0;
	model->stepTimeEnable = 0.5;


	initStruct.matlabModel = boost::dynamic_pointer_cast<ML_ModelBase, MatlabModel>(model);
	initStruct.solverParam.modelType = GA_MODEL_MATLAB;


	return initStruct;
}

FCSOptimizer_Init_t pitchTunerVer1_Init()
{
	/*-----------------------------
	* Setup all parameters for a single tuner
	*----------------------------*/
	FCSOptimizer_Init_t initStruct;
	initStruct.optimizerName = "PITCH";
	initStruct.messageQueueName = "pitchCMD";
	initStruct.logPath = "pitchVer1Log.txt";

	/*-----------------------------
	* STEP CHOOSER
	*----------------------------*/
	initStruct.solverParam.rngEngine = GA_MERSENNE_TWISTER;
	initStruct.solverParam.rngDistribution = GA_DISTRIBUTION_UNIFORM_REAL;
	initStruct.solverParam.breedType = GA_BREED_SIMPLE_CROSSOVER;
	initStruct.solverParam.fitnessType = GA_FITNESS_WEIGHTED_SUM;
	initStruct.solverParam.selectType = GA_SELECT_TOURNAMENT;
	initStruct.solverParam.filterType = GA_POPULATION_STATIC_FILTER;
	initStruct.solverParam.resolutionType = GA_RESOLUTION_2DP;
	initStruct.solverParam.mutateType = GA_MUTATE_BIT_FLIP;
	initStruct.solverParam.mutateProbabilityType = GA_MUTATE_PROBABILITY_EXPONENTIAL;

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
	initStruct.pidControlSettings.tuningLimits.Kd = { 0.0, 100.0 };

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
	initStruct.advConvergenceParam.iterations_before_refresh = 5;
	initStruct.advConvergenceParam.populationSize = 10;
	initStruct.advConvergenceParam.generationLimit = 15;
	initStruct.advConvergenceParam.limitingBehavior = FCS_LIMITER_REGENERATE_CHROMOSOME;

	/*-----------------------------
	* SIMULATION MODEL
	*----------------------------*/
	MatlabModel_sPtr model = boost::make_shared<MatlabModel>();
	model->setModelRoot("C:\\git\\GitHub\\ValkyrieAI\\Matlab");
	model->setInitFunction("quad_initialization");
	model->setModelFunction("quad_simulation");

	model->simAxis = "pitch";
	model->endTime = 10.0;
	model->stepMagnitude = 10.0;
	model->stepTimeEnable = 0.5;


	initStruct.matlabModel = boost::dynamic_pointer_cast<ML_ModelBase, MatlabModel>(model);
	initStruct.solverParam.modelType = GA_MODEL_MATLAB;


	return initStruct;
}

FCSOptimizer_Init_t yawTunerVer1_Init()
{
	/*-----------------------------
	* Setup all parameters for a single tuner
	*----------------------------*/
	FCSOptimizer_Init_t initStruct;
	initStruct.optimizerName = "YAW";
	initStruct.messageQueueName = "yawCMD";
	initStruct.logPath = "yawVer1Log.txt";

	/*-----------------------------
	* STEP CHOOSER
	*----------------------------*/
	initStruct.solverParam.rngEngine = GA_MERSENNE_TWISTER;
	initStruct.solverParam.rngDistribution = GA_DISTRIBUTION_UNIFORM_REAL;
	initStruct.solverParam.breedType = GA_BREED_SIMPLE_CROSSOVER;
	initStruct.solverParam.fitnessType = GA_FITNESS_WEIGHTED_SUM;
	initStruct.solverParam.selectType = GA_SELECT_TOURNAMENT;
	initStruct.solverParam.filterType = GA_POPULATION_STATIC_FILTER;
	initStruct.solverParam.resolutionType = GA_RESOLUTION_2DP;
	initStruct.solverParam.mutateType = GA_MUTATE_BIT_FLIP;
	initStruct.solverParam.mutateProbabilityType = GA_MUTATE_PROBABILITY_EXPONENTIAL;

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
	initStruct.pidControlSettings.tuningLimits.Kd = { 0.0, 100.0 };

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
	initStruct.advConvergenceParam.iterations_before_refresh = 5;
	initStruct.advConvergenceParam.populationSize = 10;
	initStruct.advConvergenceParam.generationLimit = 15;
	initStruct.advConvergenceParam.limitingBehavior = FCS_LIMITER_REGENERATE_CHROMOSOME;

	/*-----------------------------
	* SIMULATION MODEL
	*----------------------------*/
	MatlabModel_sPtr model = boost::make_shared<MatlabModel>();
	model->setModelRoot("C:\\git\\GitHub\\ValkyrieAI\\Matlab");
	model->setInitFunction("quad_initialization");
	model->setModelFunction("quad_simulation");

	model->simAxis = "yaw";
	model->endTime = 10.0;
	model->stepMagnitude = 10.0;
	model->stepTimeEnable = 0.5;


	initStruct.matlabModel = boost::dynamic_pointer_cast<ML_ModelBase, MatlabModel>(model);
	initStruct.solverParam.modelType = GA_MODEL_MATLAB;


	return initStruct;
}
