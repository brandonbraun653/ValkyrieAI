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
		GA_BREED_SIMPLE_CROSSOVER |
		GA_BREED_SIMULATED_BINARY_CROSSOVER |
		GA_BREED_DYNAMIC_CROSSOVER |
		GA_BREED_FIXED_POINT_CROSSOVER);

	initStruct.solverParamOptions.allowedFitnessTypes = (
		GA_FITNESS_WEIGHTED_SUM |
		GA_FITNESS_MEAN_SQUARE_ERROR);

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
		GA_BREED_SIMPLE_CROSSOVER |
		GA_BREED_SIMULATED_BINARY_CROSSOVER |
		GA_BREED_DYNAMIC_CROSSOVER |
		GA_BREED_FIXED_POINT_CROSSOVER);

	initStruct.solverParamOptions.allowedFitnessTypes = (
		GA_FITNESS_WEIGHTED_SUM |
		GA_FITNESS_MEAN_SQUARE_ERROR);

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
		GA_BREED_SIMPLE_CROSSOVER |
		GA_BREED_SIMULATED_BINARY_CROSSOVER |
		GA_BREED_DYNAMIC_CROSSOVER |
		GA_BREED_FIXED_POINT_CROSSOVER);

	initStruct.solverParamOptions.allowedFitnessTypes = (
		GA_FITNESS_WEIGHTED_SUM |
		GA_FITNESS_MEAN_SQUARE_ERROR);

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


FCSOptimizer_Init_t rollTunerVer2_Init()
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
		GA_BREED_SIMPLE_CROSSOVER |
		GA_BREED_SIMULATED_BINARY_CROSSOVER |
		GA_BREED_DYNAMIC_CROSSOVER |
		GA_BREED_FIXED_POINT_CROSSOVER);

	initStruct.solverParamOptions.allowedFitnessTypes = (
		GA_FITNESS_WEIGHTED_SUM |
		GA_FITNESS_MEAN_SQUARE_ERROR);

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
	NN_JSONModel_sPtr model = boost::make_shared<NN_JSONModel>("C:\\Users\\Valkyrie\\Desktop\\TEMP\\");

	initStruct.neuralNetModel = boost::dynamic_pointer_cast<NN_ModelBase, NN_JSONModel>(model);
	initStruct.solverParam.modelType = GA_MODEL_NEURAL_NETWORK;


	return initStruct;
}


/* Benchmarking functions! */
//All options
FCSOptimizer_Init_t dynamicInit()
{
	/*-----------------------------
	* Setup all parameters for a single tuner
	*----------------------------*/
	FCSOptimizer_Init_t initStruct;
	initStruct.optimizerName = "DynamicOptimizer";
	initStruct.messageQueueName = "dynamic";
	initStruct.logPath = "dynamicOptimizer.txt";

	/*-----------------------------
	* STEP CHOOSER
	*----------------------------*/
	initStruct.solverParam.rngEngine = GA_MERSENNE_TWISTER;
	initStruct.solverParam.rngDistribution = GA_DISTRIBUTION_UNIFORM_REAL;
	initStruct.solverParam.breedType = GA_BREED_SIMPLE_CROSSOVER;
	initStruct.solverParam.fitnessType = GA_FITNESS_WEIGHTED_SUM;
	initStruct.solverParam.selectType = GA_SELECT_RANDOM;
	initStruct.solverParam.filterType = GA_POPULATION_STATIC_FILTER;
	initStruct.solverParam.resolutionType = GA_RESOLUTION_2DP;
	initStruct.solverParam.mutateType = GA_MUTATE_BIT_FLIP;
	initStruct.solverParam.mutateProbabilityType = GA_MUTATE_PROBABILITY_EXPONENTIAL;

	/*-----------------------------
	* DYNAMIC RECONFIGURATION OPTIONS
	*----------------------------*/
	initStruct.solverParamOptions.allowedParentSelectTypes = (
		GA_SELECT_RANDOM | 
		GA_SELECT_TOURNAMENT);

	initStruct.solverParamOptions.allowedBreedTypes = (
		GA_BREED_SIMPLE_CROSSOVER |
		GA_BREED_FIXED_POINT_CROSSOVER);

	initStruct.solverParamOptions.allowedFitnessTypes = (
		GA_FITNESS_WEIGHTED_SUM );

	initStruct.solverParamOptions.allowedMutateTypes = (
		GA_MUTATE_BIT_FLIP );

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
	initStruct.advConvergenceParam.iterations_before_refresh = 2;
	initStruct.advConvergenceParam.populationSize = 10;
	initStruct.advConvergenceParam.generationLimit = 100;
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


	return initStruct;
}

//Random, Simple Crossover
FCSOptimizer_Init_t staticInit1()
{
	/*-----------------------------
	* Setup all parameters for a single tuner
	*----------------------------*/
	FCSOptimizer_Init_t initStruct;
	initStruct.optimizerName = "static1";
	initStruct.messageQueueName = "static1";
	initStruct.logPath = "static1.txt";

	/*-----------------------------
	* STEP CHOOSER
	*----------------------------*/
	initStruct.solverParam.rngEngine = GA_MERSENNE_TWISTER;
	initStruct.solverParam.rngDistribution = GA_DISTRIBUTION_UNIFORM_REAL;
	initStruct.solverParam.breedType = GA_BREED_SIMPLE_CROSSOVER;
	initStruct.solverParam.fitnessType = GA_FITNESS_WEIGHTED_SUM;
	initStruct.solverParam.selectType = GA_SELECT_RANDOM;
	initStruct.solverParam.filterType = GA_POPULATION_STATIC_FILTER;
	initStruct.solverParam.resolutionType = GA_RESOLUTION_2DP;
	initStruct.solverParam.mutateType = GA_MUTATE_BIT_FLIP;
	initStruct.solverParam.mutateProbabilityType = GA_MUTATE_PROBABILITY_EXPONENTIAL;

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
	initStruct.advConvergenceParam.iterations_before_refresh = 2;
	initStruct.advConvergenceParam.populationSize = 10;
	initStruct.advConvergenceParam.generationLimit = 100;
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


	return initStruct;
}

//Tournament, Simple Crossover
FCSOptimizer_Init_t staticInit2()
{
	/*-----------------------------
	* Setup all parameters for a single tuner
	*----------------------------*/
	FCSOptimizer_Init_t initStruct;
	initStruct.optimizerName = "static2";
	initStruct.messageQueueName = "static2";
	initStruct.logPath = "static2.txt";

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
	initStruct.advConvergenceParam.iterations_before_refresh = 2;
	initStruct.advConvergenceParam.populationSize = 10;
	initStruct.advConvergenceParam.generationLimit = 100;
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


	return initStruct;
}

//Tournament, Fixed Point Crossover
FCSOptimizer_Init_t staticInit3()
{
	/*-----------------------------
	* Setup all parameters for a single tuner
	*----------------------------*/
	FCSOptimizer_Init_t initStruct;
	initStruct.optimizerName = "static3";
	initStruct.messageQueueName = "static3";
	initStruct.logPath = "static3.txt";

	/*-----------------------------
	* STEP CHOOSER
	*----------------------------*/
	initStruct.solverParam.rngEngine = GA_MERSENNE_TWISTER;
	initStruct.solverParam.rngDistribution = GA_DISTRIBUTION_UNIFORM_REAL;
	initStruct.solverParam.breedType = GA_BREED_FIXED_POINT_CROSSOVER;
	initStruct.solverParam.fitnessType = GA_FITNESS_WEIGHTED_SUM;
	initStruct.solverParam.selectType = GA_SELECT_TOURNAMENT;
	initStruct.solverParam.filterType = GA_POPULATION_STATIC_FILTER;
	initStruct.solverParam.resolutionType = GA_RESOLUTION_2DP;
	initStruct.solverParam.mutateType = GA_MUTATE_BIT_FLIP;
	initStruct.solverParam.mutateProbabilityType = GA_MUTATE_PROBABILITY_EXPONENTIAL;

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
	initStruct.advConvergenceParam.iterations_before_refresh = 2;
	initStruct.advConvergenceParam.populationSize = 10;
	initStruct.advConvergenceParam.generationLimit = 100;
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


	return initStruct;
}