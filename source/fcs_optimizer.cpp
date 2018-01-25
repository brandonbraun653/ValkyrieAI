#include "fcs_optimizer.h"
#include <string>

using namespace boost::interprocess;

typedef boost::random::uniform_real_distribution<> urd_t;
typedef boost::random::uniform_int_distribution<> uid_t;

//////////////////////////////////////////////////////////////////
/* Helper Functions */
//////////////////////////////////////////////////////////////////
void calculateMappingCoefficients(FCSOptimizer_MappingCoeff *mapping, double lower, double upper)
{
	/* Currently hard coded precision, but I may want to template this and the output types later */
	mapping->bytes_precision = 2.0;

	mapping->x_lo = lower;
	mapping->x_hi = upper;

	mapping->x_dPow = 8 * mapping->bytes_precision - 1;
	mapping->x_bPow = ceil(log2(mapping->x_hi - mapping->x_lo));

	mapping->x_sF = pow(2, (mapping->x_dPow - mapping->x_bPow));
	mapping->x_sR = pow(2, (mapping->x_bPow - mapping->x_dPow));

	mapping->x_offset = 0.0;
	if (mapping->x_lo < 0.0 && mapping->x_hi > 0.0)
		mapping->x_offset = (mapping->x_sF*mapping->x_lo) - floor(mapping->x_sF*mapping->x_lo);
}

double enforceResolution(double in, GA_METHOD_Resolution res)
{
	double fracPart = 0.0;
	double intPart = 0.0;

	/* Decompose the data into integral and fractional parts */
	fracPart = std::modf(in, &intPart);

	/* Shift up, truncate, shift down */
	fracPart *= std::pow(10.0, (int)res);
	fracPart = floor(fracPart);
	fracPart /= std::pow(10.0, (int)res);

	return (intPart + fracPart);
}

void writeCSV_StepData(StepPerformance data, std::string filename)
{
	std::ofstream csvFile;
	csvFile.open(filename);
	size_t num_rows = data.data.rows();
	size_t num_cols = data.data.cols();

	//Write each full row
	for (int row = 0; row < num_rows; row++)
	{
		for (int col = 0; col < num_cols; col++)
			if (col == num_cols - 1)
				csvFile << std::to_string(data.data(row, col)*1.0);
			else
				csvFile << std::to_string(data.data(row, col)*1.0) << ",";

		csvFile << "\n";
	}

	//Write the performance data
	csvFile << data.finalValue_performance << "," << data.settlingTime_window << "\n";
	csvFile << data.delta_overshoot_performance << "\n";
	csvFile << data.percentOvershoot_performance << "\n";
	csvFile << data.riseTime_performance << "," << data.riseTime_Idx[0] << "," << data.riseTime_Idx[1] << "\n";
	csvFile << data.settlingTime_performance << "," << data.settlingTime_Idx << "\n";
	csvFile << data.steadyStateError_performance << "\n";

	csvFile.close();
}

//////////////////////////////////////////////////////////////////
/* CLASS: FCSOptimizer */
//////////////////////////////////////////////////////////////////
/*-----------------------------------------------
* Constructors/Destructor
*-----------------------------------------------*/
FCSOptimizer::FCSOptimizer()
{
	currentStatus = GA_IDLE;
	currentIteration = 0;
}

FCSOptimizer::~FCSOptimizer()
{
}

/*-----------------------------------------------
* Public Functions
*-----------------------------------------------*/
void FCSOptimizer::run()
{
	#if (FCS_TRACE_EXECUTION_TIME == 1)
	auto start = boost::chrono::high_resolution_clock::now();
	#endif

	std::cout << "Hello from thread " << boost::this_thread::get_id() << "!" << std::endl;
	
	//Do some initialization here 
	unsigned int commandPriority;
	message_queue::size_type received_size;
	int externCmd;


	for(;;)
	{
		/*-------------------------------------
		* Check for new commands from the main thread 
		*-------------------------------------*/
		if (commandQueue->try_receive(&externCmd, sizeof(externCmd), received_size, commandPriority))
		{
			std::cout << "I received command: " << externCmd << std::endl;

			/* Only set a flag here so that the error/message handing section below can take care of 
			commands from this block as well as commands issued by the core algorithm. */
			switch (externCmd)
			{
			case START:
				currentStatus = GA_OK;
				break;

			case PAUSE:
				currentStatus = GA_PAUSED;
				break;

			case STOP:
				currentStatus = GA_HALT;
				break;

			case REPORT_STATUS:
				//TODO: Enable passing back data to the main thread here
				//Or maybe do it in the error handling section??? doesn't quite fit either place...
				break;
			}
		}

		
		/*-------------------------------------
		* Propagate forward 1 generation
		*-------------------------------------*/
		if (currentStatus == GA_OK)
		{
			boundaryCheck();		/* Ensure the newly generated data falls within the user constraints */

			evaluateModel();		/* With the current population members, run a test scenario on the chosen model and record output data */

			evaluateFitness();		/* Use recorded performance data from last step to gauge how well it meets user goals */

			checkConvergence();		/* Given the new fitness scores, see if we have a "winner" that met user goals sufficiently */

			filterPopulation();		/* Apply random "filtering" to the population to stepResponse things like death/disaster/disease etc. */

			selectParents();		/* Of the current population, select those who will mate */

			breedGeneration();		/* Use the selected mating pairs from the last step to create new offspring */

			mutateGeneration();		/* Mutate the new generation's chromosomes based on random chance */
			
			currentIteration += 1;
		}

		/*-------------------------------------
		* Handle errors or external commands
		*-------------------------------------*/
		else
		{
			if (currentStatus == GA_HALT)
			{
				//Do any logging that might report why a halt was called. 
				#if (CONSOLE_LOGGING_ENABLED == 1)
				std::cout << "GA Solver was halted. Who knows why." << std::endl;
				#endif
				break;
			}

			if (currentStatus == GA_COMPLETE)
			{
				//Do other things?
				#if (CONSOLE_LOGGING_ENABLED == 1)
				std::cout << "GA Solver Complete. Exiting Optimizer." << std::endl;
				#endif
				break;
			}
		}

		/*-------------------------------------
		* Handle post-processing of round for reporting status to user
		*-------------------------------------*/
		//Do stuff

		/* Allow other waiting threads to run */
		boost::this_thread::yield();
	}

	#if (FCS_TRACE_EXECUTION_TIME == 1)
	auto stop = boost::chrono::high_resolution_clock::now();
	auto totalTime = boost::chrono::duration_cast<boost::chrono::milliseconds>(stop - start);

		#if (CONSOLE_LOGGING_ENABLED == 1)
		std::cout << "Total Execution Time: " << totalTime << std::endl;
		#endif
	#endif
}

void FCSOptimizer::init(FCSOptimizer_Init_t initializationSettings)
{
	settings = initializationSettings;

	/*-----------------------------
	* Order specific initialization sequence 
	*----------------------------*/
	initMemory();			/* Allocate container sizes before algorithm begins */
	initRNG();				/* Make sure the RNG is setup and well warmed up before use */
	initModel();			/* Call the specific model initializer */
	initPopulation();		/* Set up the initial population */


	//TODO: Set up the intelligent strategy selection engine here
	currentSolverParam = settings.solverParam;


	/* Create a new message queue for receiving commands from the main thread.
	This WILL destroy anything the main thread has created, so use this first and then
	create the main thread's interface after. */
	try
	{
		//Erase the previous message queue if it exists 
		message_queue::remove(settings.messageQueueName.data());

		//Create a fresh queue
		commandQueue = new message_queue(
			create_only,						/* Only create the queue */
			settings.messageQueueName.data(),	/* Grab the user specified queue name */
			10,									/* Limit the queue size to 10 */
			sizeof(int)							/* The queue will hold integers */
			);
	}
	catch (interprocess_exception &ex)
	{
		//TODO: Use the console mutex here to make sure the output is readable 
		std::cout << ex.what() << " in thread " << boost::this_thread::get_id() << std::endl;
		std::cout << "Unable to properly use the queue. Ignoring all messages from main thread." << std::endl;
	}

	/* After this function exits, the Genetic Algorithm should be ready to go */
	currentStatus = GA_OK;
}

void FCSOptimizer::requestOutput(FCSOptimizer_Output_t& output)
{
	//Might want to deprecate this? Not sure. Might just set the flag
	//here rather than checking for it in the main loop. Log the reference
	//to the output container so we can write to it later. 
}

/*-----------------------------------------------
* Private Functions
*-----------------------------------------------*/

void FCSOptimizer::initMemory()
{
	const size_t popSize = settings.advConvergenceParam.populationSize;
	population.resize(popSize);

	/* Glue memory between the GA functions */
	fcs_modelEvalutationData.resize(popSize);
	fcs_fitnessData.resize(popSize);
	fcs_parentSelections.resize(popSize);

	/* Allocate the memory needed to hold all the operation functions */
	runtimeStep.populationFilterInstances.resize(GA_POPULATION_TOTAL_OPTIONS);
	runtimeStep.evaluateModelInstances.resize(GA_MODEL_TOTAL_OPTIONS);
	runtimeStep.evaluateFitnessInstances.resize(GA_FITNESS_TOTAL_OPTIONS);
	runtimeStep.selectParentInstances.resize(GA_SELECT_TOTAL_OPTIONS);
	runtimeStep.breedInstances.resize(GA_BREED_TOTAL_OPTIONS);
	runtimeStep.mutateInstances.resize(GA_MUTATE_TOTAL_OPTIONS);
}

void FCSOptimizer::initRNG()
{
	/* Unfortunately, this is the best way I currently know how to get a single RNG interface
	to be used with many different kinds of possible engines and distributions. I tried templating
	the whole thing, but it ended up being far more trouble than it was worth. This is the solution
	found in the available time frame.
	*/

	std::cout << "Initializing the RNG Engines with time based seed..." << std::endl;

	if (settings.solverParam.rngDistribution == GA_DISTRIBUTION_UNIFORM_REAL)
	{
		auto Kp_distribution = urd_t(
			settings.pidControlSettings.tuningLimits.Kp.lower,
			settings.pidControlSettings.tuningLimits.Kp.upper);

		auto Ki_distribution = urd_t(
			settings.pidControlSettings.tuningLimits.Ki.lower,
			settings.pidControlSettings.tuningLimits.Ki.upper);

		auto Kd_distribution = urd_t(
			settings.pidControlSettings.tuningLimits.Kd.lower,
			settings.pidControlSettings.tuningLimits.Kd.upper);

		switch (settings.solverParam.rngEngine)
		{
		case GA_MERSENNE_TWISTER:
			PID_RNG.Kp = boost::make_shared<RNGInstance<boost::mt19937, urd_t>>(Kp_distribution);
			PID_RNG.Ki = boost::make_shared<RNGInstance<boost::mt19937, urd_t>>(Ki_distribution);
			PID_RNG.Kd = boost::make_shared<RNGInstance<boost::mt19937, urd_t>>(Kd_distribution);
			break;

		/*Add support for more as needed*/
		default: break;
		}
	}

	if (settings.solverParam.rngDistribution == GA_DISTRIBUTION_UNIFORM_INT)
	{
		auto Kp_distribution = uid_t(
			(int)settings.pidControlSettings.tuningLimits.Kp.lower,
			(int)settings.pidControlSettings.tuningLimits.Kp.upper);

		auto Ki_distribution = uid_t(
			(int)settings.pidControlSettings.tuningLimits.Ki.lower,
			(int)settings.pidControlSettings.tuningLimits.Ki.upper);

		auto Kd_distribution = uid_t(
			(int)settings.pidControlSettings.tuningLimits.Kd.lower,
			(int)settings.pidControlSettings.tuningLimits.Kd.upper);

		switch (settings.solverParam.rngEngine)
		{
		case GA_MERSENNE_TWISTER:
			PID_RNG.Kp = boost::make_shared<RNGInstance<boost::mt19937, uid_t>>(Kp_distribution);
			PID_RNG.Ki = boost::make_shared<RNGInstance<boost::mt19937, uid_t>>(Ki_distribution);
			PID_RNG.Kd = boost::make_shared<RNGInstance<boost::mt19937, uid_t>>(Kd_distribution);
			break;

			/*Add support for more as needed*/
		default: break;
		}
	}

	std::cout << "Done!" << std::endl;
}

void FCSOptimizer::initModel()
{
	/* State Space Model */
	

	/* Neural Network Model */
}

void FCSOptimizer::initPopulation()
{
	/* Initialize the mapping constants to enable conversions between real
	valued data and unsigned integers (bit manipulation) */
	calculateMappingCoefficients(
		&mapCoefficients_Kp, 
		settings.pidControlSettings.tuningLimits.Kp.lower, 
		settings.pidControlSettings.tuningLimits.Kp.upper);

	calculateMappingCoefficients(
		&mapCoefficients_Ki, 
		settings.pidControlSettings.tuningLimits.Ki.lower,
		settings.pidControlSettings.tuningLimits.Ki.upper);

	calculateMappingCoefficients(
		&mapCoefficients_Kd,
		settings.pidControlSettings.tuningLimits.Kd.lower,
		settings.pidControlSettings.tuningLimits.Kd.upper);

	/* Initialize the PID values for every population member */
	PID_RNG.Kp->acquireEngine();
	PID_RNG.Ki->acquireEngine();
	PID_RNG.Kd->acquireEngine();

	for (size_t i = 0; i < settings.advConvergenceParam.populationSize; i++)
	{
		/* Get some random PID values to start off with */
		population[i].realPID.Kp = PID_RNG.Kp->getDouble();
		population[i].realPID.Ki = PID_RNG.Ki->getDouble();
		population[i].realPID.Kd = PID_RNG.Kd->getDouble();
	}

	PID_RNG.Kp->releaseEngine();
	PID_RNG.Ki->releaseEngine();
	PID_RNG.Kd->releaseEngine();
}


void FCSOptimizer::evaluateModel()
{
	const GA_METHOD_ModelEvaluation modelType = currentSolverParam.modelType;

	/* Ensure an instance of the mutation type exists before evaluation */
	if (!runtimeStep.evaluateModelInstances[modelType])
	{
		switch (modelType)
		{
		case GA_MODEL_STATE_SPACE:
			//TODO: refactor this ridiculous "runtimeStep"...I forgot what it meant and
			// didn't understand at a glance
			runtimeStep.evaluateModelInstances[modelType] = boost::make_shared<StateSpaceEvaluator>();
			break;

		case GA_MODEL_NEURAL_NETWORK:
			runtimeStep.evaluateModelInstances[modelType] = boost::make_shared<NeuralNetworkEvaluator>();
			break;

		//Add more as needed here
		default:
			std::cout << "Evaluation model not configured or is unknown. You are about to crash." << std::endl;
			break;
		}
	}

	if (modelType == GA_MODEL_STATE_SPACE)
	{
		StateSpaceModelInput input;
		StateSpaceModelOutput output;

		/* Simulation time constraints */
		input.dt = 0.01;			//TODO: Replace with un-hardcoded version
		input.startTime = 0.0;
		input.endTime = 1.5;

		/* Simulation Model */
		input.simulationType = STEP;
		input.model = settings.stateSpaceModel;

		/* This is incredibly slow (~6mS, eT3.0, dT0.1) per member...ridiculous. More than likely
		this is due to being single threaded and not taking up enough logical cores. */
		//TODO: Enforce machine specific multi threading based on number of logical cores 
		for (int member = 0; member < settings.advConvergenceParam.populationSize; member++)
		{
			/* FUTURE NOTE:
			When doing the multi threaded version, it probably would be extremely helpful to build 
			up all the data needed to run the simulation and pass it into a thread that is waiting 
			to consume some data and spit out results before sleeping again.
			
			Create these threads during initialization so that thread allocation and destruction does
			not consume any runtime resources. There will be as many threads as the number of logical 
			processor cores. */

			/* Swap out the PID values and simulate */
			input.pid = population[member].realPID;
			//input.pid = { 12.05, 68.38, 0.0 };

			runtimeStep.evaluateModelInstances[modelType]->evaluate(input, output);

			/* Pass on the resulting data to the optimizer's record of everything */
			fcs_modelEvalutationData[member].modelType = GA_MODEL_STATE_SPACE;
			fcs_modelEvalutationData[member].ss_output = output;
		}
	}
	
	else if (modelType == GA_MODEL_NEURAL_NETWORK)
	{
		NeuralNetworkModelInput input;
		NeuralNetworkModelOutput output;

		runtimeStep.evaluateModelInstances[modelType]->evaluate(input, output);
	} 
}

void FCSOptimizer::evaluateFitness()
{
	const GA_METHOD_ModelEvaluation modelType = currentSolverParam.modelType;
	const GA_METHOD_FitnessEvaluation fitnessType = currentSolverParam.fitnessType;

	GA_EvaluateFitnessDataInput input;
	GA_EvaluateFitnessDataOutput output;

	/* Assign the static fitness evaluation parameters */
	input.goals = settings.pidControlSettings.performanceGoals;
	input.tolerance = settings.pidControlSettings.performanceTolerance;

	/*-----------------------------
	* Ensure a fitness evaluation instance exists 
	*----------------------------*/
	if (!runtimeStep.evaluateFitnessInstances[fitnessType])
	{
		switch (fitnessType)
		{
		case GA_FITNESS_WEIGHTED_SUM:
			runtimeStep.evaluateFitnessInstances[fitnessType] = boost::make_shared<WeightedSum>();
			break;

		case GA_FITNESS_NON_DOMINATED_SORT:
			runtimeStep.evaluateFitnessInstances[fitnessType] = boost::make_shared<NonDominatedSort>();
			break;

		//Add more as needed here
		default:
			std::cout << "Fitness method not configured or is unknown. You are about to crash." << std::endl;
			break;
		}
	}

	/*-----------------------------
	* Calculate the fitness 
	*----------------------------*/
	if (fitnessType == GA_FITNESS_WEIGHTED_SUM)
	{
		if (modelType == GA_MODEL_STATE_SPACE)
		{
			for (int member = 0; member < settings.advConvergenceParam.populationSize; member++)
			{
				/* Update the simulation results to be passed in */
				input.POS = fcs_modelEvalutationData[member].ss_output.stepPerformance.percentOvershoot_performance;
				input.SSER = fcs_modelEvalutationData[member].ss_output.stepPerformance.steadyStateError_performance;
				input.TR = fcs_modelEvalutationData[member].ss_output.stepPerformance.percentOvershoot_performance;
				input.TS = fcs_modelEvalutationData[member].ss_output.stepPerformance.settlingTime_performance;
				input.FV = fcs_modelEvalutationData[member].ss_output.stepPerformance.finalValue_performance;

				input.simulationData = fcs_modelEvalutationData[member].ss_output.stepPerformance.data;
				//NOTE: Eventually I am going to need to specify what kind of simulation was performed...

				/* Evaluate fitness */
				runtimeStep.evaluateFitnessInstances[fitnessType]->evaluateFitness(input, output);

				/* Copy the results of the step performance analyzer into the output struct */
				output.fit.stepPerformanceData = fcs_modelEvalutationData[member].ss_output.stepPerformance;

				/* Pass on the resulting data to the optimizer's record of everything */
				fcs_fitnessData[member].modelType = modelType;
				fcs_fitnessData[member].fit = output.fit;

				#if (DEBUGGING_ENABLED == 1)
				double memberFit = 0.0;
				memberFit = output.fit.global_fitness;
				memberFit += 1.0;
				#endif

				#if (FILE_LOGGING_ENABLED == 1) && (GA_FILELOG_MEMBER_FIT_DATA == 1)
				std::string filename = "Member" + std::to_string(member) + "_StepPerformanceFitnessData.csv";
				writeCSV_StepData(output.fit.stepPerformanceData, filename);
				#endif
			}
		}
		else if (modelType == GA_MODEL_NEURAL_NETWORK)
		{

		}
	}

	else if (fitnessType == GA_FITNESS_NON_DOMINATED_SORT)
	{

	}
}

void FCSOptimizer::checkConvergence()
{
	/*-----------------------------------------------
	* Save best result from each round
	*-----------------------------------------------*/
	double iteration_bestFitScore = 0.0;			/* Keep track of highest fitness score found */
	FCSOptimizer_FitnessData iteration_bestFit;		/* Copy of associated best fitness results */

	for (unsigned int member = 0; member < settings.advConvergenceParam.populationSize; member++)
	{
		if (fcs_fitnessData[member].fit.global_fitness > iteration_bestFitScore)
		{
			iteration_bestFitScore = fcs_fitnessData[member].fit.global_fitness;
			iteration_bestFit = fcs_fitnessData[member];
		}
	}

	fcs_generationalBestFitnessData.push_back(iteration_bestFit);

	#if (CONSOLE_LOGGING_ENABLED == 1)
	//print_to_console_mutex->lock();
	//TODO: Add which generation this is out of total, ie (1 of 40)
	std::cout << "Iteration Best Fit: " << iteration_bestFitScore << std::endl;
	//print_to_console_mutex->unlock();
	#endif

	/*-----------------------------------------------
	* Break based on generation limits
	*-----------------------------------------------*/
	if (currentIteration >= settings.advConvergenceParam.generationLimit)
	{
		#ifdef _DEBUG
		if (settings.advConvergenceParam.generationLimit == 0)
			std::cout << "You forgot to specify a generational limit for the algorithm." << std::endl;
		#endif

		#if (CONSOLE_LOGGING_ENABLED == 1)
		//print_to_console_mutex->lock();
		std::cout << "Generational Limit Reached! Done." << std::endl;
		//print_to_console_mutex->unlock();
		#endif
		currentStatus = GA_COMPLETE;
	}

	/*-----------------------------------------------
	* Break based on finding a good solution 
	*-----------------------------------------------*/
	if (iteration_bestFitScore > 0.95)
	{
		#if (CONSOLE_LOGGING_ENABLED == 1)
		//print_to_console_mutex->lock();
		std::cout << "\nFound a good enough solution. Done." << std::endl;
		//print_to_console_mutex->unlock();
		#endif

		currentStatus = GA_COMPLETE;
	}
		
}

void FCSOptimizer::filterPopulation()
{
	const GA_METHOD_PopulationFilter filterType = currentSolverParam.filterType;

	GA_PopulationFilterDataInput input;
	GA_PopulationFilterDataOutput output;

	/*-----------------------------
	* Ensure a filtering evaluation instance exists
	*----------------------------*/
	if (!runtimeStep.populationFilterInstances[filterType])
	{
		switch (filterType)
		{
		case GA_POPULATION_STATIC_FILTER:
			runtimeStep.populationFilterInstances[filterType] = boost::make_shared<StaticFilter>(
				settings.pidControlSettings.tuningLimits.Kp.upper,
				settings.pidControlSettings.tuningLimits.Kp.lower,
				settings.pidControlSettings.tuningLimits.Ki.upper,
				settings.pidControlSettings.tuningLimits.Ki.lower,
				settings.pidControlSettings.tuningLimits.Kd.upper,
				settings.pidControlSettings.tuningLimits.Kd.lower);
			break;

		case GA_POPULATION_DYNAMIC_FILTER:
			runtimeStep.populationFilterInstances[filterType] = boost::make_shared<DynamicFilter>();
			break;

			//Add more as needed here
		default:
			std::cout << "Filtering method not configured or is unknown. You are about to crash." << std::endl;
			break;
		}
	}

	/*-----------------------------
	* Perform the filtering operation
	*----------------------------*/
	
	/* Populate the input struct */
	std::transform(fcs_fitnessData.begin(), fcs_fitnessData.end(), std::back_inserter(input.currentGlobalFitScores),
		[](const FCSOptimizer_FitnessData& input) { return input.fit.global_fitness; });

	input.static_performanceThreshold = 0.4;

	runtimeStep.populationFilterInstances[filterType]->filter(input, output);


	/* Replace the indicated population members */
	for (int member = 0; member < output.replacedMemberIndexes.size(); member++)
		population[output.replacedMemberIndexes[member]].realPID = output.replacementPIDValues[member];
}

void FCSOptimizer::selectParents()
{
	const GA_METHOD_ParentSelection selectType = currentSolverParam.selectType;

	GA_SelectParentDataInput input; 
	GA_SelectParentDataOutput output; 

	input.populationSize = settings.advConvergenceParam.populationSize;
	output.parentSelections.resize(settings.advConvergenceParam.populationSize);

	/* Ensure an instance of the selection type exists before calling the selection function */
	if (!runtimeStep.selectParentInstances[selectType])
	{
		switch (selectType)
		{
		case GA_SELECT_RANDOM:
			runtimeStep.selectParentInstances[selectType] = boost::make_shared<RandomSelection>(settings.advConvergenceParam.populationSize);
			break;

		case GA_SELECT_RANKED:
			runtimeStep.selectParentInstances[selectType] = boost::make_shared<RankedSelection>();
			break;

		case GA_SELECT_ROULETTE:
			runtimeStep.selectParentInstances[selectType] = boost::make_shared<RouletteSelection>();
			break;

		case GA_SELECT_STOCHASTIC_SAMPLING:
			runtimeStep.selectParentInstances[selectType] = boost::make_shared<StochasticSelection>();
			break;

		case GA_SELECT_TOURNAMENT:
			runtimeStep.selectParentInstances[selectType] = boost::make_shared<TournamentSelection>();
			break;

		case GA_SELECT_ELITIST:
			runtimeStep.selectParentInstances[selectType] = boost::make_shared<ElitistSelection>();
			break;

		//Add more as needed here
		default:
			std::cout << "Parent selection method not configured or is unknown. You are about to crash." << std::endl;
			break;
		}
	}

	/* Fill in the input data structure based on what kind of selectType is used */
	if (selectType == GA_SELECT_RANDOM)
	{
	}
	else if (selectType == GA_SELECT_RANKED)
	{
	}
	else if (selectType == GA_SELECT_ROULETTE)
	{
	}
	else if (selectType == GA_SELECT_STOCHASTIC_SAMPLING)
	{
	}
	else if (selectType == GA_SELECT_TOURNAMENT)
	{
		std::transform(fcs_fitnessData.begin(), fcs_fitnessData.end(), std::back_inserter(input.popGlobalFitScores),
			[](const FCSOptimizer_FitnessData& input) { return input.fit.global_fitness; });
	}
	else if (selectType == GA_SELECT_ELITIST)
	{
	}
	else
	{
		std::cout << "You have not chosen a valid Parent Selection strategy" << std::endl;
	}

	
	runtimeStep.selectParentInstances[selectType]->selectParent(input, output);

	
	fcs_parentSelections = output.parentSelections;
}

void FCSOptimizer::breedGeneration()
{
	const GA_METHOD_Breed breedType = currentSolverParam.breedType;

	GA_BreedingDataInput input;
	GA_BreedingDataOutput output;

	/* Ensure an instance of the breeding type exists before calling the breed function */
	if (!runtimeStep.breedInstances[breedType])
	{
		switch (breedType)
		{
		case GA_BREED_SIMPLE_CROSSOVER:
			runtimeStep.breedInstances[breedType] = boost::make_shared<SimpleCrossover>();
			break;

		case GA_BREED_DYNAMIC_CROSSOVER:
			runtimeStep.breedInstances[breedType] = boost::make_shared<DynamicCrossover>();
			break;

		case GA_BREED_FIXED_POINT_CROSSOVER:
			runtimeStep.breedInstances[breedType] = boost::make_shared<FixedPointCrossover>();
			break;

		case GA_BREED_SIMULATED_BINARY_CROSSOVER:
			runtimeStep.breedInstances[breedType] = boost::make_shared<SimulatedBinaryCrossover>();
			break;

		//Add more as needed here
		default: 
			std::cout << "Breeding method not configured or is unknown. You are about to crash." << std::endl;
			break;
		}
	}


	/*-----------------------------------------------
	* Assign data to the input struct. For now, force
	* chromosomes to be real valued.
	*-----------------------------------------------*/
	input.parentSelections = fcs_parentSelections;

	if (currentSolverParam.chromType == MAPPING_TYPE_REAL)
	{
		input.chromType = currentSolverParam.chromType;

		//TODO: I need the user to specify these variable values in some kind of advanced initializer....
		input.crossoverPoint = 9;
		input.swap_both_chrom_halves = false;
		input.swap_lower_chrom_half = false;

		/* Fancy iterator to copy over the relevant genetic material. TODO: Switch out population real/mapped types into chrom type */
		std::transform(population.begin(), population.end(), std::back_inserter(input.d_chrom),
			[](const FCSOptimizer_PopulationMember& member)
		{
			return GA_PIDChromosome<double>() = { member.realPID.Kp, member.realPID.Ki, member.realPID.Kd };
		}
		);
	}

	else if (currentSolverParam.chromType == MAPPING_TYPE_BIT_FIELD)
	{
		std::cout << "Global bit field mapping type not supported yet." << std::endl;
	}
	

	input.mapCoeff_Kp = &mapCoefficients_Kp;
	input.mapCoeff_Ki = &mapCoefficients_Ki;
	input.mapCoeff_Kd = &mapCoefficients_Kd;

	runtimeStep.breedInstances[breedType]->breed(input, output);

	fcs_bredChromosomes = output;
}

void FCSOptimizer::mutateGeneration()
{
	const GA_METHOD_MutateType mutateType = currentSolverParam.mutateType;

	GA_MutateDataInput input;
	GA_MutateDataOutput output;

	/* Ensure an instance of the mutation type exists before calling the mutate function */
	if (!runtimeStep.mutateInstances[mutateType])
	{
		switch (mutateType)
		{
		case GA_MUTATE_BIT_FLIP:
			runtimeStep.mutateInstances[mutateType] = boost::make_shared<BitFlipMutator>();
			break;

		case GA_MUTATE_ADD_SUB:
			runtimeStep.mutateInstances[mutateType] = boost::make_shared<AddSubMutator>();
			break;

		//Add more as needed here
		default:
			std::cout << "Mutate method not configured or is unknown. You are about to crash." << std::endl;
			break;
		}
	}


	input.mutateProbType = currentSolverParam.mutateProbabilityType;
	input.optimizerChromType = currentSolverParam.chromType;			/* Final desired chromosome type */
	input.chromType = fcs_bredChromosomes.chromType;					/* Current input data chromosome type*/

	//TODO: I need the user to specify these parameters!!!
	input.mutationProbabilityThreshold = 0.25;

	if (input.chromType == MAPPING_TYPE_REAL)
		input.d_chrom = fcs_bredChromosomes.d_chrom;

	else if (input.chromType == MAPPING_TYPE_BIT_FIELD)
		input.u16_chrom = fcs_bredChromosomes.u16_chrom;

	input.mapCoeff_Kp = &mapCoefficients_Kp;
	input.mapCoeff_Ki = &mapCoefficients_Ki;
	input.mapCoeff_Kd = &mapCoefficients_Kd;


	runtimeStep.mutateInstances[mutateType]->mutate(input, output);


	//Put the mutated data back into the population
	if (currentSolverParam.chromType == MAPPING_TYPE_REAL && output.chromType == MAPPING_TYPE_REAL)
	{
		for (int i = 0; i < population.size(); i++)
		{
			population[i].realPID.Kp = output.d_chrom[i].Kp;
			population[i].realPID.Ki = output.d_chrom[i].Ki;
			population[i].realPID.Kd = output.d_chrom[i].Kd;
		}
	}
	else if (currentSolverParam.chromType == MAPPING_TYPE_BIT_FIELD && output.chromType == MAPPING_TYPE_BIT_FIELD)
	{
		for (int i = 0; i < population.size(); i++)
		{
			population[i].mappedPID.Kp = output.u16_chrom[i].Kp;
			population[i].mappedPID.Ki = output.u16_chrom[i].Ki;
			population[i].mappedPID.Kd = output.u16_chrom[i].Kd;
		}
	}
}

void FCSOptimizer::boundaryCheck()
{
	/* Force decimal point resolution so we aren't searching through an ENORMOUS space */
	enforceResolution();

	//Do other things
}


void FCSOptimizer::enforceResolution()
{
	double fracPart = 0.0;
	double intPart = 0.0;

	for (int member = 0; member < settings.advConvergenceParam.populationSize; member++)
	{
		/* Decompose the data into integral and fractional parts */
		fracPart = std::modf(population[member].realPID.Kp, &intPart);

		/* Shift up, truncate, shift down data */
		fracPart *= std::pow(10.0, (int)currentSolverParam.resolutionType);
		fracPart = floor(fracPart);
		fracPart /= std::pow(10.0, (int)currentSolverParam.resolutionType);

		/* Reassign truncated value */
		population[member].realPID.Kp = intPart + fracPart;

		/*KI*/
		fracPart = std::modf(population[member].realPID.Ki, &intPart);

		fracPart *= std::pow(10.0, (int)currentSolverParam.resolutionType);
		fracPart = floor(fracPart);
		fracPart /= std::pow(10.0, (int)currentSolverParam.resolutionType);

		population[member].realPID.Ki = intPart + fracPart;

		/*KD*/
		fracPart = std::modf(population[member].realPID.Kd, &intPart);

		fracPart *= std::pow(10.0, (int)currentSolverParam.resolutionType);
		fracPart = floor(fracPart);
		fracPart /= std::pow(10.0, (int)currentSolverParam.resolutionType);

		population[member].realPID.Kd = intPart + fracPart;
	}
}

// int FCSOptimizer::reportResults(int trialNum)
// {
// 	/* Find the highest performer */
// 	double highestPerformerVal = 0.0;
// 	int highestPerformerIndex = 0;
// 	for (int i = 0; i < GA_BestFitnessValues.size(); i++)
// 	{
// 		if (GA_BestFitnessValues.data()[i].global_fitness > highestPerformerVal)
// 		{
// 			highestPerformerVal = GA_BestFitnessValues.data()[i].global_fitness;
// 			highestPerformerIndex = i;
// 		}
// 	}
// 
// 	/*-----------------------------------------------
// 	* Write results to output screen
// 	*-----------------------------------------------*/
// 	printResultHighlights(highestPerformerVal, highestPerformerIndex);
// 
// 	/*-----------------------------------------------
// 	* Write results to file
// 	*-----------------------------------------------*/
// 	#if defined(GA_CPU_SINGLE_THREADED) && !defined(GA_CPU_MULTI_THREADED)
// 	std::string filename = "results/singlethreaded/Trial" + std::to_string(trialNum) + "_" + ga_instance_optimizer_name + "_bestPerformer.csv";
// 	#endif
// 
// 	#if defined(GA_CPU_MULTI_THREADED) && !defined(GA_CPU_SINGLE_THREADED)
// 	std::string filename = "results/multithreaded/Trial" + std::to_string(trialNum) + "_" + ga_instance_optimizer_name + "_bestPerformer.csv";
// 	#endif
// 
// 	writeCSV_StepData(GA_BestFitnessValues.data()[highestPerformerIndex].fitness_performance, filename);
// 
// 	return highestPerformerIndex;
// }

// void FCSOptimizer::printResultHighlights(double best_fit, int best_fit_idx)
// {
// 	print_to_console_mutex->lock();
// 	/*-----------------------------------------------
// 	* Best Fit
// 	*-----------------------------------------------*/
// 	std::cout << "=============================================================\n"
// 		<< "\t\tPERFORMANCE RESULTS: " + ga_instance_optimizer_name + "\n"
// 		<< "============================================================="
// 		<< std::endl;
// 	std::cout << "Best Fit:\t\t" << best_fit << std::endl;
// 
// 	/*-----------------------------------------------
// 	* Final Value
// 	*-----------------------------------------------*/
// 	std::cout << "\n" << std::endl;
// 	std::cout << "Final Value:\t\t"
// 		<< GA_BestFitnessValues.data()[best_fit_idx].fitness_performance.finalValue_performance
// 		<< std::endl;
// 
// 	/*-----------------------------------------------
// 	* Percent Overshoot
// 	*-----------------------------------------------*/
// 	std::cout << "\n" << std::endl;
// 	std::cout << "POS Fitness:\t\t"
// 		<< GA_BestFitnessValues.data()[best_fit_idx].percentOvershoot_fitness
// 		<< std::endl;
// 	std::cout << "POS Actual:\t\t"
// 		<< GA_BestFitnessValues.data()[best_fit_idx].fitness_performance.percentOvershoot_performance
// 		<< "%"
// 		<< std::endl;
// 	std::cout << "POS Delta:\t\t"
// 		<< GA_BestFitnessValues.data()[best_fit_idx].fitness_performance.delta_overshoot_performance
// 		<< ""
// 		<< std::endl;
// 	std::cout << "POS Target:\t\t"
// 		<< ga_instance_pid_config_data->performanceGoals.percentOvershoot_goal*100.0
// 		<< "%"
// 		<< std::endl;
// 
// 	/*-----------------------------------------------
// 	* Rise Time
// 	*-----------------------------------------------*/
// 	std::cout << "\n" << std::endl;
// 	std::cout << "RiseTime Fitness:\t"
// 		<< GA_BestFitnessValues.data()[best_fit_idx].riseTime_fitness
// 		<< std::endl;
// 	std::cout << "RiseTime Actual:\t"
// 		<< GA_BestFitnessValues.data()[best_fit_idx].fitness_performance.riseTime_performance
// 		<< " sec"
// 		<< std::endl;
// 	std::cout << "RiseTime Target:\t"
// 		<< ga_instance_pid_config_data->performanceGoals.riseTime_goal
// 		<< " sec"
// 		<< std::endl;
// 
// 	/*-----------------------------------------------
// 	* Settling Time
// 	*-----------------------------------------------*/
// 	std::cout << "\n" << std::endl;
// 	std::cout << "SettleTime Fitness:\t"
// 		<< GA_BestFitnessValues.data()[best_fit_idx].settlingTime_fitness
// 		<< std::endl;
// 	std::cout << "SettleTime Actual:\t"
// 		<< GA_BestFitnessValues.data()[best_fit_idx].fitness_performance.settlingTime_performance
// 		<< " sec"
// 		<< std::endl;
// 	std::cout << "SettleTime Target:\t"
// 		<< ga_instance_pid_config_data->performanceGoals.settlingTime_goal
// 		<< " sec"
// 		<< std::endl;
// 
// 	/*-----------------------------------------------
// 	* Steady State Error
// 	*-----------------------------------------------*/
// 	std::cout << "\n" << std::endl;
// 	std::cout << "SSErr Fitness:\t\t"
// 		<< GA_BestFitnessValues.data()[best_fit_idx].steadyStateError_fitness
// 		<< std::endl;
// 	std::cout << "SSErr Actual:\t\t"
// 		<< GA_BestFitnessValues.data()[best_fit_idx].fitness_performance.steadyStateError_performance
// 		<< ""
// 		<< std::endl;
// 	std::cout << "SSErr Target:\t\t"
// 		<< "+/- "
// 		<< ga_instance_pid_config_data->performanceGoals.steadyStateError_goal
// 		<< ""
// 		<< std::endl;
// 
// 	/*-----------------------------------------------
// 	* PID Solution
// 	*-----------------------------------------------*/
// 	std::cout << "\n" << std::endl;
// 	std::cout << "Kp:\t\t\t" << GA_BestFitnessValues.data()[best_fit_idx].fitness_performance.Kp << std::endl;
// 	std::cout << "Ki:\t\t\t" << GA_BestFitnessValues.data()[best_fit_idx].fitness_performance.Ki << std::endl;
// 	std::cout << "Kd:\t\t\t" << GA_BestFitnessValues.data()[best_fit_idx].fitness_performance.Kd << std::endl;
// 
// 	print_to_console_mutex->unlock();
// }