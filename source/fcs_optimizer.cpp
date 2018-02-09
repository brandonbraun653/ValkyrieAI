#include "fcs_optimizer.h"
#include <string>
#include <iostream>
#include <fstream>

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
	csvFile << data.steadyStateValue_performance << "," << data.settlingPcntRange << "\n";
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
	currentGeneration = 0;
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

	refreshCounter = settings.advConvergenceParam.iterations_before_refresh;


	/* Generate the first population */
	evaluateModel(parents);
	evaluateFitness(parents);
	parents = sortPopulation(&parents, NULL);

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
				//Or maybe do it in the error handling section??? doesn't quite fitness either place...
				break;
			}
		}

		
		/*-------------------------------------
		* Propagate forward 1 generation
		*-------------------------------------*/
		if (currentStatus == GA_OK)
		{
			selectParents(parents);								/* Of the current population, select those who will mate */

			breedGeneration(parents);							/* Breed the parents and stores the output as bit-mapped chromosomes that
																are ready for the mutation stage.*/

			mutateGeneration(children);							/* Mutate the chromosomes resulting from the breedGeneration step. Output the results
																into the children population struct. */

			boundaryCheck(children);							/* Ensure the new children chromosomes fall within the user constraints */

			evaluateModel(children);							/* Evaluate the children in the system model */

			evaluateFitness(children);							/* Use recorded performance data from last step to gauge how well it meets user goals */

			parents = sortPopulation(&parents, &children);		/* Given the two populations (parent and child), sort them for the best performers 
																and select new parents from the top N solutions, where N = user defined population size. */

			checkConvergence();									/* Given the new sorted population, see if we have a "winner" that met user goals sufficiently */
			

			/* Once the refresh counter reaches zero, the current state of the algorithm will be 
			analyzed to determine if the solution is converging decently well. Assuming progress
			has slowed down, a different combination of GA operators (mutation/breeding/etc) will
			be selected at random and the algorithm will proceed as normal until then next refresh.*/
			if (--refreshCounter == 0)
			{
				if (algorithmNeedsNewSteps())
					updateAlgorithmSteps();
				 
				refreshCounter = settings.advConvergenceParam.iterations_before_refresh;
				std::cout << "I hit the refresh target!" << std::endl;
			}
		} else {

		/*-------------------------------------
		* Handle errors or external commands
		*-------------------------------------*/
			if (currentStatus == GA_HALT)
			{
				#if (CONSOLE_LOGGING_ENABLED == 1)
				std::cout << "GA Solver was halted. Who knows why." << std::endl;
				#endif
				break;
			}

			if (currentStatus == GA_COMPLETE)
			{
				#if (CONSOLE_LOGGING_ENABLED == 1)
				std::cout << "GA Solver Complete. Exiting Optimizer." << std::endl;
				#endif
				break;
			}
		}

		/*-------------------------------------
		* Do post-processing for reporting status to user
		*-------------------------------------*/
		gatherStatisticalData();

		/* Allow other waiting threads to run */
		currentGeneration += 1;
		boost::this_thread::yield();
	}

	logData();

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
	initStatistics();

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
	const size_t genLimit = settings.advConvergenceParam.generationLimit;
	parents.resize(popSize);
	children.resize(popSize);

	/* Glue memory between the GA functions */
	fcs_parentSelections.resize(popSize);

	/* Allocate the memory needed to hold all the operation functions */
	runtimeStep.populationFilterInstances.resize(GA_POPULATION_TOTAL_OPTIONS);
	runtimeStep.evaluateModelInstances.resize(GA_MODEL_TOTAL_OPTIONS);
	runtimeStep.evaluateFitnessInstances.resize(GA_FITNESS_TOTAL_OPTIONS);
	runtimeStep.selectParentInstances.resize(GA_SELECT_TOTAL_OPTIONS);
	runtimeStep.breedInstances.resize(GA_BREED_TOTAL_OPTIONS);
	runtimeStep.mutateInstances.resize(GA_MUTATE_TOTAL_OPTIONS);
	runtimeStep.sortingInstances.resize(GA_SORT_TOTAL_OPTIONS);

	/* Statistics Containers */
	GenerationalChromStats.Kp.resize(genLimit+1);
	GenerationalChromStats.Ki.resize(genLimit+1);
	GenerationalChromStats.Kd.resize(genLimit+1);

	ChromOccurance.Kp.resize((int)(settings.pidControlSettings.tuningLimits.Kp.upper - settings.pidControlSettings.tuningLimits.Kp.lower)+1);
	ChromOccurance.Ki.resize((int)(settings.pidControlSettings.tuningLimits.Ki.upper - settings.pidControlSettings.tuningLimits.Ki.lower)+1);
	ChromOccurance.Kd.resize((int)(settings.pidControlSettings.tuningLimits.Kd.upper - settings.pidControlSettings.tuningLimits.Kd.lower)+1);
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
		parents[i].dChrom.Kp = PID_RNG.Kp->getDouble();
		parents[i].dChrom.Ki = PID_RNG.Ki->getDouble();
		parents[i].dChrom.Kd = PID_RNG.Kd->getDouble();
	}

	PID_RNG.Kp->releaseEngine();
	PID_RNG.Ki->releaseEngine();
	PID_RNG.Kd->releaseEngine();
}

void FCSOptimizer::initStatistics()
{
	/* Initialize the count of each possible PID term to zero */
	std::for_each(ChromOccurance.Kp.begin(), ChromOccurance.Kp.end(), [](int& val) {val = 0; });
	std::for_each(ChromOccurance.Ki.begin(), ChromOccurance.Ki.end(), [](int& val) {val = 0; });
	std::for_each(ChromOccurance.Kd.begin(), ChromOccurance.Kd.end(), [](int& val) {val = 0; });

	/* Calculate the ideal variance of the chromosomes given a uniform distribution (hard-coded) 
	See: https://en.wikibooks.org/wiki/Statistics/Distributions/Uniform */

	IdealChromVariance.Kp = pow((settings.pidControlSettings.tuningLimits.Kp.upper - settings.pidControlSettings.tuningLimits.Kp.lower), 2) / 12.0;
	IdealChromVariance.Ki = pow((settings.pidControlSettings.tuningLimits.Ki.upper - settings.pidControlSettings.tuningLimits.Ki.lower), 2) / 12.0;
	IdealChromVariance.Kd = pow((settings.pidControlSettings.tuningLimits.Kd.upper - settings.pidControlSettings.tuningLimits.Kd.lower), 2) / 12.0;

}

void FCSOptimizer::evaluateModel(PopulationType& population)
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
			input.pid = population[member].dChrom;
			//input.pid.Kp = 12.0;
			//input.pid.Ki = 68.0;
			//input.pid.Kd = 0.0;

			runtimeStep.evaluateModelInstances[modelType]->evaluate(input, output);

			/* Pass on the resulting data to the optimizer's record of everything */
			population[member].modelType = modelType;
			population[member].evaluationPerformance.stepPerformanceData = output.stepPerformance;
		}
	}
	
	else if (modelType == GA_MODEL_NEURAL_NETWORK)
	{
	} 
}

void FCSOptimizer::evaluateFitness(PopulationType& population)
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

		case GA_FITNESS_MEAN_SQUARE_ERROR:
			runtimeStep.evaluateFitnessInstances[fitnessType] = boost::make_shared<MeanSquareError>();
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
		for (int member = 0; member < settings.advConvergenceParam.populationSize; member++)
		{
			/* Update the simulation results to be passed in */
			input.POS	= population[member].evaluationPerformance.stepPerformanceData->percentOvershoot_performance;
			input.SSER	= population[member].evaluationPerformance.stepPerformanceData->steadyStateError_performance;
			input.TR	= population[member].evaluationPerformance.stepPerformanceData->percentOvershoot_performance;
			input.TS	= population[member].evaluationPerformance.stepPerformanceData->settlingTime_performance;
			input.FV	= population[member].evaluationPerformance.stepPerformanceData->steadyStateValue_performance;

			input.simulationData = population[member].evaluationPerformance.stepPerformanceData->data;
			//NOTE: Eventually I am going to need to specify what kind of simulation was performed...

			/* Evaluate fitness */
			runtimeStep.evaluateFitnessInstances[fitnessType]->evaluateFitness(input, output);

				
			/* Copy the output results */
			population[member].fitnessScores = output.fitness;

			#if (DEBUGGING_ENABLED == 1)
			std::cout
				<< "Fitness: "	<< output.fitness.fitness_total
				<< "\t\tKp: "	<< population[member].dChrom.Kp
				<< "\tKi: "		<< population[member].dChrom.Ki
				<< "\tKd: "		<< population[member].dChrom.Kd
				<< std::endl;
			#endif

			#if (FILE_LOGGING_ENABLED == 1) && (GA_FILELOG_MEMBER_FIT_DATA == 1)
			std::string filename = "Member" + std::to_string(member) + "_StepPerformanceFitnessData.csv";
			writeCSV_StepData(output.fitness.stepPerformanceData, filename);
			#endif
		}
	}

	else if (fitnessType == GA_FITNESS_MEAN_SQUARE_ERROR)
	{

	}
}

void FCSOptimizer::checkConvergence()
{
	/*-----------------------------------------------
	* Save best result from each round
	*-----------------------------------------------*/
	double iteration_bestFitScore = 0.0;		/* Keep track of highest fitness score found */
	PID_FitnessScores iteration_bestFit;		/* Copy of associated best fitness results */

	for (unsigned int member = 0; member < settings.advConvergenceParam.populationSize; member++)
	{
		if (parents[member].fitnessScores.fitness_total > iteration_bestFitScore)
		{
			iteration_bestFitScore = parents[member].fitnessScores.fitness_total;
			iteration_bestFit = parents[member].fitnessScores;
		}
	}

	fcs_generationalBestFitnessData.push_back(iteration_bestFit);

	#if (CONSOLE_LOGGING_ENABLED == 1)
	std::cout << "Iteration Best Fit: " << iteration_bestFitScore << "\n" << std::endl;
	#endif

	/*-----------------------------------------------
	* Break based on generation limits
	*-----------------------------------------------*/
	if (currentGeneration >= settings.advConvergenceParam.generationLimit)
	{
		#ifdef _DEBUG
		if (settings.advConvergenceParam.generationLimit == 0)
			std::cout << "You forgot to specify a generational limit for the algorithm." << std::endl;
		#endif

		#if (CONSOLE_LOGGING_ENABLED == 1)
		std::cout << "Generational Limit Reached! Done." << std::endl;
		#endif
		currentStatus = GA_COMPLETE;
	}

	/*-----------------------------------------------
	* Break based on finding a good solution 
	*-----------------------------------------------*/
	if (iteration_bestFitScore > 0.95) //TODO: Make this a user parameter
	{
		#if (CONSOLE_LOGGING_ENABLED == 1)
		std::cout << "\nFound a good enough solution. Done." << std::endl;
		#endif

		currentStatus = GA_COMPLETE;
	}		
}

void FCSOptimizer::filterPopulation()
{
	//const GA_METHOD_PopulationFilter filterType = currentSolverParam.filterType;

	//GA_PopulationFilterDataInput input;
	//GA_PopulationFilterDataOutput output;

	///*-----------------------------
	//* Ensure a filtering evaluation instance exists
	//*----------------------------*/
	//if (!runtimeStep.populationFilterInstances[filterType])
	//{
	//	switch (filterType)
	//	{
	//	case GA_POPULATION_STATIC_FILTER:
	//		runtimeStep.populationFilterInstances[filterType] = boost::make_shared<StaticFilter>(
	//			settings.pidControlSettings.tuningLimits.Kp.upper,
	//			settings.pidControlSettings.tuningLimits.Kp.lower,
	//			settings.pidControlSettings.tuningLimits.Ki.upper,
	//			settings.pidControlSettings.tuningLimits.Ki.lower,
	//			settings.pidControlSettings.tuningLimits.Kd.upper,
	//			settings.pidControlSettings.tuningLimits.Kd.lower);
	//		break;

	//	case GA_POPULATION_DYNAMIC_FILTER:
	//		runtimeStep.populationFilterInstances[filterType] = boost::make_shared<DynamicFilter>();
	//		break;

	//		//Add more as needed here
	//	default:
	//		std::cout << "Filtering method not configured or is unknown. You are about to crash." << std::endl;
	//		break;
	//	}
	//}

	///*-----------------------------
	//* Perform the filtering operation
	//*----------------------------*/
	//
	///* Populate the input struct */
	//std::transform(population.begin(), population.end(), std::back_inserter(input.currentGlobalFitScores),
	//	[](const FCSOptimizer_PopulationMember& input) { return input.fitnessScores.fitness_total; });

	//input.static_performanceThreshold = 0.4;

	///* Filter */
	//runtimeStep.populationFilterInstances[filterType]->filter(input, output);


	///* Replace the indicated population members */
	//for (int member = 0; member < output.replacedMemberIndexes.size(); member++)
	//	population[output.replacedMemberIndexes[member]].dChrom = output.replacementPIDValues[member];
}

void FCSOptimizer::selectParents(PopulationType& population)
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
			runtimeStep.selectParentInstances[selectType] = boost::make_shared<TournamentSelection>(settings.advConvergenceParam.populationSize);
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
		std::transform(population.begin(), population.end(), std::back_inserter(input.popGlobalFitScores),
			[](const FCSOptimizer_PopulationMember& input) { return input.fitnessScores.fitness_total; });
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

void FCSOptimizer::breedGeneration(PopulationType& population)
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
		input.crossoverPoint = 13;
		input.swap_both_chrom_halves = true;
		input.swap_lower_chrom_half = false;

		/* Fancy iterator to copy over the relevant genetic material */
		std::transform(population.begin(), population.end(), std::back_inserter(input.d_chrom),
			[](const FCSOptimizer_PopulationMember& member){ return member.dChrom; });
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

void FCSOptimizer::mutateGeneration(PopulationType& population)
{
	const GA_METHOD_MutateType mutateType = currentSolverParam.mutateType;
	const size_t popSize = settings.advConvergenceParam.populationSize;

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


	//Put the mutated data back into the population (should be children)
	if (currentSolverParam.chromType == MAPPING_TYPE_REAL && output.chromType == MAPPING_TYPE_REAL)
	{
		for (int i = 0; i < popSize; i++)
			population[i].dChrom = output.d_chrom[i];
	}
	else if (currentSolverParam.chromType == MAPPING_TYPE_BIT_FIELD && output.chromType == MAPPING_TYPE_BIT_FIELD)
	{
		for (int i = 0; i < popSize; i++)
			population[i].u16Chrom = output.u16_chrom[i];
	}
}

void FCSOptimizer::boundaryCheck(PopulationType& population)
{
	/* Force decimal point resolution so we aren't searching through an ENORMOUS space */
	enforceResolution(population);

	/* Ensure we aren't accidentally exceeding the tuning limits */
	enforceTunerLimits(population);
}

void FCSOptimizer::enforceResolution(PopulationType& population)
{
	double fracPart = 0.0;
	double intPart = 0.0;

	for (int member = 0; member < settings.advConvergenceParam.populationSize; member++)
	{
		/* Decompose the data into integral and fractional parts */
		fracPart = std::modf(population[member].dChrom.Kp, &intPart);

		/* Shift up, truncate, shift down data */
		fracPart *= std::pow(10.0, (int)currentSolverParam.resolutionType);
		fracPart = floor(fracPart);
		fracPart /= std::pow(10.0, (int)currentSolverParam.resolutionType);

		/* Reassign truncated value */
		population[member].dChrom.Kp = intPart + fracPart;

		/*KI*/
		fracPart = std::modf(population[member].dChrom.Ki, &intPart);

		fracPart *= std::pow(10.0, (int)currentSolverParam.resolutionType);
		fracPart = floor(fracPart);
		fracPart /= std::pow(10.0, (int)currentSolverParam.resolutionType);

		population[member].dChrom.Ki = intPart + fracPart;

		/*KD*/
		fracPart = std::modf(population[member].dChrom.Kd, &intPart);

		fracPart *= std::pow(10.0, (int)currentSolverParam.resolutionType);
		fracPart = floor(fracPart);
		fracPart /= std::pow(10.0, (int)currentSolverParam.resolutionType);

		population[member].dChrom.Kd = intPart + fracPart;
	}
}

void FCSOptimizer::enforceTunerLimits(PopulationType& population)
{
	if (settings.advConvergenceParam.limitingBehavior == FCS_LIMITER_FORCE_TO_BOUNDARY)
	{
		for (int i = 0; i < settings.advConvergenceParam.populationSize; i++)
		{
			//Kp
			if (population[i].dChrom.Kp > settings.pidControlSettings.tuningLimits.Kp.upper)
				population[i].dChrom.Kp = settings.pidControlSettings.tuningLimits.Kp.upper;

			else if (population[i].dChrom.Kp < settings.pidControlSettings.tuningLimits.Kp.lower)
				population[i].dChrom.Kp = settings.pidControlSettings.tuningLimits.Kp.lower;

			//Ki
			if (population[i].dChrom.Ki > settings.pidControlSettings.tuningLimits.Ki.upper)
				population[i].dChrom.Ki = settings.pidControlSettings.tuningLimits.Ki.upper;

			else if (population[i].dChrom.Ki < settings.pidControlSettings.tuningLimits.Ki.lower)
				population[i].dChrom.Ki = settings.pidControlSettings.tuningLimits.Ki.lower;

			//Kd
			if (population[i].dChrom.Kd > settings.pidControlSettings.tuningLimits.Kd.upper)
				population[i].dChrom.Kd = settings.pidControlSettings.tuningLimits.Kd.upper;

			else if (population[i].dChrom.Kd < settings.pidControlSettings.tuningLimits.Kd.lower)
				population[i].dChrom.Kd = settings.pidControlSettings.tuningLimits.Kd.lower;
		}
	}

	else if (settings.advConvergenceParam.limitingBehavior == FCS_LIMITER_REGENERATE_CHROMOSOME)
	{
		PID_RNG.Kp->acquireEngine();
		PID_RNG.Ki->acquireEngine();
		PID_RNG.Kd->acquireEngine();

		for (int i = 0; i < settings.advConvergenceParam.populationSize; i++)
		{
			//Kp
			if (population[i].dChrom.Kp > settings.pidControlSettings.tuningLimits.Kp.upper ||
				population[i].dChrom.Kp < settings.pidControlSettings.tuningLimits.Kp.lower)
				population[i].dChrom.Kp = PID_RNG.Kp->getDouble();


			//Ki
			if (population[i].dChrom.Ki > settings.pidControlSettings.tuningLimits.Ki.upper ||
				population[i].dChrom.Ki < settings.pidControlSettings.tuningLimits.Ki.lower)
				population[i].dChrom.Ki = PID_RNG.Ki->getDouble();

			//Kd
			if (population[i].dChrom.Kd > settings.pidControlSettings.tuningLimits.Kd.upper ||
				population[i].dChrom.Kd < settings.pidControlSettings.tuningLimits.Kd.lower)
				population[i].dChrom.Kd = PID_RNG.Kd->getDouble();
		}

		PID_RNG.Kp->releaseEngine();
		PID_RNG.Ki->releaseEngine();
		PID_RNG.Kd->releaseEngine();
	}
}

PopulationType FCSOptimizer::sortPopulation(PopulationType* parents, PopulationType* children)
{
	const GA_METHOD_Sorting sortType = currentSolverParam.sortingType;
	const int popSize = settings.advConvergenceParam.populationSize;

	GA_SortingInput input;
	GA_SortingOutput output;

	/* Ensure an instance of the mutation type exists before calling the mutate function */
	if (!runtimeStep.mutateInstances[sortType])
	{
		switch (sortType)
		{
		case GA_SORT_FAST_NONDOMINATED:
			runtimeStep.sortingInstances[sortType] = boost::make_shared<FastNonDominatedSort>();
			break;

			//Add more as needed here
		default:
			std::cout << "Sorting method not configured or is unknown. You are about to crash." << std::endl;
			break;
		}
	}

	if (parents == NULL)
		throw std::invalid_argument("You must input a parent population. Child population is optional.");
	else
	{
		for (int member = 0; member < parents->size(); member++)
			input.parentChildFitScores.push_back((*parents)[member].fitnessScores);
	}
	
	/* This is an optional argument */
	if (children != NULL)
	{
		for (int member = 0; member < children->size(); member++)
			input.parentChildFitScores.push_back((*children)[member].fitnessScores);
	}

	/* Execute the sorting algorithm */
	runtimeStep.sortingInstances[sortType]->sort(input, output);


	/* Assign new parents from output. Only the top 'popSize' results will be selected. */
	PopulationType newParents;
	if (children == NULL)
	{
		//New parents ONLY come from the function's 'parent' input. Don't bother checking 
		//if there is a size mismatch as the output is guaranteed to be of size 'popSize'
		for (int i = 0; i < popSize; i++)
		{
			int memberIdx = output.sortedPopulation[i];
			
			newParents.push_back((*parents)[memberIdx]);
		}
	}
	else
	{
		//New parents come from BOTH 'parent' and 'children' inputs 
		for (int i = 0; i < popSize; i++)
		{
			int memberIdx = output.sortedPopulation[i];

			if (memberIdx > parents->size() - 1)
				newParents.push_back((*children)[memberIdx - parents->size()]);
			else
				newParents.push_back((*parents)[memberIdx]);
		}
	}

	return newParents;
}


bool FCSOptimizer::algorithmNeedsNewSteps()
{
	//Calculate average slope of global fitness and sub scores??
	return false;
}

void FCSOptimizer::updateAlgorithmSteps()
{
	
}


void FCSOptimizer::gatherStatisticalData()
{
	double KP_Val, KI_Val, KD_Val;

	double avgGlobFit = 0.0;
	double avgPOSFit = 0.0;
	double avgSSERFit = 0.0; 
	double avgTSFit = 0.0;
	double avgTRFit = 0.0;

	for (int member = 0; member < settings.advConvergenceParam.populationSize; member++)
	{
		KP_Val = parents[member].dChrom.Kp;
		KI_Val = parents[member].dChrom.Ki;
		KD_Val = parents[member].dChrom.Kd;

		/*---------------------------------------------
		* Update the running tally of global solution mean/variance
		*----------------------------------------------*/
		GlobalChromStats.Kp(KP_Val);
		GlobalChromStats.Ki(KI_Val);
		GlobalChromStats.Kd(KD_Val);

		/*---------------------------------------------
		* Update the generational tally of mean/variance
		*----------------------------------------------*/
		GenerationalChromStats.Kp[currentGeneration](KP_Val);
		GenerationalChromStats.Ki[currentGeneration](KI_Val);
		GenerationalChromStats.Kd[currentGeneration](KD_Val);

		/*---------------------------------------------
		* Increment how many times a given PID value has been used, rounded down to the nearest integer
		*----------------------------------------------*/
		ChromOccurance.Kp[(int)std::floor(KP_Val)]++;
		ChromOccurance.Ki[(int)std::floor(KI_Val)]++;
		ChromOccurance.Kd[(int)std::floor(KD_Val)]++;


		avgGlobFit	+= parents[member].fitnessScores.fitness_total;
		avgPOSFit	+= parents[member].fitnessScores.fitness_POS;
		avgSSERFit	+= parents[member].fitnessScores.fitness_SSER;
		avgTSFit	+= parents[member].fitnessScores.fitness_TS;
		avgTRFit	+= parents[member].fitnessScores.fitness_TR;
	}

	avgGlobFit /= settings.advConvergenceParam.populationSize;
	avgPOSFit /= settings.advConvergenceParam.populationSize;
	avgSSERFit /= settings.advConvergenceParam.populationSize;
	avgTSFit /= settings.advConvergenceParam.populationSize;
	avgTRFit /= settings.advConvergenceParam.populationSize;

	avgFitness += avgGlobFit;
	avgPOSFit += avgGlobFit;
	avgSSERFit += avgGlobFit;
	avgTSFit += avgGlobFit;
	avgTRFit += avgGlobFit;

	avgFitness /= 2.0;
	avgPOSFit /= 2.0;
	avgSSERFit /= 2.0;
	avgTSFit /= 2.0;
	avgTRFit /= 2.0;
	

	averageFitness.push_back(avgFitness);
	averagePOS.push_back(avgPOSFit);
	averageSSER.push_back(avgSSERFit);
	averageTS.push_back(avgTSFit);
	averageTR.push_back(avgTRFit);

	
	//std::cout
	//	<< "Kp Variance: " << boost::accumulators::variance(GenerationalChromStats.Kp[currentGeneration]) << ", " << IdealChromVariance.Kp << "\n"
	//	<< "Ki Variance: " << boost::accumulators::variance(GenerationalChromStats.Ki[currentGeneration]) << ", " << IdealChromVariance.Ki << "\n"
	//	<< "Kd Variance: " << boost::accumulators::variance(GenerationalChromStats.Kd[currentGeneration]) << ", " << IdealChromVariance.Kd
	//	<< std::endl;

	std::cout
		<< "Kp Variance: " << boost::accumulators::variance(GlobalChromStats.Kp) << ", " << IdealChromVariance.Kp << "\n"
		<< "Ki Variance: " << boost::accumulators::variance(GlobalChromStats.Ki) << ", " << IdealChromVariance.Ki << "\n"
		<< "Kd Variance: " << boost::accumulators::variance(GlobalChromStats.Kd) << ", " << IdealChromVariance.Kd
		<< std::endl;
}

void FCSOptimizer::logData()
{
	std::ofstream file;

	std::string filename = LOG_PATH_CONVERGENCE_DATA + "avgFit.csv";
	file.open(filename);

	for (int i = 0; i < averageFitness.size()-1; i++)
		file << std::to_string(averageFitness[i]) << ",";

	file << std::to_string(averageFitness[averageFitness.size() - 1]) << "\n";
	file.close();


	filename = LOG_PATH_CONVERGENCE_DATA + "avgPOSFit.csv";
	file.open(filename);

	for (int i = 0; i < averagePOS.size()-1; i++)
		file << std::to_string(averagePOS[i]) << ",";

	file << std::to_string(averagePOS[averagePOS.size() - 1]) << "\n";
	file.close();


	filename = LOG_PATH_CONVERGENCE_DATA + "avgSSERFit.csv";
	file.open(filename);

	for (int i = 0; i < averageSSER.size()-1; i++)
		file << std::to_string(averageSSER[i]) << ",";

	file << std::to_string(averageSSER[averageSSER.size() - 1]) << "\n";
	file.close();


	filename = LOG_PATH_CONVERGENCE_DATA + "avgTSFit.csv";
	file.open(filename);

	for (int i = 0; i < averageTS.size()-1; i++)
		file << std::to_string(averageTS[i]) << ",";

	file << std::to_string(averageTS[averageTS.size() - 1]) << "\n";
	file.close();


	filename = LOG_PATH_CONVERGENCE_DATA + "avgTRFit.csv";
	file.open(filename);

	for (int i = 0; i < averageTR.size()-1; i++)
		file << std::to_string(averageTR[i]) << ",";

	file << std::to_string(averageTR[averageTR.size() - 1]) << "\n";
	file.close();
}