#include "fcs_optimizer.h"

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
	size_t num_rows = data.performance_simulation_data.rows();
	size_t num_cols = data.performance_simulation_data.cols();

	//Write each full row
	for (int row = 0; row < num_rows; row++)
	{
		for (int col = 0; col < num_cols; col++)
			if (col == num_cols - 1)
				csvFile << std::to_string(data.performance_simulation_data(row, col)*1.0);
			else
				csvFile << std::to_string(data.performance_simulation_data(row, col)*1.0) << ",";

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
/* CLASS: GACLASS_MOP */
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
			evaluateModel();		/* With the current population members, run a test scenario on the chosen model and record output data */

			evaluateFitness();		/* Use recorded performance data from last step to gauge how well it meets user goals */

			checkConvergence();		/* Given the new fitness scores, see if we have a "winner" that met user goals sufficiently */

			filterPopulation();		/* Apply random "filtering" to the population to simulate things like death/disaster/disease etc. */

			selectParents();		/* Of the current population, select those who will mate */

			breedGeneration();		/* Use the selected mating pairs from the last step to create new offspring */
		}


		/*-------------------------------------
		* Handle errors or external commands
		*-------------------------------------*/
		else
		{
			if (currentStatus == GA_HALT)
			{
				//Do any logging that might report why a halt was called. 
				break;
			}



		}


		/* Allow other waiting threads to run */
		boost::this_thread::yield();
	}

	/* Do some minor post processing here */
	//What should I do? 
	std::cout << "Exiting spawned thread" << std::endl;
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

// void FCSOptimizer::run(boost::mutex* resultsMutex, GAEngineStatistics_Vec* resultsStatistics, 
// 	GAEngineStatistics_Vec* avgResultsStatistics, int threadIndex, int trialNum)
// {
// 	currentStatus = GA_OK;
// 
// 	auto start_time = boost::chrono::high_resolution_clock::now();
// 	initPopulation();
// 
// 	while (currentStatus == GA_OK)
// 	{
// 		evaluateModel();
// 		evaluateFitness();
// 		filterPopulation();
// 		selectParents();
// 		breedGeneration();
// 		mutateGeneration();
// 
// 		checkConvergence();
// 
// 		currentIteration += 1u;
// 	}
// 
// 	auto end_time = boost::chrono::high_resolution_clock::now();
// 	boost::chrono::duration<double> execution_time = boost::chrono::duration_cast<boost::chrono::duration<double>>(end_time - start_time);
// 	
// 	print_to_console_mutex->lock();
// 	std::cout << "\nFinished optimizing " << ga_instance_optimizer_name << " in " << execution_time << std::endl;
// 	print_to_console_mutex->unlock();
// 
// 	int topPerformerIndex = reportResults(trialNum);
// 
// 	resultsMutex->lock();
// 	
// 	/* Log results for a one off trial run */
// 	resultsStatistics->data()[threadIndex].totalRunTime = execution_time;
// 	resultsStatistics->data()[threadIndex].topPerformer = GA_BestFitnessValues.data()[topPerformerIndex];
// 
// 	/* Log results for many runs on a single trial. To be averaged later. */
// 	avgResultsStatistics->push_back(resultsStatistics->data()[threadIndex]);
// 
// 	resultsMutex->unlock();
// }

/*-----------------------------------------------
* Private Functions
*-----------------------------------------------*/

void FCSOptimizer::initMemory()
{
// 	/* Expand the results logging containers */
// 	SS_StepPerformance.resize(hData.sim_data.population_size);
// 	SS_FitnessValues.resize(hData.sim_data.population_size);
// 	parentSelection.resize(hData.sim_data.population_size);

	population.resize(settings.advConvergenceParam.populationSize);

	/* Allocate the memory needed to hold all the operation functions */
	breedInstances.resize(GA_BREED_TOTAL_OPTIONS);

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


void FCSOptimizer::checkConvergence()
{
	
}

void FCSOptimizer::evaluateModel()
{
	

}

void FCSOptimizer::evaluateFitness()
{
	/* Switch how the fitness is evaluated based on current settings */
	if (settings.solverParam.fitnessType == GA_FITNESS_WEIGHTED_SUM)
	{

	}

	else if (settings.solverParam.fitnessType == GA_FITNESS_NON_DOMINATED_SORT)
	{

	}
	else
	{
		std::cout << "No valid solverParam for FitnessType!!" << std::endl;
	}
}

void FCSOptimizer::filterPopulation()
{
	if (settings.solverParam.filterType == GA_POPULATION_STATIC_FILTER)
	{

	}
	else if (settings.solverParam.filterType == GA_POPULATION_DYNAMIC_FILTER)
	{

	}
	else
	{
		std::cout << "No valid solverParam for FilterType!!" << std::endl;
	}
}

void FCSOptimizer::selectParents()
{

}

void FCSOptimizer::breedGeneration()
{
	//TODO: FUTURE WORK BELOW
	/* What if you logged a base class instance instead of filtering through all this
	junk? This works ok for smaller stuff, but as the number of methods grow, it would be
	extremely useful to have a common interface for everything. You could handle different
	kinds of settings by passing in a common struct that contains every field needed...each
	specific sub class could then use only the fields it really needs. That sounds like a pretty
	good idea. */

	GA_BreedingDataInput input;
	GA_BreedingDataOutput output;

	/* Execute the proper function, assuming it has already been bound properly */
	
	breedFunction(input, output);

	/* If not, bind the correct implementation function and then execute it */
// 	else
// 	{
// 
// 	}
}

void FCSOptimizer::mutateGeneration()
{
}


// void FCSOptimizer::evaluateModel()
// {
// 	#ifdef GA_TRACE_EVALUATE_MODEL
// 	auto start = boost::chrono::high_resolution_clock::now();
// 	#endif
// 
// 	const size_t popSize = hData.sim_data.population_size;
// 	if (ga_instance_model_name == SSModelName)
// 	{
// 		double new_pid_vals[3] = { 0.0, 0.0, 0.0 };
// 
// 		SS_NLTIV_Dynamics ss_system_dynamics(
// 			ss_user_system_model->getNumInputs(),
// 			ss_user_system_model->getNumOutputs(),
// 			ss_user_system_model->getNumStates());
// 
// 		/* Assign the simulation integration parameters */
// 		ss_system_dynamics.integrator_dt = modelSS->sim_dt;
// 		ss_system_dynamics.integrator_time_start = modelSS->sim_start_time;
// 		ss_system_dynamics.integrator_time_end = modelSS->sim_end_time;
// 
// 		/* Assign the state matrix equation model. The input U is 
// 		   given by the simulator. */
// 		ss_system_dynamics.A = ss_user_system_model->getA();
// 		ss_system_dynamics.B = ss_user_system_model->getB();
// 		ss_system_dynamics.C = ss_user_system_model->getC();
// 		ss_system_dynamics.D = ss_user_system_model->getD();
// 		ss_system_dynamics.X0 = ss_user_system_model->getX0();
// 
// 	#if defined(GA_CPU_SINGLE_THREADED) && !defined(GA_CPU_MULTI_THREADED)
// 		for (int member = 0; member < popSize; member++)
// 		{
// 			/* Update PID control values  */
// 			#ifndef FIXED_PID_VALUES
// 			new_pid_vals[0] = hData.pid_data.Kp[member];
// 			new_pid_vals[1] = hData.pid_data.Ki[member];
// 			new_pid_vals[2] = hData.pid_data.Kd[member];
// 			#else
// 			new_pid_vals[0] = 12.05;
// 			new_pid_vals[1] = 68.38;
// 			new_pid_vals[2] = 0.0;
// 			#endif
// 
// 			/* Compute Step Response and Assign Results */
// 			SS_StepPerformance[member] = modelSS->stepResponseSingleThreaded(member, ss_system_dynamics, new_pid_vals[0], new_pid_vals[1], new_pid_vals[2]);
// 		}
// 	#endif
// 
// 	#if defined(GA_CPU_MULTI_THREADED) && !defined(GA_CPU_SINGLE_THREADED)
// 		boost::thread_group tgroup;
// 
// 		for (int member = 0; member < popSize; member++)
// 		{
// 			/* Update PID control values  */
// 			#ifndef FIXED_PID_VALUES
// 			new_pid_vals[0] = hData.pid_data.Kp[member];
// 			new_pid_vals[1] = hData.pid_data.Ki[member];
// 			new_pid_vals[2] = hData.pid_data.Kd[member];
// 			#else
// 			new_pid_vals[0] = 12.05;
// 			new_pid_vals[1] = 68.38;
// 			new_pid_vals[2] = 0.0;
// 			#endif
// 
// 			/* Spawn a new thread to compute system step response for each member */
// 			tgroup.create_thread(boost::bind(&StateSpaceModel::stepResponseMultiThreaded, modelSS.get(),
// 				member,
// 				ss_system_dynamics,
// 				boost::ref(SS_StepPerformance),
// 				boost::ref(SS_StepPerformance_mutex),
// 				new_pid_vals[0],
// 				new_pid_vals[1],
// 				new_pid_vals[2]
// 			));
// 		}
// 
// 		tgroup.join_all();
// 	#endif
// 	}
// 
// 	#ifdef GA_TRACE_EVALUATE_MODEL
// 	auto stop = boost::chrono::high_resolution_clock::now();
// 	for (int member = 0; member < popSize; member++)
// 	{
// 		std::string filename = "trace/evaluate_model/" + ga_instance_optimizer_name + "_member_data_" + std::to_string(member) + ".csv";
// 		writeCSV_StepData(SS_StepPerformance[member], filename);
// 	}
// 	#endif
// }
// 
// void FCSOptimizer::evaluateFitness()
// {
// 	#ifdef GA_TRACE_EVALUATE_FITNESS
// 	auto start = boost::chrono::high_resolution_clock::now();
// 	#endif
// 
// 	if (ga_instance_model_name == SSModelName)
// 	{
// 	#if defined(GA_CPU_SINGLE_THREADED) && !defined(GA_CPU_MULTI_THREADED)
// 		if (ga_instance_step_methods.fitnessType == GA_FITNESS_WEIGHTED_SUM)
// 		{
// 			WeightedSum ws(SINGLE_THREADED);
// 			ws.calculateFitness(SS_StepPerformance, ga_instance_pid_config_data, &SS_FitnessValues);
// 		}
// 	#endif
// 
// 	#if defined(GA_CPU_MULTI_THREADED) && !defined(GA_CPU_SINGLE_THREADED)
// 		if (ga_instance_step_methods.fitnessType == GA_FITNESS_WEIGHTED_SUM)
// 		{
// 			WeightedSum ws(MULTI_THREADED);
// 			ws.calculateFitness(SS_StepPerformance, ga_instance_pid_config_data, &SS_FitnessValues);
// 		}
// 	#endif
// 	}
// 
// 	#ifdef GA_TRACE_EVALUATE_FITNESS
// 	auto stop = boost::chrono::high_resolution_clock::now();
// 	#endif
// }
// 
// void FCSOptimizer::filterPopulation()
// {
// 	#ifdef GA_TRACE_FILTER_POPULATION
// 	auto start = boost::chrono::high_resolution_clock::now();
// 	#endif
// 
// 	if (ga_instance_model_name == SSModelName)
// 	{
// 		#if defined(GA_CPU_SINGLE_THREADED) || defined(GA_CPU_MULTI_THREADED)
// 		/*-----------------------------------
// 		* Elite Performance Set
// 		*-----------------------------------*/
// 		if (GA_ElitistSolutions.EliteFitSolutions.size() < 10)
// 		{
// 			/* First time around, just update the variables directly */
// 			if (GA_BestFitnessValues.empty())
// 			{
// 				GA_ElitistSolutions.worst_performer_index = 0;
// 				GA_ElitistSolutions.worst_performer_value = 1.0;
// 			}
// 			/* Otherwise, push solution to the end */
// 			else
// 				GA_ElitistSolutions.EliteFitSolutions.push_back(GA_BestFitnessValues.back());
// 		}
// 
// 		/* Assuming the buffer is full, decide whether or not to add a solution in. */
// 		else
// 		{
// 			/* IF the latest "BestFitnessValue" is higher than the current worst fitness, this means 
// 			the Elitist solution set must be updated to include the new value. The goal is for the elite 
// 			set to update to higher and higher performing members each time, assuming such a member exists. */
// 			if (GA_ElitistSolutions.worst_performer_value < GA_BestFitnessValues.back().global_fitness)
// 			{
// 				GA_ElitistSolutions.EliteFitSolutions.data()[GA_ElitistSolutions.worst_performer_index] =
// 					GA_BestFitnessValues.back();
// 			}
// 		}
// 
// 		/* Find the new lowest performing individual out of the full set */
// 		double worst_performer = 2.0;
// 		for (int i = 0; i < GA_ElitistSolutions.EliteFitSolutions.size(); i++)
// 		{
// 			if (GA_ElitistSolutions.EliteFitSolutions.data()[i].global_fitness < worst_performer)
// 			{
// 				GA_ElitistSolutions.worst_performer_value = GA_ElitistSolutions.EliteFitSolutions.data()[i].global_fitness;
// 				GA_ElitistSolutions.worst_performer_index = i;
// 
// 				worst_performer = GA_ElitistSolutions.worst_performer_value;
// 			}
// 		}
// 
// 		/*-----------------------------------
// 		* Rejection of Population Members:
// 		* In this case, it is a simulation of natural disaster/survival of fittest
// 		*-----------------------------------*/
// 		boost::container::vector<int> rejectionIdxs;
// 
// 		/* Create a random number generator that is used to simulate how 
// 		   many population members are killed off at each generation. */
// 		int maxReplacements = (int)floor(hData.sim_data.population_size*0.4);
// 		std::random_device rd1;
// 		std::mt19937 rng1(rd1());
// 		std::uniform_int_distribution<uint32_t> uniform_int(0, maxReplacements);
// 
// 		/* Use a static threshold to filter through member performance */
// 		if (ga_instance_step_methods.filterType == GA_POPULATION_STATIC_FILTER)
// 		{
// 			int totalReplacements = uniform_int(rng1);
// 			double filter_threshold = 0.1;
// 
// 			for (int i = 0; i < totalReplacements; i++)
// 			{
// 				if (SS_FitnessValues.data()[i].global_fitness < filter_threshold)
// 					rejectionIdxs.push_back(i);
// 			}
// 		}
// 
// 		/* Dynamically choose a threshold to filter through member performance */
// 		if (ga_instance_step_methods.filterType == GA_POPULATION_DYNAMIC_FILTER)
// 		{
// 			std::random_device rd2;
// 			std::mt19937 rng2(rd2());
// 			std::uniform_real_distribution<double> uniform_dbl(0.1, 0.7);
// 
// 			int totalReplacements = uniform_int(rng1);
// 			double filter_threshold = uniform_dbl(rng2);
// 
// 			for (int i = 0; i < totalReplacements; i++)
// 			{
// 				if (SS_FitnessValues.data()[i].global_fitness < filter_threshold)
// 					rejectionIdxs.push_back(i);
// 			}
// 		}
// 
// 		/*-----------------------------------
// 		* Creation of New Members
// 		*-----------------------------------*/
// 		//Random, elite, copy of existing solutions?
// 
// 		double kpl = ga_instance_pid_config_data->tuningLimits.Kp_limits_lower;
// 		double kpu = ga_instance_pid_config_data->tuningLimits.Kp_limits_upper;
// 
// 		double kil = ga_instance_pid_config_data->tuningLimits.Ki_limits_lower;
// 		double kiu = ga_instance_pid_config_data->tuningLimits.Ki_limits_upper;
// 
// 		double kdl = ga_instance_pid_config_data->tuningLimits.Kd_limits_lower;
// 		double kdu = ga_instance_pid_config_data->tuningLimits.Kd_limits_upper;
// 
// 		/* Pure random generation of replacement members */
// 		for (int i = 0; i < rejectionIdxs.size(); i++)
// 		{
// 			hData.pid_data.Kp[rejectionIdxs[i]] = uniformRandomNumber(kpl, kpu);
// 			hData.pid_data.Ki[rejectionIdxs[i]] = uniformRandomNumber(kil, kiu);
// 			hData.pid_data.Kd[rejectionIdxs[i]] = uniformRandomNumber(kdl, kdu);
// 		}
// 		#endif
// 	}
// 
// 	#ifdef GA_TRACE_FILTER_POPULATION
// 	auto stop = boost::chrono::high_resolution_clock::now();
// 	#endif
// }
// 
// void FCSOptimizer::selectParents()
// {
// 	#ifdef GA_TRACE_SELECT_PARENTS
// 	auto start = boost::chrono::high_resolution_clock::now();
// 	#endif
// 
// 	const size_t popSize = hData.sim_data.population_size;
// 
// 	if (ga_instance_model_name == SSModelName)
// 	{
// 		/* For now, only use single threaded parent selection due to small
// 		   population sizes. Switching over to a multi-threaded version won't
// 		   give much of a speed up.*/
// 	#if defined(GA_CPU_SINGLE_THREADED) || defined(GA_CPU_MULTI_THREADED)
// 		if (ga_instance_step_methods.selectType == GA_SELECT_RANDOM)
// 		{
// 			RandomSelection RS(SINGLE_THREADED);
// 			RS.selectParents(&parentSelection);
// 		}
// 
// 		if (ga_instance_step_methods.selectType == GA_SELECT_TOURNAMENT)
// 		{
// 			TournamentSelection TS(SINGLE_THREADED);
// 			TS.selectParents(SS_FitnessValues, &parentSelection);
// 		}
// 	#endif
// 	}
// 
// 	#ifdef GA_TRACE_SELECT_PARENTS
// 	auto stop = boost::chrono::high_resolution_clock::now();
// 	#endif
// }
// 
// void FCSOptimizer::breedGeneration()
// {
// 	#ifdef GA_TRACE_BREED_GENERATION
// 	auto trace_start_time = boost::chrono::high_resolution_clock::now();
// 	#endif
// 
// 	const size_t popSize = hData.sim_data.population_size;
// 	if (ga_instance_model_name == SSModelName)
// 	{
// 	#ifdef GA_ENFORCE_RESOLUTION_BG
// 	enforceResolution();
// 	#endif
// 
// 	#if defined(GA_CPU_SINGLE_THREADED) && !defined(GA_CPU_MULTI_THREADED)
// 	
// 		/*-----------------------------------------------
// 		* Simple Crossover 
// 		*-----------------------------------------------*/
// 		if (ga_instance_step_methods.breedType == GA_BREED_SIMPLE_CROSSOVER)
// 		{
// 			SimpleCrossover SC(SINGLE_THREADED);
// 			SC.breed(parentSelection, &hData.pid_data, &bredChromosomes, &mapCoefficients_Kp, &mapCoefficients_Ki, &mapCoefficients_Kd);
// 		}
// 
// 		/*-----------------------------------------------
// 		* Dynamic Crossover
// 		*-----------------------------------------------*/
// 		if (ga_instance_step_methods.breedType == GA_BREED_DYNAMIC_CROSSOVER)
// 		{
// 			DynamicCrossover DC(SINGLE_THREADED);
// 			DC.breed(parentSelection, &hData.pid_data, &bredChromosomes, &mapCoefficients_Kp, &mapCoefficients_Ki, &mapCoefficients_Kd);
// 		}
// 
// 		/*-----------------------------------------------
// 		* Fixed Ratio Crossover
// 		*-----------------------------------------------*/
// 		if (ga_instance_step_methods.breedType == GA_BREED_FIXED_RATIO_CROSSOVER)
// 		{
// 			FixedRatioCrossover FRC(SINGLE_THREADED, 0.25);
// 			FRC.breed(parentSelection, &hData.pid_data, &bredChromosomes, &mapCoefficients_Kp, &mapCoefficients_Ki, &mapCoefficients_Kd);
// 		}
// 
// 	
// 	#endif
// 
// 	#if defined(GA_CPU_MULTI_THREADED) && !defined(GA_CPU_SINGLE_THREADED)
// 		/*-----------------------------------------------
// 		* Simple Crossover
// 		*-----------------------------------------------*/
// 		if (ga_instance_step_methods.breedType == GA_BREED_SIMPLE_CROSSOVER)
// 		{
// 			#ifdef DEBUGGING_ENABLED
// 			size_t pidSize = hData.pid_data.Kd.size();
// 			#endif
// 
// 			SimpleCrossover SC(MULTI_THREADED);
// 			SC.breed(parentSelection, &hData.pid_data, &bredChromosomes, &mapCoefficients_Kp, &mapCoefficients_Ki, &mapCoefficients_Kd);
// 		}
// 
// 		/*-----------------------------------------------
// 		* Dynamic Crossover
// 		*-----------------------------------------------*/
// 		if (ga_instance_step_methods.breedType == GA_BREED_DYNAMIC_CROSSOVER)
// 		{
// 			DynamicCrossover DC(MULTI_THREADED);
// 			DC.breed(parentSelection, &hData.pid_data, &bredChromosomes, &mapCoefficients_Kp, &mapCoefficients_Ki, &mapCoefficients_Kd);
// 		}
// 
// 		/*-----------------------------------------------
// 		* Fixed Ratio Crossover
// 		*-----------------------------------------------*/
// 		if (ga_instance_step_methods.breedType == GA_BREED_FIXED_RATIO_CROSSOVER)
// 		{
// 			FixedRatioCrossover FRC(MULTI_THREADED, 0.25);
// 			FRC.breed(parentSelection, &hData.pid_data, &bredChromosomes, &mapCoefficients_Kp, &mapCoefficients_Ki, &mapCoefficients_Kd);
// 		}
// 	#endif
// 	}
// 	
// 	#ifdef GA_TRACE_BREED_GENERATION
// 	auto trace_end_time = boost::chrono::high_resolution_clock::now();
// 	#endif
// }
// 
// void FCSOptimizer::mutateGeneration()
// {
// 	#ifdef GA_TRACE_MUTATE_GENERATION
// 	auto start = boost::chrono::high_resolution_clock::now();
// 	#endif
// 
// 	const size_t popSize = hData.sim_data.population_size;
// 	if (ga_instance_model_name == SSModelName)
// 	{
// 	#if defined(GA_CPU_SINGLE_THREADED) && !defined(GA_CPU_MULTI_THREADED)
// 		/*-----------------------------------------------
// 		* Bit Flipping Mutator
// 		*-----------------------------------------------*/
// 		if (ga_instance_step_methods.mutateType == GA_MUTATE_BIT_FLIP)
// 		{
// 			BitFlipMutator mutator(SINGLE_THREADED, ga_instance_step_methods.mutateProbabilityType);
// 
// 			mutator.mutate(&bredChromosomes, ga_instance_convergence_criteria, ga_instance_pid_config_data,
// 				&hData.pid_data, &mapCoefficients_Kp, &mapCoefficients_Ki, &mapCoefficients_Kd);
// 		}
// 
// 		/*-----------------------------------------------
// 		* Add-Subtract Mutator
// 		*-----------------------------------------------*/
// 		if (ga_instance_step_methods.mutateType == GA_MUTATE_ADD_SUB)
// 		{
// 			AddSubMutator mutator(SINGLE_THREADED,
// 				ga_instance_step_methods.mutateProbabilityType,
// 				ga_instance_step_methods.resolutionType);
// 
// 			mutator.mutate(&bredChromosomes,
// 				ga_instance_convergence_criteria,
// 				ga_instance_pid_config_data,
// 				&hData.pid_data,
// 				&mapCoefficients_Kp,
// 				&mapCoefficients_Ki,
// 				&mapCoefficients_Kd);
// 		}
// 	#endif
// 
// 	#if defined(GA_CPU_MULTI_THREADED) && !defined(GA_CPU_SINGLE_THREADED)
// 		/*-----------------------------------------------
// 		* Bit Flipping Mutator
// 		*-----------------------------------------------*/
// 		if (ga_instance_step_methods.mutateType == GA_MUTATE_BIT_FLIP)
// 		{
// 			BitFlipMutator mutator(MULTI_THREADED, ga_instance_step_methods.mutateProbabilityType);
// 
// 			mutator.mutate(&bredChromosomes, ga_instance_convergence_criteria, ga_instance_pid_config_data,
// 				&hData.pid_data, &mapCoefficients_Kp, &mapCoefficients_Ki, &mapCoefficients_Kd);
// 		}
// 
// 		/*-----------------------------------------------
// 		* Add-Subtract Mutator
// 		*-----------------------------------------------*/
// 		if (ga_instance_step_methods.mutateType == GA_MUTATE_ADD_SUB)
// 		{
// 			AddSubMutator mutator(MULTI_THREADED,
// 				ga_instance_step_methods.mutateProbabilityType,
// 				ga_instance_step_methods.resolutionType);
// 
// 			mutator.mutate(&bredChromosomes,
// 				ga_instance_convergence_criteria,
// 				ga_instance_pid_config_data,
// 				&hData.pid_data,
// 				&mapCoefficients_Kp,
// 				&mapCoefficients_Ki,
// 				&mapCoefficients_Kd);
// 		}
// 	#endif
// 
// 	#ifdef GA_ENFORCE_RESOLUTION_MG
// 		enforceResolution();
// 	#endif
// 	}
// 	
// 	#ifdef GA_TRACE_MUTATE_GENERATION
// 	auto stop = boost::chrono::high_resolution_clock::now();
// 	#endif
// }
// 
// void FCSOptimizer::checkConvergence()
// {
// 	/*-----------------------------------------------
// 	* Save best result from each round
// 	*-----------------------------------------------*/
// 	double iteration_bestFitGlobal = 0.0;
// 	PID_FitnessScores iteration_bestFit;
// 
// 	for (int member = 0; member < hData.sim_data.population_size; member++)
// 	{
// 		if (SS_FitnessValues.data()[member].global_fitness > iteration_bestFitGlobal)
// 		{
// 			iteration_bestFit = SS_FitnessValues.data()[member];
// 			iteration_bestFitGlobal = SS_FitnessValues.data()[member].global_fitness;
// 		}
// 	}
// 
// 	GA_BestFitnessValues.push_back(iteration_bestFit);
// 	
// 	#ifdef GA_REPORT_DATA_CHECK_CONVERGENCE
// 	print_to_console_mutex->lock();
// 	std::cout << "Iteration Best Fit: " << iteration_bestFitGlobal << std::endl;
// 	print_to_console_mutex->unlock();
// 	#endif
// 
// 	/*-----------------------------------------------
// 	* Break based on generation limits
// 	*-----------------------------------------------*/
// 	if (currentIteration >= ga_instance_convergence_criteria->generationLimit)
// 		currentStatus = GA_COMPLETE;
// 	
// 	/*-----------------------------------------------
// 	* Break based on finding a good solution 
// 	*-----------------------------------------------*/
// 	if (iteration_bestFitGlobal > 0.9)
// 	{
// 		print_to_console_mutex->lock();
// 		std::cout << "\nFound a good enough solution. Done." << std::endl;
// 		print_to_console_mutex->unlock();
// 		currentStatus = GA_COMPLETE;
// 	}
// 		
// }
// 
// void FCSOptimizer::enforceResolution()
// {
// 	for (int member = 0; member < hData.sim_data.population_size; member++)
// 	{
// 		double fracPart = 0.0;
// 		double intPart = 0.0;
// 
// 		/* Decompose the data into integral and fractional parts */
// 		fracPart = std::modf(hData.pid_data.Kp.data()[member], &intPart);
// 
// 		/* Shift up, truncate, shift down data */
// 		fracPart *= std::pow(10.0, (int)ga_instance_step_methods.resolutionType);
// 		fracPart = floor(fracPart);
// 		fracPart /= std::pow(10.0, (int)ga_instance_step_methods.resolutionType);
// 
// 		/* Reassign truncated value */
// 		hData.pid_data.Kp.data()[member] = intPart + fracPart;
// 	}
// }
// 
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
// 
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