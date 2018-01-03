#include "ga_engine.h"

//////////////////////////////////////////////////////////////////
/* CLASS: GAEngine_PID */
//////////////////////////////////////////////////////////////////
/*-----------------------------------------------
* Constructor/Destructor
*-----------------------------------------------*/
GAEngine_PID::GAEngine_PID(GAEngineStatistics_Vec *statisticsData)
{
	avgThreadStatistics = statisticsData;
}

GAEngine_PID::~GAEngine_PID()
{
}

/*-----------------------------------------------
* Public Functions
*-----------------------------------------------*/
void GAEngine_PID::createOptimizer(std::string optimizer, GA_AlgorithmMethods alg_step_types)
{
	//Create a new Multi-Objective Optimizer instance and assign to the vector
	ga_optimizer.push_back(std::make_shared<GA_MOP>(alg_step_types));

	//Based on the latest index value of the above, map the optimizer to the accessIndex
	size_t accessIndex = (ga_optimizer.size() - 1);
	thread_mapping.insert(std::make_pair(optimizer, accessIndex));

	//Assign the optimizer string to the Engine's vector
	size_t opNameVecSize = optimizer_names.size();

	if (accessIndex < opNameVecSize)
		optimizer_names[accessIndex] = optimizer;

	else if (accessIndex == opNameVecSize)
		optimizer_names.push_back(optimizer);

	else if (accessIndex > opNameVecSize)
	{
		optimizer_names.resize((accessIndex + 1), std::string(""));
		optimizer_names[accessIndex] = optimizer;
	}

	//Tell the actual optimizer instance what its name is.
	assignOptimizer_name(optimizer, optimizer);

	//Give the optimizer a reference to the console print mutex
	ga_optimizer.back()->registerOutputMutex(&printMutex);

	//Resize the results container to be the same size as the optimizer
	threadStatistics.resize(ga_optimizer.size());
}

void GAEngine_PID::attachModel(std::string optimizer, GAModel_sPtr model_type, SS_NLTIV_ModelBase_sPtr model_system, std::string model_name, int processor)
{
	uint32_t accessIndex = thread_mapping[optimizer];
	size_t modelVecSize = model.size();
	size_t modelNameVecSize = model_names.size();

	/*-------------------------
	* Handle model_type placement
	*-------------------------*/
	/* Write directly to the vector if enough room is allocated. */
	if (accessIndex < modelVecSize)
		model[accessIndex] = model_type;

	/* If this is the next value to be inserted, do it. */
	else if (accessIndex == modelVecSize)
		model.push_back(model_type);

	/* If the accessIndex is way ahead of the vector size, resize
	the vector so the new data can fit. This can happen occasionally
	non-sequential model assignments. */
	else if (accessIndex > modelVecSize)
	{
		model.resize((accessIndex + 1), nullptr);
		model[accessIndex] = model_type;
	}

	/*-------------------------
	* Handle model_name placement
	*-------------------------*/
	if (accessIndex < modelNameVecSize)
		model_names[accessIndex] = model_name;

	else if (accessIndex == modelNameVecSize)
		model_names.push_back(model_name);

	else if (accessIndex > modelNameVecSize)
	{
		model_names.resize((accessIndex + 1), std::string(""));
		model_names[accessIndex] = model_name;
	}

	/*-------------------------
	* Inform the actual optimizer of the correct model to use
	*-------------------------*/
	assignOptimizer_model(optimizer, model[accessIndex], model_system, processor);
}

void GAEngine_PID::attachObjective(std::string optimizer, PID_ControlGoals_sPtr pid_parameters, std::string pid_name)
{
	uint32_t accessIndex = thread_mapping[optimizer];
	size_t pid_paramVecSize = pid.size();
	size_t pid_nameVecSize = pid_names.size();

	/*-------------------------
	* Handle pid_parameters placement
	*-------------------------*/
	if (accessIndex < pid_paramVecSize)
		pid[accessIndex] = pid_parameters;

	else if (accessIndex == pid_paramVecSize)
		pid.push_back(pid_parameters);

	else if (accessIndex > pid_paramVecSize)
	{
		pid.resize((accessIndex + 1), nullptr);
		pid[accessIndex] = pid_parameters;
	}

	/*-------------------------
	* Handle pid_name placement
	*-------------------------*/
	if (accessIndex < pid_nameVecSize)
		pid_names[accessIndex] = pid_name;

	else if (accessIndex == pid_nameVecSize)
		pid_names.push_back(pid_name);

	else if (accessIndex > pid_nameVecSize)
	{
		pid_names.resize((accessIndex + 1), std::string(""));
		pid_names[accessIndex] = pid_name;
	}

	/*-------------------------
	* Inform the actual optimizer of the correct PID data to use
	*-------------------------*/
	assignOptimizer_objective(optimizer, pid[accessIndex]);
}

void GAEngine_PID::attachConvergence(std::string optimizer, GA_ConverganceCriteria_sPtr convergence_param, std::string convergence_name)
{
	uint32_t accessIndex = thread_mapping[optimizer];
	size_t conv_paramVecSize = ga_convergence.size();
	size_t conv_nameVecSize = convergence_names.size();

	/*-------------------------
	* Handle pid_parameters placement
	*-------------------------*/
	if (accessIndex < conv_paramVecSize)
		ga_convergence[accessIndex] = convergence_param;

	else if (accessIndex == conv_paramVecSize)
		ga_convergence.push_back(convergence_param);

	else if (accessIndex > conv_paramVecSize)
	{
		ga_convergence.resize((accessIndex + 1), nullptr);
		ga_convergence[accessIndex] = convergence_param;
	}

	/*-------------------------
	* Handle pid_name placement
	*-------------------------*/
	if (accessIndex < conv_nameVecSize)
		convergence_names[accessIndex] = convergence_name;

	else if (accessIndex == conv_nameVecSize)
		convergence_names.push_back(convergence_name);

	else if (accessIndex > conv_nameVecSize)
	{
		convergence_names.resize((accessIndex + 1), std::string(""));
		convergence_names[accessIndex] = convergence_name;
	}

	assignOptimizer_convergence(optimizer, ga_convergence[accessIndex]);
}

void GAEngine_PID::join(std::string optimizer)
{
	uint32_t accessIndex = thread_mapping[optimizer];
	if (thread[accessIndex] != nullptr)
		if (thread[accessIndex]->joinable())
			thread[accessIndex]->join();
}

void GAEngine_PID::joinAll()
{
	//Well this is nifty. Range based loop!
	for (auto & i : thread)
	{
		if (i == nullptr)
			continue;

		else if (i->joinable())
			i->join();
	}
}

void GAEngine_PID::init(std::string optimizer)
{
	uint32_t accessIndex = thread_mapping[optimizer];
	if (ga_optimizer[accessIndex] != nullptr)
		ga_optimizer[accessIndex]->init();
}

void GAEngine_PID::initAll()
{
	/* Grab container sizes for key parameters */
	size_t optimizer_names_size = optimizer_names.size();
	size_t convergence_size = ga_convergence.size();
	size_t model_size = model.size();
	size_t pid_size = pid.size();

	/* Run an algorithm IFF it has the required parameters set */
	for (int i = 0; i < optimizer_names.size(); i++)
	{
		/* Check for valid convergence criteria input */
		if (i < convergence_size)
		{
			if (ga_convergence[i] == nullptr)
				continue;
		}
		else
			continue;

		/* Check for a valid model input */
		if (i < model_size)
		{
			if (model[i] == nullptr)
				continue;
		}
		else
			continue;

		/* Check for a valid PID goal input*/
		if (i < pid_size)
		{
			if (pid[i] == nullptr)
				continue;
		}
		else
			continue;

		//If the loop makes it this far, the optimizer was set up correctly
		init(optimizer_names[i]);
	}
}

void GAEngine_PID::run(std::string optimizer, int trialNum)
{
	uint32_t accessIndex = thread_mapping[optimizer];
	size_t threadVecSize = thread.size();

	/*----------------------------
	* Create the function the thread will execute
	*----------------------------*/
	auto bindOp = boost::bind(
		&GA_MOP::run,					//Template of the function to run
		ga_optimizer[accessIndex],		//Specific instance to bind
		&statisticsMutex,
		&threadStatistics,
		avgThreadStatistics,
		accessIndex,
		trialNum);		


	/*----------------------------
	* Assign the thread to the management vector.
	* This automatically executes the thread
	*----------------------------*/
	if (accessIndex < threadVecSize)
	{
		/* Make sure that an accessed index will not override a current thread.*/
//		if (thread[accessIndex] == nullptr)
			thread[accessIndex] = std::make_shared<boost::thread>(bindOp);
// 		else
// 			std::cout << "Attempted to overwrite existing thread. Not allowed." << std::endl;
	}

	else if (accessIndex == threadVecSize)
	{
		thread.push_back(std::make_shared<boost::thread>(bindOp));
	}

	else if (accessIndex > threadVecSize)
	{
		thread.resize((accessIndex + 1), nullptr);

		//Don't have to check for nullptr here. Expanding the size of the
		//thread vector guarantees nullptr existence.
		thread[accessIndex] = std::make_shared<boost::thread>(bindOp);
	}
}

void GAEngine_PID::runAll(int trialNum)
{
	/* Grab container sizes for key parameters */
	size_t optimizer_names_size = optimizer_names.size();
	size_t convergence_size = ga_convergence.size();
	size_t model_size = model.size();
	size_t pid_size = pid.size();

	/* Run an algorithm IFF it has the required parameters set */
	for (int i = 0; i < optimizer_names.size(); i++)
	{
		if (i < convergence_size)
		{
			if (ga_convergence[i] == nullptr)
				continue;
		}
		else
			continue;

		if (i < model_size)
		{
			if (model[i] == nullptr)
				continue;
		}
		else
			continue;

		if (i < pid_size)
		{
			if (pid[i] == nullptr)
				continue;
		}
		else
			continue;

		//If the loop makes it this far, the optimizer was set up correctly
		run(optimizer_names[i], trialNum);
	}
}

void GAEngine_PID::reportResultsToConsole()
{
	std::cout << "=============================================================\n"
		<< "\t OVERALL PERFORMANCE RESULTS:\n"
		<< "============================================================="
		<< std::endl;

	/*-----------------------------------------------
	* Runtime 
	*-----------------------------------------------*/
	for (int i = 0; i < thread.size(); i++)
	{
		std::string optimizer = optimizer_names[i];
		uint32_t accessIndex = thread_mapping[optimizer];

		std::cout << optimizer << " Runtime: " << threadStatistics.data()[accessIndex].totalRunTime << std::endl;
	}

	/*-----------------------------------------------
	* Best Fit
	*-----------------------------------------------*/
	std::cout << "\n" << std::endl;
	for (int i = 0; i < thread.size(); i++)
	{
		std::string optimizer = optimizer_names[i];
		uint32_t accessIndex = thread_mapping[optimizer];

		std::cout << optimizer << " Fitness: " << threadStatistics.data()[accessIndex].topPerformer.global_fitness << std::endl;
	}

}

void GAEngine_PID::reportResultsToFile(int trialNum, double full_execution_time)
{
	#if defined(GA_CPU_SINGLE_THREADED) && !defined(GA_CPU_MULTI_THREADED)
	std::string filename = "SingleThreadedTrial" + std::to_string(trialNum) + "_Results.csv";
	#endif

	#if defined(GA_CPU_MULTI_THREADED) && !defined(GA_CPU_SINGLE_THREADED)
	std::string filename = "MultiThreadedTrial" + std::to_string(trialNum) + "_Results.csv";
	#endif

	std::ofstream csvFile;
	csvFile.open(filename);


	csvFile << "TotalTime:," << std::to_string(full_execution_time) << "\n";

	/*-----------------------------------------------
	* Runtime
	*-----------------------------------------------*/
	for (int i = 0; i < thread.size(); i++)
	{
		std::string optimizer = optimizer_names[i];
		uint32_t accessIndex = thread_mapping[optimizer];

		csvFile << optimizer << " Runtime:," << std::to_string(threadStatistics.data()[accessIndex].totalRunTime.count()) << "\n";
	}

	/*-----------------------------------------------
	* Best Fit
	*-----------------------------------------------*/
	for (int i = 0; i < thread.size(); i++)
	{
		std::string optimizer = optimizer_names[i];
		uint32_t accessIndex = thread_mapping[optimizer];

		csvFile << optimizer << " Fitness:," << std::to_string(threadStatistics.data()[accessIndex].topPerformer.global_fitness) << "\n";
	}
}

void GAEngine_PID::reportAvgResultsToFile(int trialNum, double avg_execution_time)
{
	#if defined(GA_CPU_SINGLE_THREADED) && !defined(GA_CPU_MULTI_THREADED)
	std::string filename = "SingleThreadedTrial" + std::to_string(trialNum) + "_AVGResults.csv";
	#endif

	#if defined(GA_CPU_MULTI_THREADED) && !defined(GA_CPU_SINGLE_THREADED)
	std::string filename = "MultiThreadedTrial" + std::to_string(trialNum) + "_AVGResults.csv";
	#endif

	std::ofstream csvFile;
	csvFile.open(filename);


	csvFile << "Total Time:," << std::to_string(avg_execution_time) << "\n";

	/*-----------------------------------------------
	* Runtime
	*-----------------------------------------------*/
	double avgRunTime = 0.0;
	for (int i = 0; i < avgThreadStatistics->size(); i++)
		avgRunTime += avgThreadStatistics->data()[i].totalRunTime.count();

	avgRunTime /= (double)avgThreadStatistics->size();

	csvFile << "AVG Runtime:," << std::to_string(avgRunTime) << "\n";


	/*-----------------------------------------------
	* Best Fit
	*-----------------------------------------------*/
	double avgFitness = 0.0;
	for (int i = 0; i < avgThreadStatistics->size(); i++)
		avgFitness += avgThreadStatistics->data()[i].topPerformer.global_fitness;

	avgFitness /= (double)avgThreadStatistics->size();

	csvFile << "AVG Fitness:," << std::to_string(avgFitness) << "\n";
}


/*-----------------------------------------------
* Private Functions
*-----------------------------------------------*/
void GAEngine_PID::assignOptimizer_name(std::string optimizer, std::string optimizer_name)
{
	uint32_t accessIndex = thread_mapping[optimizer];
	ga_optimizer[accessIndex]->registerName(optimizer_name);
}

void GAEngine_PID::assignOptimizer_model(std::string optimizer, GAModel_sPtr model, SS_NLTIV_ModelBase_sPtr model_system, int processor)
{
	uint32_t accessIndex = thread_mapping[optimizer];
	ga_optimizer[accessIndex]->registerModel(model, model_system, model_names[accessIndex], processor);
}

void GAEngine_PID::assignOptimizer_objective(std::string optimizer, PID_ControlGoals_sPtr objective)
{
	uint32_t accessIndex = thread_mapping[optimizer];
	ga_optimizer[accessIndex]->registerObjective(objective, pid_names[accessIndex]);
}

void GAEngine_PID::assignOptimizer_convergence(std::string optimizer, GA_ConverganceCriteria_sPtr convergence)
{
	uint32_t accessIndex = thread_mapping[optimizer];
	ga_optimizer[accessIndex]->registerConvergence(convergence, convergence_names[accessIndex]);
}