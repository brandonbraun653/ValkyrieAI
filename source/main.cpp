#include <stdlib.h>
#include <iostream>
#include "ga_config.h"
#include "ga_engine.h"
#include "logger.h"
#include "model.h"
#include "data.h"
#include "user_sim_models.h"

/* Global control flags */
#define POPULATION_SIZE 500
#define GENERATION_LIMIT 25
#define MUTATE_THRESHOLD 0.1
#define MUTATE_SEVERITY 12

#define RISE_TIME 0.16
#define SETTLING_TIME 0.25
#define PERCENT_OVERSHOOT 0.1
#define STEADY_STATE_ERROR 0.05

#define RISE_TIME_TOLERANCE 0.2
#define SETTLING_TIME_TOLERANCE 0.5
#define PERCENT_OVERSHOOT_TOLERANCE 0.25
#define STEADY_STATE_ERROR_TOLERANCE 0.25

#define NUM_OPTIMIZERS 3
#define MAX_TRIALS 6
#define AVERAGING_CYCLES 10
int trial_popSizes[MAX_TRIALS] = { 5, 10, 25, 50, 100, 250 };
// #define MAX_TRIALS 1
// int trial_popSizes[MAX_TRIALS] = { 50 };

void reportAvgResultsToFile(int trialNum, double avg_execution_time, GAEngineStatistics_Vec statisticsData);

int main()
{
	boost::container::vector<GA_ConverganceCriteria_sPtr> trial_convergence_goals(NUM_OPTIMIZERS);
	boost::container::vector<PID_ControlGoals_sPtr> trial_pid_goals(NUM_OPTIMIZERS);
	boost::container::vector<GA_AlgorithmMethods> trial_steps(NUM_OPTIMIZERS);

	boost::container::vector<SS_NLTIV_ModelBase_sPtr> trial_eval_models(7);
	trial_eval_models.data()[0] = std::make_shared<Example_SS_PID>();
	trial_eval_models.data()[1] = std::make_shared<Example_SS_PID>();
	trial_eval_models.data()[2] = std::make_shared<Example_SS_PID>();
	trial_eval_models.data()[3] = std::make_shared<Example_SS_PID>();
	trial_eval_models.data()[4] = std::make_shared<Example_SS_PID>();
	trial_eval_models.data()[5] = std::make_shared<Example_SS_PID>();
	trial_eval_models.data()[6] = std::make_shared<Example_SS_PID>();



	/* Initialize all the data that will stay the SAME for each trial */
	for (int opNum = 0; opNum < NUM_OPTIMIZERS; opNum++)
	{
		/*------------------------
		* GA Convergence Criteria
		*------------------------*/
		trial_convergence_goals.data()[opNum] = std::make_shared<GA_ConverganceCriteria>();
		trial_convergence_goals.data()[opNum]->generationLimit = GENERATION_LIMIT;
		trial_convergence_goals.data()[opNum]->mutation_severity = MUTATE_SEVERITY;
		trial_convergence_goals.data()[opNum]->mutation_threshold = MUTATE_THRESHOLD;

		/*------------------------
		* GA Execution Steps
		*------------------------*/
		trial_steps.data()[opNum].breedType =				GA_BREED_DYNAMIC_CROSSOVER;
		trial_steps.data()[opNum].fitnessType =				GA_FITNESS_WEIGHTED_SUM;
		trial_steps.data()[opNum].mutateProbabilityType =	GA_MUTATE_PROBABILITY_EXPONENTIAL;
		trial_steps.data()[opNum].mutateType =				GA_MUTATE_ADD_SUB;
		trial_steps.data()[opNum].filterType =				GA_POPULATION_STATIC_FILTER;
		trial_steps.data()[opNum].selectType =				GA_SELECT_RANDOM;
		trial_steps.data()[opNum].resolutionType =			GA_RESOLUTION_3DP;

		/*------------------------
		* PID Performance Goals
		*------------------------*/
		trial_pid_goals.data()[opNum] = std::make_shared<PID_ControlGoals>();
		trial_pid_goals.data()[opNum]->performance_goals.percentOvershoot_goal = PERCENT_OVERSHOOT;
		trial_pid_goals.data()[opNum]->performance_goals.riseTime_goal = RISE_TIME;
		trial_pid_goals.data()[opNum]->performance_goals.settlingTime_goal = SETTLING_TIME;
		trial_pid_goals.data()[opNum]->performance_goals.steadyStateError_goal = STEADY_STATE_ERROR;

		trial_pid_goals.data()[opNum]->performance_tolerance.percentOvershoot_pcntTol = PERCENT_OVERSHOOT_TOLERANCE;
		trial_pid_goals.data()[opNum]->performance_tolerance.riseTime_pcntTol = RISE_TIME_TOLERANCE;
		trial_pid_goals.data()[opNum]->performance_tolerance.settlingTime_pcntTol = SETTLING_TIME_TOLERANCE;
		trial_pid_goals.data()[opNum]->performance_tolerance.steadyStateError_pcntTol = STEADY_STATE_ERROR_TOLERANCE;
	}


	/* The outer loop will cycle through the many configured options */
	for (int trialNum = 0; trialNum < MAX_TRIALS; trialNum++)
	{
		/* Update any variables that need updating */
		for (int opNum = 0; opNum < NUM_OPTIMIZERS; opNum++)
		{
			trial_convergence_goals.data()[opNum]->populationSize = trial_popSizes[trialNum];
		}


		
		GAEngineStatistics_Vec avgStatistics;
		 
		/* The inner loop runs a trial many times to get an average performance number */
		auto avg_start_time = boost::chrono::high_resolution_clock::now();
		for (int avgCycle = 0; avgCycle < AVERAGING_CYCLES; avgCycle++)
		{ 
			GAEngine_PID engine(&avgStatistics);

			/* Give names to the optimizers being created */
			std::string optimizer1("OP1");
			std::string optimizer2("OP2");
			std::string optimizer3("OP3");

			engine.createOptimizer(optimizer1, trial_steps.data()[0]);
			engine.createOptimizer(optimizer2, trial_steps.data()[1]);
			engine.createOptimizer(optimizer3, trial_steps.data()[2]);

			GAModel_sPtr SS1 = std::make_shared<SSModel>();
			GAModel_sPtr SS2 = std::make_shared<SSModel>();
			GAModel_sPtr SS3 = std::make_shared<SSModel>();


			/*-------------------------
			* Attach everything to the engine
			*-------------------------*/
			engine.attachModel(optimizer1, SS1, trial_eval_models.data()[trialNum], SSModelName, USE_CPU);
			engine.attachObjective(optimizer1, trial_pid_goals.data()[0], std::string("PID1"));
			engine.attachConvergence(optimizer1, trial_convergence_goals.data()[0], std::string("CONVG1"));

			engine.attachModel(optimizer2, SS2, trial_eval_models.data()[trialNum], SSModelName, USE_CPU);
			engine.attachObjective(optimizer2, trial_pid_goals.data()[1], std::string("PID2"));
			engine.attachConvergence(optimizer2, trial_convergence_goals.data()[1], std::string("CONVG2"));

			engine.attachModel(optimizer3, SS3, trial_eval_models.data()[trialNum], SSModelName, USE_CPU);
			engine.attachObjective(optimizer3, trial_pid_goals.data()[2], std::string("PID3"));
			engine.attachConvergence(optimizer3, trial_convergence_goals.data()[2], std::string("CONVG3"));

			/*-------------------------
			* Execute in either serial or parallel
			*-------------------------*/
			engine.initAll();

			auto start_time = boost::chrono::high_resolution_clock::now();
			#if defined(GA_CPU_SINGLE_THREADED) && !defined(GA_CPU_MULTI_THREADED)
			engine.run(optimizer1, trialNum);
			engine.join(optimizer1);

			engine.run(optimizer2, trialNum);
			engine.join(optimizer2);

			engine.run(optimizer3, trialNum);
			engine.join(optimizer3);
			#endif

			#if defined(GA_CPU_MULTI_THREADED) && !defined(GA_CPU_SINGLE_THREADED)
			engine.runAll(trialNum);
			engine.joinAll();
			#endif
			auto end_time = boost::chrono::high_resolution_clock::now();
			boost::chrono::duration<double> total_time = boost::chrono::duration_cast<boost::chrono::duration<double>>(end_time - start_time);
			std::cout << "Total Execution Time: " << total_time << std::endl;

			engine.reportResultsToConsole();
			engine.reportResultsToFile(trialNum, total_time.count());
		}
		auto avg_end_time = boost::chrono::high_resolution_clock::now();
		boost::chrono::duration<double> avg_total_time = boost::chrono::duration_cast<boost::chrono::duration<double>>(avg_end_time - avg_start_time);
		
		reportAvgResultsToFile(trialNum, avg_total_time.count(), avgStatistics);
		

		std::cout << "FINISHED" << std::endl;

	}

	return 0;
}

void reportAvgResultsToFile(int trialNum, double avg_execution_time, GAEngineStatistics_Vec statisticsData)
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
	for (int i = 0; i < statisticsData.size(); i++)
		avgRunTime += statisticsData.data()[i].totalRunTime.count();

	avgRunTime /= (double)statisticsData.size();

	csvFile << "AVG Runtime:," << std::to_string(avgRunTime) << "\n";


	/*-----------------------------------------------
	* Best Fit
	*-----------------------------------------------*/
	double avgFitness = 0.0;
	for (int i = 0; i < statisticsData.size(); i++)
		avgFitness += statisticsData.data()[i].topPerformer.global_fitness;

	avgFitness /= (double)statisticsData.size();

	csvFile << "AVG Fitness:," << std::to_string(avgFitness) << "\n";
}