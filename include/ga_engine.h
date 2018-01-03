#pragma once
#ifndef GA_ENGINE_H_
#define GA_ENGINE_H_

/* C/C++ Includes */
#include <stdlib.h>
#include <stdio.h>

/* Boost Includes */
#include <boost/thread.hpp>
#include <boost/container/vector.hpp>

/* Local Includes */
#include "data.h"
#include "fcs_optimizer.h"
#include "model.h"


struct FCSOptimizer_Init_t //TODO: Rename this struct 
{
	std::string optimizerName;			/* Use this to give a user friendly unique name to the optimizer */
	
	SSModel_sPtr stateSpaceModel;		/* Possible reference to a State Space Model implementation */
	NNModel_sPtr neuralNetModel;		/* Possible reference to a TensorFlow Neural Network Graph  */
	
	

	PID_ControlGoals_sPtr pid;
	GA_ConverganceCriteria_sPtr ga_convergence;
	
};

struct FCSOptimizer_Handle_t
{
	FCSOptimizer_Init_t Init;

	FCSOptimizer_sPtr optimizerEngine;				/* Instance of an optimization engine */
	boost::shared_ptr<boost::thread> threadHandle;	/* Reference to the thread running the optimizerEngine */
};
typedef boost::shared_ptr<FCSOptimizer_Handle_t> FCSOptimizer_Handle;


class ValkyrieEngine
{
public:
	/*-------------------------------
	* Initialization and Setup
	*------------------------------*/
	FCSOptimizer_Handle newOptimizer(FCSOptimizer_Init_t init);


	/*-------------------------------
	* Basic Runtime Flow Control
	*------------------------------*/
	void init(FCSOptimizer_Handle optimizer);
	void run(FCSOptimizer_Handle optimizer);
	void pause(FCSOptimizer_Handle optimizer);
	void stop();

	void initAll();
	void runAll();
	void pauseAll();
	void stopAll();

	void waitForCompletion();


private:
	boost::container::vector<FCSOptimizer_Handle> optimizerInstances;

};


class GAEngine_PID
{
public:
	void createOptimizer(std::string optimizer, GA_AlgorithmMethods alg_step_types);
	void attachModel(std::string optimizer, GAModel_sPtr model_type, SS_NLTIV_ModelBase_sPtr model_system, std::string model_name, int processor);
	void attachObjective(std::string optimizer, PID_ControlGoals_sPtr pid_parameters, std::string pid_name);
	void attachConvergence(std::string optimizer, GA_ConverganceCriteria_sPtr convergence_param, std::string convergence_name);

	void join(std::string optimizer);
	void joinAll();

	void init(std::string optimizer);
	void initAll();

	void run(std::string optimizer, int trialNum);
	void runAll(int trialNum);

	void reportResultsToConsole();
	void reportResultsToFile(int trialNum, double full_execution_time);
	void reportAvgResultsToFile(int trialNum, double avg_execution_time);

	GAEngine_PID(GAEngineStatistics_Vec *statisticsData);
	~GAEngine_PID();

private:
	boost::mutex printMutex;
	
	boost::mutex statisticsMutex;
	GAEngineStatistics_Vec threadStatistics;
	GAEngineStatistics_Vec* avgThreadStatistics;

	std::map<std::string, uint32_t> thread_mapping;
	std::vector< Thread_sPtr > thread;

	/*----------------------------------------------
	* Containers for handling an arbitrary number of optimizers.
	* Note: shared_ptr chosen to prevent memory leaks
	*----------------------------------------------*/
	std::vector< GAModel_sPtr > model;
	std::vector< PID_ControlGoals_sPtr > pid;
	std::vector< FCSOptimizer_sPtr > ga_optimizer;
	std::vector< GA_ConverganceCriteria_sPtr> ga_convergence;

	/*----------------------------------------------
	* Containers for various kinds of meta-data. Usually these
	* are for user-friendly reporting of performance characteristics.
	*----------------------------------------------*/
	std::vector< std::string > optimizer_names;
	std::vector< std::string > convergence_names;
	std::vector< std::string > pid_names;
	std::vector< std::string > model_names;

	void assignOptimizer_name(std::string optimizer, std::string optimizer_name);
	void assignOptimizer_model(std::string optimizer, GAModel_sPtr model, SS_NLTIV_ModelBase_sPtr model_system, int processor);
	void assignOptimizer_objective(std::string optimizer, PID_ControlGoals_sPtr objective);
	void assignOptimizer_convergence(std::string optimizer, GA_ConverganceCriteria_sPtr convergence);
};

#endif /* !GA_ENGINE_H_*/