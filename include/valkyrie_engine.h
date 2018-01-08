#pragma once
#ifndef GA_ENGINE_H_
#define GA_ENGINE_H_

/* C/C++ Includes */
#include <stdlib.h>
#include <stdio.h>

/* Boost Includes */
#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/container/vector.hpp>
#include <boost/interprocess/ipc/message_queue.hpp>

/* Local Includes */
#include "types.h"
#include "fcs_optimizer.h"
#include "model.h"

/**
* Management engine for handling initialization, running, and collecting
* performance results from an arbitrary number of FCSOptimizers.
**/
class ValkyrieEngine
{
public:
	/** 
	* \brief Creates an optimizer
	* Create a new Flight Control System (FCS) tuner instance that utilizes
	* artificial intelligence to try and quickly tune a quadcopter control 
	* loop.
	* 
	* \param init 
	* Initialization struct for the optimizer
	* \return 
	* Reference handle to the newly created object
	*/
	FCSOptimizer_Handle newOptimizer(FCSOptimizer_Init_t init);

	/**
	* \brief Initialize an FCSOptimizer
	* Executes the initialization function tied to the optimizer instance
	*
	* \param optimizer
	* Reference to an optimizer handle
	* 
	* \return void
	*/
	void initialize(const FCSOptimizer_Handle& optimizer);

	/**
	* \brief Start an FCSOptimizer
	* Begin the main thread of the optimizer instance that performs a single 
	* round of Flight Control System (FCS) tuning.
	*
	* \param optimizer
	* Reference to an optimizer handle
	*
	* \return void
	*/
	void start(const FCSOptimizer_Handle& optimizer);

	/**
	* \brief Pause an FCSOptimizer
	* Pause (block) execution of an optimizer instance
	*
	* \param optimizer
	* Reference to an optimizer handle
	*
	* \return void
	*/
	void pause(const FCSOptimizer_Handle& optimizer);

	/**
	* \brief Stop an FCSOptimizer
	* Stop execution of an optimizer instance and destroy all progress made.
	* This irreversibly halts the optimizer thread and exits upon resource 
	* cleanup. 
	*
	* \param optimizer
	* Reference to an optimizer handle
	*
	* \return void
	*/
	void stop(const FCSOptimizer_Handle& optimizer);

	/**
	* \brief Initialize all valid instances of an FCSOptimizer
	*
	* \return void
	*/
	void initAll();

	/**
	* \brief Start execution of all valid instances of an FCSOptimizer
	*
	* \return void
	*/
	void startAll();

	/**
	* \brief Pause execution of all valid instances of an FCSOptimizer
	*
	* \return void
	*/
	void pauseAll();

	/**
	* \brief Halt and exit all valid instances of an FCSOptimizer
	*
	* \return void
	*/
	void stopAll();

	/**
	* \brief Wait for all valid instances of an FCSOptimizer to complete
	*
	* \return void
	*/
	void waitForCompletion();

	ValkyrieEngine();
	~ValkyrieEngine();

private:
	boost::container::vector<FCSOptimizer_Handle> optimizerInstances;

	/*-----------------------------
	* Resource Protection 
	*----------------------------*/

	boost::shared_ptr<boost::mutex> ostream_mtx_sPtr;	/* Print to console "cout" mutex */

};

// 
// class GAEngine_PID
// {
// public:
// 	void createOptimizer(std::string optimizer, GA_AlgorithmMethods alg_step_types);
// 	void attachModel(std::string optimizer, GAModel_sPtr model_type, SS_NLTIV_ModelBase_sPtr model_system, std::string model_name, int processor);
// 	void attachObjective(std::string optimizer, PID_ControlSettings_sPtr pid_parameters, std::string pid_name);
// 	void attachConvergence(std::string optimizer, FCSOptimizer_AdvConstraints_sPtr convergence_param, std::string convergence_name);
// 
// 	void join(std::string optimizer);
// 	void joinAll();
// 
// 	void init(std::string optimizer);
// 	void initAll();
// 
// 	void run(std::string optimizer, int trialNum);
// 	void runAll(int trialNum);
// 
// 	void reportResultsToConsole();
// 	void reportResultsToFile(int trialNum, double full_execution_time);
// 	void reportAvgResultsToFile(int trialNum, double avg_execution_time);
// 
// 	GAEngine_PID(GAEngineStatistics_Vec *statisticsData);
// 	~GAEngine_PID();
// 
// private:
// 	boost::mutex printMutex;
// 	
// 	boost::mutex statisticsMutex;
// 	GAEngineStatistics_Vec threadStatistics;
// 	GAEngineStatistics_Vec* avgThreadStatistics;
// 
// 	std::map<std::string, uint32_t> thread_mapping;
// 	std::vector< Thread_sPtr > thread;
// 
// 	/*----------------------------------------------
// 	* Containers for handling an arbitrary number of optimizers.
// 	* Note: shared_ptr chosen to prevent memory leaks
// 	*----------------------------------------------*/
// 	std::vector< GAModel_sPtr > model;
// 	std::vector< PID_ControlSettings_sPtr > pid;
// 	std::vector< FCSOptimizerClass_sPtr > ga_optimizer;
// 	std::vector< FCSOptimizer_AdvConstraints_sPtr> ga_convergence;
// 
// 	/*----------------------------------------------
// 	* Containers for various kinds of meta-data. Usually these
// 	* are for user-friendly reporting of performance characteristics.
// 	*----------------------------------------------*/
// 	std::vector< std::string > optimizer_names;
// 	std::vector< std::string > convergence_names;
// 	std::vector< std::string > pid_names;
// 	std::vector< std::string > model_names;
// 
// 	void assignOptimizer_name(std::string optimizer, std::string optimizer_name);
// 	void assignOptimizer_model(std::string optimizer, GAModel_sPtr model, SS_NLTIV_ModelBase_sPtr model_system, int processor);
// 	void assignOptimizer_objective(std::string optimizer, PID_ControlSettings_sPtr objective);
// 	void assignOptimizer_convergence(std::string optimizer, FCSOptimizer_AdvConstraints_sPtr convergence);
// };

#endif /* !GA_ENGINE_H_*/