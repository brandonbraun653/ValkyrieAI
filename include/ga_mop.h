#ifndef GA_H_
#define GA_H_

/* C/C++ Includes */
#include <iostream>
#include <random>
#include <algorithm>
#include <numeric>
#include <windows.h>
#include "stdint.h"
#include "math.h"

/* Boost Includes */
#include <boost/thread.hpp>
#include <boost/chrono.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/pointer_cast.hpp>  //TODO: Remove header when finished with removal of old model referencing 

/* Local Includes */
#include "ga_config.h"
#include "ga_mop_steps.h"
#include "model.h"
#include "logger.h"
#include "data.h"
#include "host_memory.h"

/*-----------------------------------------------
* Global Genetic Algorithm "Stuff"
*-----------------------------------------------*/
struct GA_AlgorithmMethods
{
	GA_METHOD_Breed breedType;
	GA_METHOD_PopulationFilter filterType;
	GA_METHOD_ParentSelection selectType;
	GA_METHOD_MutateProbability mutateProbabilityType;
	GA_METHOD_MutateType mutateType;
	GA_METHOD_FitnessEvaluation fitnessType;
	GA_METHOD_Resolution resolutionType;
};

//////////////////////////////////////////////////////////////////
/* Helper Functions */
//////////////////////////////////////////////////////////////////
extern void calculateMappingCoefficients(mapCoeff_t *mapping, double lower, double upper);

//////////////////////////////////////////////////////////////////
/* CLASS: GA_MOP */
//////////////////////////////////////////////////////////////////
class GA_MOP
{
public:
	void init();
	void reset();

	void run(boost::mutex* resultsMutex, GAEngineStatistics_Vec* resultsStatistics, 
		GAEngineStatistics_Vec* avgResultsStatistics, int threadIndex, int trialNum);

	/* Configuration Functions */
	void registerOutputMutex(boost::mutex* printMutex);
	void registerName(std::string optimizerName);
	void registerObjective(PID_ControlGoals_sPtr pidGoals, std::string pidName);
	void registerModel(GAModel_sPtr modelType, SS_NLTIV_ModelBase_sPtr model, std::string name, int processor);
	void registerConvergence(GA_ConverganceCriteria_sPtr convLimits, std::string convName);

	GA_MOP(GA_AlgorithmMethods alg_methods);
	~GA_MOP();

private:
	/*-----------------------------
	* User Input Data
	*----------------------------*/
	boost::mutex* print_to_console_mutex;

	int compute_processor;
	std::string ga_instance_optimizer_name,
		ga_instance_pid_name,
		ga_instance_model_name,
		ga_instance_convergence_name;

	/* Allows the user to specify how the algorithm runs at each step */
	GA_AlgorithmMethods ga_instance_step_methods;

	/* Contains the user's pid goals, constraints, and weightings */
	PID_ControlGoals_sPtr ga_instance_pid_config_data;

	/* Contains the user's desired convergence criteria */
	GA_ConverganceCriteria_sPtr ga_instance_convergence_criteria;

	/*-----------------------------
	* Runtime Processing Data
	*----------------------------*/
	GAMOP_hPID_Data hData;

	/* Step Performance */
	boost::mutex SS_StepPerformance_mutex;
	StepPerformance_Vec SS_StepPerformance;

	/* Fitness Values */
	PIDFitness_Vec SS_FitnessValues;
	PIDFitness_Vec GA_BestFitnessValues;
	PIDElitist GA_ElitistSolutions;

	/* Parent Selections */
	iVec parentSelection;

	/* Bred Chromosomes */
	hPID_Chromosomes bredChromosomes;

	/*-----------------------------
	* Runtime Flags
	*----------------------------*/
	bool optimizer_initialized;
	int currentIteration;
	GA_Status currentStatus;


	SSModel_sPtr modelSS;


	SS_NLTIV_ModelBase_sPtr ss_user_system_model;

	/*-----------------------------
	* Constants for Mapping Conversions
	*-----------------------------*/
	mapCoeff_t mapCoefficients_Kp;
	mapCoeff_t mapCoefficients_Ki;
	mapCoeff_t mapCoefficients_Kd;

	/*-----------------------------
	* Setup Functions
	*----------------------------*/
	void initMemory();
	void initModel();
	void initPopulation();
	void deInitMemory();

	/*-----------------------------
	* Primary Algorithm Functions
	*----------------------------*/
	void evaluateModel();
	void evaluateFitness();
	void filterPopulation();
	void selectParents();
	void breedGeneration();
	void mutateGeneration();
	void checkConvergence();

	/*-----------------------------
	* Useful Helper Functions
	*----------------------------*/
	void enforceResolution();
	int reportResults(int trialNum);
	void printResultHighlights(double best_fit, int best_fit_idx);

};
typedef std::shared_ptr<GA_MOP> GAMOP_sPtr;

#endif /* !GA_H_ */