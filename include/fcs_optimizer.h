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
enum GA_Status
{
	GA_IDLE,
	GA_OK,
	GA_READY,
	GA_HALT,
	GA_BUSY,
	GA_SETUP,
	GA_PAUSED,
	GA_INPROGRESS,
	GA_COMPLETE,
	GA_ERROR
};

enum GA_METHOD_Breed
{
	GA_BREED_SIMPLE_CROSSOVER,
	GA_BREED_DYNAMIC_CROSSOVER,
	GA_BREED_FIXED_RATIO_CROSSOVER,
	GA_BREED_SIMULATED_BINARY_CROSSOVER
};

enum GA_METHOD_PopulationFilter
{
	GA_POPULATION_STATIC_FILTER,
	GA_POPULATION_DYNAMIC_FILTER
};

enum GA_METHOD_ParentSelection
{
	GA_SELECT_RANDOM,
	GA_SELECT_RANKED,
	GA_SELECT_ROULETTE,
	GA_SELECT_STOCHASTIC_SAMPLING,
	GA_SELECT_TOURNAMENT,
	GA_SELECT_ELITIST
};

enum GA_METHOD_MutateProbability
{
	GA_MUTATE_PROBABILITY_POISSON,
	GA_MUTATE_PROBABILITY_EXPONENTIAL,
	GA_MUTATE_PROBABILITY_GAMMA,
	GA_MUTATE_PROBABILITY_WEIBULL,
	GA_MUTATE_PROBABILITY_CHI_SQUARED
};

enum GA_METHOD_MutateType
{
	GA_MUTATE_BIT_FLIP,
	GA_MUTATE_ADD_SUB
};

enum GA_METHOD_FitnessEvaluation
{
	GA_FITNESS_WEIGHTED_SUM,
	GA_FITNESS_NON_DOMINATED_SORT
};

enum GA_METHOD_Resolution
{
	GA_RESOLUTION_0DP, // DP: Decimal Places
	GA_RESOLUTION_1DP,
	GA_RESOLUTION_2DP,
	GA_RESOLUTION_3DP,
	GA_RESOLUTION_4DP,
	GA_RESOLUTION_5DP,
	GA_RESOLUTION_6DP,
	GA_RESOLUTION_7DP,
	GA_RESOLUTION_8DP,
	GA_RESOLUTION_9DP,
	GA_RESOLUTION_10DP
};

struct GA_RuntimeConfig
{
	GA_METHOD_Breed					breedType;
	GA_METHOD_PopulationFilter		filterType;
	GA_METHOD_ParentSelection		selectType;
	GA_METHOD_MutateProbability		mutateProbabilityType;
	GA_METHOD_MutateType			mutateType;
	GA_METHOD_FitnessEvaluation		fitnessType;
	GA_METHOD_Resolution			resolutionType;
};

//////////////////////////////////////////////////////////////////
/* Helper Functions */
//////////////////////////////////////////////////////////////////
extern void calculateMappingCoefficients(mapCoeff_t *mapping, double lower, double upper);

//////////////////////////////////////////////////////////////////
/* CLASS: FCSOptimizer */
//////////////////////////////////////////////////////////////////
class FCSOptimizer
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

	FCSOptimizer(GA_RuntimeConfig alg_methods);
	~FCSOptimizer();

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
	GA_RuntimeConfig ga_instance_step_methods;

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
typedef std::shared_ptr<FCSOptimizer> FCSOptimizer_sPtr;

#endif /* !GA_H_ */