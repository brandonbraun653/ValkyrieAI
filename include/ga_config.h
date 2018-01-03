#pragma once
#ifndef GA_CONFIG_H_
#define GA_CONFIG_H_

/*-----------------------------------------------
* Global Genetic Algorithm "Stuff"
*-----------------------------------------------*/
#define DEBUGGING_ENABLED
//#define GA_CPU_SINGLE_THREADED
#define GA_CPU_MULTI_THREADED

/*-----------------------------------------------
* Runtime Options
*-----------------------------------------------*/
/* Initialize Population Function */
//#define GA_TRACE_INITIALIZE_POPULATION

/* Evaluate Model Function */
#define GA_TRACE_EVALUATE_MODEL
//#define FIXED_PID_VALUES

/* Population Filter Function */
//#define GA_TRACE_FILTER_POPULATION

/* Evaluate Fitness Function */
//#define GA_TRACE_EVALUATE_FITNESS
//#define GA_TRACE_SELECT_PARENTS

/* Breed Generation Function */
#define GA_TRACE_BREED_GENERATION
#define GA_ENFORCE_RESOLUTION_BG

/* Mutate Generation Function */
#define GA_TRACE_MUTATE_GENERATION
#define GA_ENFORCE_RESOLUTION_MG

/* Check Convergence Function */
//#define GA_TRACE_CHECK_CONVERGENCE
#define GA_REPORT_DATA_CHECK_CONVERGENCE

/*-----------------------------------------------
* model.h
*-----------------------------------------------*/
/* Debugging Statements */
//#define SS_TRACE
//#define SS_TRACE_LOG_COUT
//#define SS_TRACE_LOG_CSV
//#define SS_THREAD_TRACE_LOG_COUT
//#define SS_THREAD_TRACE_LOG_CSV


/*-----------------------------------------------
* Useful Enums
*-----------------------------------------------*/
enum GA_RunMode
{
	SINGLE_THREADED,
	MULTI_THREADED,
	SINGLE_THREADED_WITH_CUDA,
	MULTI_THREADED_WITH_CUDA
};

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

enum GA_BreedMethod
{
	GA_BREED_SIMPLE_CROSSOVER,
	GA_BREED_DYNAMIC_CROSSOVER,
	GA_BREED_FIXED_RATIO_CROSSOVER,
	GA_BREED_SIMULATED_BINARY_CROSSOVER
};

enum GA_PopulationFilterMethod
{
	GA_POPULATION_STATIC_FILTER,
	GA_POPULATION_DYNAMIC_FILTER
};

enum GA_ParentSelectMethod
{
	GA_SELECT_RANDOM,
	GA_SELECT_RANKED,
	GA_SELECT_ROULETTE,
	GA_SELECT_STOCHASTIC_SAMPLING,
	GA_SELECT_TOURNAMENT,
	GA_SELECT_ELITIST
};

enum GA_MutateProbabilityMethod
{
	GA_MUTATE_PROBABILITY_POISSON,
	GA_MUTATE_PROBABILITY_EXPONENTIAL,
	GA_MUTATE_PROBABILITY_GAMMA,
	GA_MUTATE_PROBABILITY_WEIBULL,
	GA_MUTATE_PROBABILITY_CHI_SQUARED
};

enum GA_MutateMethod
{
	GA_MUTATE_BIT_FLIP,
	GA_MUTATE_ADD_SUB
};

enum GA_EvaluateFitnessMethod
{
	GA_FITNESS_WEIGHTED_SUM,
	GA_FITNESS_NON_DOMINATED_SORT
};

enum GA_Resolution
{
	// DP: Decimal Places
	GA_RESOLUTION_0DP,
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
#endif