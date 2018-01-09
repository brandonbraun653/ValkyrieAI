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
#include <boost/any.hpp>
#include <boost/function.hpp>
#include <boost/thread.hpp>
#include <boost/chrono.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/interprocess/ipc/message_queue.hpp>

/* Local Includes */
#include "ga_config.h"
#include "ga_steps.h"
#include "rng.hpp"
#include "model.h"
#include "logger.h"
#include "types.h"

/* Forward Declarations */
class FCSOptimizer;
typedef boost::shared_ptr<FCSOptimizer> FCSOptimizerClass_sPtr;
typedef boost::shared_ptr<boost::thread> Thread_sPtr;

/*-----------------------------------------------
* Useful Enumerations 
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

enum GA_RNG_Engine
{
	GA_MERSENNE_TWISTER
};

enum GA_RNG_Distribution
{
	GA_DISTRIBUTION_UNIFORM_REAL,
	GA_DISTRIBUTION_UNIFORM_INT
};

enum GA_METHOD_Breed
{
	GA_BREED_SIMPLE_CROSSOVER,
	GA_BREED_DYNAMIC_CROSSOVER,
	GA_BREED_FIXED_RATIO_CROSSOVER,
	GA_BREED_SIMULATED_BINARY_CROSSOVER,
	GA_BREED_TOTAL_OPTIONS,

	GA_BREED_SIMPLE_CROSSOVER_MSK = (1u << GA_BREED_SIMPLE_CROSSOVER),
	GA_BREED_DYNAMIC_CROSSOVER_MSK = (1u << GA_BREED_DYNAMIC_CROSSOVER),
	GA_BREED_FIXED_RATIO_CROSSOVER_MSK = (1u << GA_BREED_FIXED_RATIO_CROSSOVER),
	GA_BREED_SIMULATED_BINARY_CROSSOVER_MSK = (1u << GA_BREED_SIMULATED_BINARY_CROSSOVER)
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

struct FCSOptimizer_RuntimeConfig
{
	GA_METHOD_Breed					breedType;					/* Breed/Mate step operation */
	GA_METHOD_PopulationFilter		filterType;					/* Population Filtering step operation */
	GA_METHOD_ParentSelection		selectType;					/* Parent Selection step operation */
	GA_METHOD_MutateProbability		mutateProbabilityType;		/* Mutation step operation */
	GA_METHOD_MutateType			mutateType;					/* Mutation step operation */
	GA_METHOD_FitnessEvaluation		fitnessType;				/* Fitness Evaluation step operation */
	GA_METHOD_Resolution			resolutionType;				/* Sets mathematical resolution for all chromosomes */

	GA_RNG_Engine					rngEngine;					/* Random Number Generator Engine type */
	GA_RNG_Distribution				rngDistribution;			/* Random Number Generator Distribution type */
};


/**
* \brief Constrains the choices for the real time reconfiguration
* //TODO: Add a better description once this idea is fully fleshed out 
*/
struct FCSOptimizer_AllowedRuntimeConfig
{
	uint32_t allowedBreedTypes;
};

enum FCSOptimizer_Commands
{
	START,
	PAUSE,
	STOP,
	REPORT_STATUS
};

/*-----------------------------------------------
* Input/Output/Referencing Containers 
*-----------------------------------------------*/
struct FCSOptimizer_Init_t
{
	std::string optimizerName;								/* This gives a user friendly name to the optimizer */

	std::string resultsPath;								/* Path to directory for results reporting */

	std::string messageQueueName;							/* Name used to create a message queue between the main thread and an optimizer thread */


	SSModel_sPtr stateSpaceModel;							/* Possible reference to a State Space Model implementation. Leave empty if unused. */

	NNModel_sPtr neuralNetModel;							/* Possible reference to a TensorFlow Neural Network Graph. Leave empty if unused.  */


	ControlResponseJargon responseFeel;						/* A non-engineering way of describing how the tuner should optimize */


	PID_ControlSettings pidControlSettings;					/* A full engineering description of desired PID controller performance and limitations */


	FCSOptimizer_BasicConstraints basicConvergenceParam;	/* Simplified global convergence parameters */

	FCSOptimizer_AdvConstraints advConvergenceParam;		/* Convergence parameters that allow tuning how the underlying Genetic Algorithm software executes */

	FCSOptimizer_RuntimeConfig solverParam;					/* Configures how the optimizer internals operate at each step, to be updated dynamically by
															the internal reconfiguration engine */

	FCSOptimizer_AllowedRuntimeConfig solverParamOptions;	/* Possible options to select from for each step of the algorithm */
};

struct FCSOptimizer_Output_t
{
	//TODO: Put all the possible output data needed in here

	/* All this data is technically being passed by another thread, so it is going to need a mutex */
};

struct FCSOptimizer_Handle_t
{
	FCSOptimizer_Init_t Init;								/* Initialization settings for the engine */

	FCSOptimizer_Output_t Output;							/* Output data metrics for the optimization run */

	FCSOptimizerClass_sPtr Engine;							/* Instance of an optimization engine */

	Thread_sPtr Thread;										/* Reference to the thread running the optimizerEngine */

	GA_Status Status;										/* Current status of the tuner algorithm */

	boost::interprocess::message_queue* CommandQueue;		/* Interface to send commands to the optimizer thread */
};
typedef boost::shared_ptr<FCSOptimizer_Handle_t> FCSOptimizer_Handle;


//////////////////////////////////////////////////////////////////
/* Helper Functions */
//////////////////////////////////////////////////////////////////
extern void calculateMappingCoefficients(FCSOptimizer_MappingCoeff *mapping, double lower, double upper);
extern double enforceResolution(double in, GA_METHOD_Resolution res);


//////////////////////////////////////////////////////////////////
/* CLASS: FCSOptimizer */
//////////////////////////////////////////////////////////////////
class FCSOptimizer
{
public:
	/** 
	*	\brief Initialize the optimizer according to input settings 
	*/
	void init(FCSOptimizer_Init_t initializationSettings);

	/**
	*	\brief Run the optimizer until a solution is found or convergence criteria are met
	*/
	void run();

	/**
	*	\brief External request for results
	* The function is designed with the intent of being called from an instance of the
	* Valkyrie Engine to report result data back into an optimizer handle
	*/
	void requestOutput(FCSOptimizer_Output_t& output);

	FCSOptimizer();
	~FCSOptimizer();

private:
	/*-----------------------------
	* Runtime Flags
	*----------------------------*/
	int currentIteration;
	GA_Status currentStatus;

	/*-----------------------------
	* Threading Related Variables
	*----------------------------*/
	boost::mutex* print_to_console_mutex;
	boost::interprocess::message_queue* commandQueue;

	/*-----------------------------
	* Shared Objects 
	*----------------------------*/

	/**
	* Custom random number generators for PID values that generate values 
	* within the upper/lower bounds of user specified pidControlSettings.tuningLimits
	*/
	struct RNGEngines
	{
		RNGManager_sPtr Kp;
		RNGManager_sPtr Ki;
		RNGManager_sPtr Kd;
	} PID_RNG;
	


	boost::container::vector<FCSOptimizer_PopulationMember> population;
	FCSOptimizer_Init_t settings;

	/*-----------------------------
	* Runtime Processing Data
	*----------------------------*/
	/* Step Performance */

	/* Fitness Values */

	/* Parent Selections */

	/* Bred Chromosomes */

	boost::container::vector<GA_BreedBase_sPtr> breedInstances;								/* Implementations of unique breeding approaches */
	boost::function<void(GA_BreedingDataInput, GA_BreedingDataOutput&)> breedFunction;		/* Reference to the current breed function (fast access)*/


	/*-----------------------------
	* Constants for Mapping Conversions
	*-----------------------------*/
	FCSOptimizer_MappingCoeff mapCoefficients_Kp;
	FCSOptimizer_MappingCoeff mapCoefficients_Ki;
	FCSOptimizer_MappingCoeff mapCoefficients_Kd;

	/*-----------------------------
	* Setup Functions
	*----------------------------*/
	void initMemory();
	void initRNG();
	void initModel();
	void initPopulation();

	/*-----------------------------
	* Primary Algorithm Functions
	*----------------------------*/
	void checkConvergence();
	void evaluateModel();
	void evaluateFitness();
	void filterPopulation();
	void selectParents();
	void breedGeneration();
	void mutateGeneration();
	
	/*-----------------------------
	* Useful Helper Functions
	*----------------------------*/
	void enforceResolution();
	int reportResults(int trialNum);
	void printResultHighlights(double best_fit, int best_fit_idx);
};


#endif /* !GA_H_ */