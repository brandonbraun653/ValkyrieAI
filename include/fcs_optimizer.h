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
#include "config.h"
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

struct FCSOptimizer_RuntimeConfig
{
	GA_ChromosomeMappingType		chromType = MAPPING_TYPE_REAL;							/* Global data type used for chromosomes */

	GA_METHOD_Breed					breedType = GA_BREED_DEFAULT;							/* Breed/Mate step operation */
	GA_METHOD_PopulationFilter		filterType = GA_POPULATION_FILTER_DEFAULT;				/* Population Filtering step operation */
	GA_METHOD_ParentSelection		selectType = GA_SELECT_DEFAULT;							/* Parent Selection step operation */
	GA_METHOD_MutateProbability		mutateProbabilityType = GA_MUTATE_PROBABILITY_DEFAULT;	/* Mutation step operation */
	GA_METHOD_MutateType			mutateType = GA_MUTATE_DEFAULT;							/* Mutation step operation */
	GA_METHOD_ModelEvaluation		modelType;												/* Model Evaluation step operation */
	GA_METHOD_FitnessEvaluation		fitnessType = GA_FITNESS_DEFAULT;						/* Fitness Evaluation step operation */
	GA_METHOD_Resolution			resolutionType = GA_RESOLUTION_DEFAULT;					/* Sets mathematical resolution for all chromosomes */

	GA_RNG_Engine					rngEngine;												/* Random Number Generator Engine type */
	GA_RNG_Distribution				rngDistribution;										/* Random Number Generator Distribution type */
};


/**
* \brief Constrains the choices for the real time reconfiguration
* //TODO: Add a better description once this idea is fully fleshed out 
*/
struct FCSOptimizer_AllowedRuntimeConfig
{
	uint32_t allowedBreedTypes;
	uint32_t allowedFitnessTypes;
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

	SS_ModelBase_sPtr stateSpaceModel;						/* Possible reference to a State Space Model implementation. Leave empty if unused. */

	//NNModel_sPtr neuralNetModel;							/* Possible reference to a TensorFlow Neural Network Graph. Leave empty if unused.  */

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
	* Custom structures for passing information around to internal functions. These
	* essentially work as glue and help keep the actual algorithm as abstract as possible.
	* Each structure corresponds to the output data for each step and should only reference
	* a single population member. See "Runtime Processing Data" section below for actual 
	* implementation of this concept. 
	*----------------------------*/

	struct FCSOptimzer_SimulationData
	{
		/* Generic information about the results of the simulation step. Not all of 
		 * these fields will be used depending on the type of simulation that has been
		 * performed by the user.
		 */
		GA_METHOD_ModelEvaluation modelType;

		StateSpaceModelOutput ss_output;
		NeuralNetworkModelOutput nn_output;

		/* Add more support here as necessary */
	};

	struct FCSOptimizer_FitnessData
	{
		GA_METHOD_ModelEvaluation modelType;

		PID_FitnessScores fitness;
	};

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


	struct RuntimeApproaches
	{
		boost::container::vector<GA_PopulationFilterBase_sPtr> populationFilterInstances;	/* Implementations of unique population filtering approaches */
		boost::container::vector<GA_EvaluateModelBase_sPtr> evaluateModelInstances;			/* Implementations of unique model evaluation approaches */
		boost::container::vector<GA_EvaluateFitnessBase_sPtr> evaluateFitnessInstances;		/* Implementations of unique fitness evaluation approaches */
		boost::container::vector<GA_SelectParentBase_sPtr> selectParentInstances;			/* Implementations of unique parent selection approaches */
		boost::container::vector<GA_BreedBase_sPtr> breedInstances;							/* Implementations of unique breeding approaches */
		boost::container::vector<GA_MutateBase_sPtr> mutateInstances;						/* Implementations of unique mutation approaches */
	} runtimeStep;
	
	FCSOptimizer_RuntimeConfig currentSolverParam;	/* Keeps track of the current execution style implemented */

	/* A highly generic description of the population members */
	boost::container::vector<FCSOptimizer_PopulationMember> population;


	FCSOptimizer_Init_t settings;


	/*-----------------------------
	* Runtime Processing Data
	*----------------------------*/
	/* Model Evaluation Data */
	boost::container::vector<FCSOptimzer_SimulationData> fcs_modelEvalutationData;

	/* Fitness Data */
	boost::container::vector<FCSOptimizer_FitnessData> fcs_fitnessData;

	/* Parent Selections */
	boost::container::vector<int> fcs_parentSelections;

	/* Bred Chromosomes */
	GA_BreedingDataOutput fcs_bredChromosomes;

	/* Mutated Chromosomes */



	/* Log of best fitness performance data for each generation.
	Do not pre-allocate memory as it is "pushed back" to add new data. */
	boost::container::vector<FCSOptimizer_FitnessData> fcs_generationalBestFitnessData;
	
	
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
	void boundaryCheck();
	
	/*-----------------------------
	* Useful Helper Functions
	*----------------------------*/
	void enforceResolution();
	void enforceTunerLimits();
	int reportResults(int trialNum);
	void printResultHighlights(double best_fit, int best_fit_idx);
};


#endif /* !GA_H_ */