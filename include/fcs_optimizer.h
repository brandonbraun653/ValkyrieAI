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
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

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
	GA_METHOD_Sorting				sortingType = GA_SORT_DEFAULT;	

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


typedef boost::accumulators::accumulator_set<
	double,
	boost::accumulators::features<
	boost::accumulators::tag::mean,
	boost::accumulators::tag::variance>> StatisticsType;

typedef boost::container::vector<FCSOptimizer_PopulationMember> PopulationType;

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
	~FCSOptimizer() = default;

private:
	/*-----------------------------
	* Runtime Flags
	*----------------------------*/
	int currentGeneration;
	int refreshCounter;
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
		boost::container::vector<GA_SortBase_sPtr> sortingInstances;			
	} runtimeStep;
	
	FCSOptimizer_RuntimeConfig currentSolverParam;	/* Keeps track of the current execution style implemented */

	/* A highly generic description of the population members */
	boost::container::vector<FCSOptimizer_PopulationMember> parents, children;

	FCSOptimizer_Init_t settings;


	/*-----------------------------
	* Runtime Processing Data
	*----------------------------*/
	/* Parent Selections */
	boost::container::vector<int> fcs_parentSelections;

	/* Bred Chromosomes */
	GA_BreedingDataOutput fcs_bredChromosomes;




	/* Log of best fitness performance data for each generation.
	Do not pre-allocate memory as it is "pushed back" to add new data. */
	boost::container::vector<PID_FitnessScores> fcs_generationalBestFitnessData;
	
	
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
	void initStatistics();

	/*-----------------------------
	* Primary Algorithm Functions
	*----------------------------*/
	void checkConvergence();
	void evaluateModel(PopulationType& population);
	void evaluateFitness(PopulationType& population);
	void filterPopulation();
	void selectParents(PopulationType& population);
	void breedGeneration(PopulationType& population);
	void mutateGeneration(PopulationType& population);
	void boundaryCheck(PopulationType& population);
	void enforceResolution(PopulationType& population);
	void enforceTunerLimits(PopulationType& population);

	PopulationType sortPopulation(PopulationType* parents, PopulationType* children);

	bool algorithmNeedsNewSteps(/*statistical_data*/);
	void updateAlgorithmSteps(/*options, currentAlgSteps*/);

	/*-----------------------------
	* Data Collection and Reporting
	*----------------------------*/
	void gatherStatisticalData();
	void logData();
	

	struct GlobalStatistics
	{
		StatisticsType Kp;
		StatisticsType Ki;
		StatisticsType Kd;
	} GlobalChromStats;

	struct GenerationalStatistics
	{
		boost::container::vector<StatisticsType> Kp;
		boost::container::vector<StatisticsType> Ki;
		boost::container::vector<StatisticsType> Kd;
	} GenerationalChromStats;

	struct ChromosomeFrequency
	{
		boost::container::vector<int> Kp;
		boost::container::vector<int> Ki;
		boost::container::vector<int> Kd;
	} ChromOccurance;

	struct IdealStats
	{
		double Kp;
		double Ki;
		double Kd;
	} IdealChromVariance;

	double avgFitness = 0.0;
	boost::container::vector<double> averageFitness;	/* Creates snapshots of the average fitness at each generation */
	boost::container::vector<double> averagePOS;
	boost::container::vector<double> averageSSER;
	boost::container::vector<double> averageTS;
	boost::container::vector<double> averageTR;
};


#endif /* !GA_H_ */