#pragma once
#ifndef TYPES_H_
#define TYPES_H_

/* C/C++ Includes */
#include <stdlib.h>
#include <stdint.h>
#include <string>
#include <iostream>

/* Boost Includes */
#include <boost/shared_ptr.hpp>
#include <boost/container/vector.hpp>

/* Eigen Includes */
#include <eigen/Eigen>



/*-----------------------------------------------
* INPUT DATA
*-----------------------------------------------*/

/**
* \brief Non-Engineering perspective on tuner goals.
* Setting this single value in the Init field for an FCS Optimizer instance
* will apply a set of PID tuning performance goals that should approximate
* the described feeling.
*/
enum ControlResponseJargon
{
	SNAPPY,
	LOOSE,
	SAFE,
	TIGHT,
	RESPONSIVE,
	AGGRESSIVE,
	MOLASSES,
	SMOOTH,
	SHARP,
	DIALED_IN
};


/*-----------------------------------------------
* PID Data Types
*-----------------------------------------------*/

/** 
* \brief Engineering perspective of tuner goals for a SISO system
*/
struct PID_PerformanceGoals
{
	double percentOvershoot_goal = 0.0;	/* Units: none, absolute percentage -> 0.02 == 2% */
	double steadyStateError_goal = 0.0;	/* Units: none, absolute in whatever units the simulation output uses */
	double settlingTime_goal = 0.0;		/* Units: seconds */
	double riseTime_goal = 0.0;			/* Units: seconds */
};

/** 
* \brief Specified tolerance (% +/- of goal) for each parameter
*/
struct PID_PerformanceTolerance
{
	double percentOvershoot_pcntTol = 0.0;	/* Percentage -> 0.02 == 2% */
	double steadyStateError_pcntTol = 0.0;	/* Percentage -> 0.02 == 2% */
	double settlingTime_pcntTol = 0.0;		/* Percentage -> 0.02 == 2% */
	double riseTime_pcntTol = 0.0;			/* Percentage -> 0.02 == 2% */
};

/** 
* \brief Weights performance metrics according to user importance
* The FCS optimizer will use these values to help give it a direction
* to try and tune in. Each parameter can have a value in the range [0-1], 
* with 0 as unimportant and 1 as high importance.
*/
struct PID_PerformanceWeights
{
	double percentOvershoot_weight = 0.0;
	double steadyStateError_weight = 0.0;
	double settlingTime_weight = 0.0;
	double riseTime_weight = 0.0;
};

/** 
* \brief Holds a single set of Kp, Ki, Kd values 
*/
struct PID_Values
{
	double Kp = 0.0;
	double Ki = 0.0;
	double Kd = 0.0;
};

/** 
* \brief Specify the upper and lower Kp, Ki, Kd tuning boundaries
*/
struct PID_TuningLimits
{
	struct PID_Limits
	{
		double lower;
		double upper;
	};

	PID_Limits Kp, Ki, Kd;
};

/** 
* \brief All Possible FCS Optimizer PID Settings
* A conglomeration of various structs that fully characterize
* performance metrics and limits for a single PID control loop.
*/
struct PID_ControlSettings
{
	PID_PerformanceGoals		performanceGoals;
	PID_PerformanceTolerance	performanceTolerance;
	PID_PerformanceWeights		performanceWeights;
	PID_TuningLimits			tuningLimits;
	PID_Values					tuningInitialValues;
};


/*-----------------------------------------------
* Genetic Algorithm Data Types 
*-----------------------------------------------*/

template<typename ChromType>
struct GA_PIDChromosome
{
	ChromType Kp;
	ChromType Ki;
	ChromType Kd;
};

enum GA_ChromosomeMappingType
{
	MAPPING_TYPE_REAL,
	MAPPING_TYPE_BIT_FIELD
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

enum GA_RNG_Engine
{
	GA_MERSENNE_TWISTER
};

enum GA_RNG_Distribution
{
	GA_DISTRIBUTION_UNIFORM_REAL,
	GA_DISTRIBUTION_UNIFORM_INT
};

enum GA_METHOD_Sorting
{
	GA_SORT_FAST_NONDOMINATED,

	GA_SORT_TOTAL_OPTIONS,
	GA_SORT_DEFAULT = GA_SORT_FAST_NONDOMINATED,

	GA_SORT_FAST_NONDOMINATED_MSK = (1u << GA_SORT_FAST_NONDOMINATED)
};

enum GA_METHOD_Breed
{
	GA_BREED_SIMPLE_CROSSOVER,
	GA_BREED_DYNAMIC_CROSSOVER,
	GA_BREED_FIXED_POINT_CROSSOVER,
	GA_BREED_SIMULATED_BINARY_CROSSOVER,
	GA_BREED_TOTAL_OPTIONS,
	GA_BREED_DEFAULT = GA_BREED_SIMPLE_CROSSOVER,

	GA_BREED_SIMPLE_CROSSOVER_MSK = (1u << GA_BREED_SIMPLE_CROSSOVER),
	GA_BREED_DYNAMIC_CROSSOVER_MSK = (1u << GA_BREED_DYNAMIC_CROSSOVER),
	GA_BREED_FIXED_RATIO_CROSSOVER_MSK = (1u << GA_BREED_FIXED_POINT_CROSSOVER),
	GA_BREED_SIMULATED_BINARY_CROSSOVER_MSK = (1u << GA_BREED_SIMULATED_BINARY_CROSSOVER)
};

enum GA_METHOD_PopulationFilter
{
	GA_POPULATION_STATIC_FILTER,
	GA_POPULATION_DYNAMIC_FILTER,
	GA_POPULATION_TOTAL_OPTIONS,
	GA_POPULATION_FILTER_DEFAULT = GA_POPULATION_STATIC_FILTER,

	GA_POPULATION_STATIC_FILTER_MSK = (1u << GA_POPULATION_STATIC_FILTER),
	GA_POPULATION_DYNAMIC_FILTER_MSK = (1u << GA_POPULATION_DYNAMIC_FILTER)
};

enum GA_METHOD_ParentSelection
{
	GA_SELECT_RANDOM,
	GA_SELECT_RANKED,
	GA_SELECT_ROULETTE,
	GA_SELECT_STOCHASTIC_SAMPLING,
	GA_SELECT_TOURNAMENT,
	GA_SELECT_ELITIST,
	GA_SELECT_TOTAL_OPTIONS,
	GA_SELECT_DEFAULT = GA_SELECT_RANDOM,

	GA_SELECT_RANDOM_MSK = (1u << GA_SELECT_RANDOM),
	GA_SELECT_RANKED_MSK = (1u << GA_SELECT_RANKED),
	GA_SELECT_ROULETTE_MSK = (1u << GA_SELECT_ROULETTE),
	GA_SELECT_STOCHASTIC_SAMPLING_MSK = (1u << GA_SELECT_STOCHASTIC_SAMPLING),
	GA_SELECT_TOURNAMENT_MSK = (1u << GA_SELECT_TOURNAMENT),
	GA_SELECT_ELITIST_MSK = (1u << GA_SELECT_ELITIST)
};

enum GA_METHOD_MutateProbability
{
	GA_MUTATE_PROBABILITY_POISSON,
	GA_MUTATE_PROBABILITY_EXPONENTIAL,
	GA_MUTATE_PROBABILITY_GAMMA,
	GA_MUTATE_PROBABILITY_WEIBULL,
	GA_MUTATE_PROBABILITY_CHI_SQUARED,
	GA_MUTATE_PROBABILITY_TOTAL_OPTIONS,
	GA_MUTATE_PROBABILITY_DEFAULT = GA_MUTATE_PROBABILITY_EXPONENTIAL,

	GA_MUTATE_PROBABILITY_POISSON_MSK = (1u << GA_MUTATE_PROBABILITY_POISSON),
	GA_MUTATE_PROBABILITY_EXPONENTIAL_MSK = (1u << GA_MUTATE_PROBABILITY_EXPONENTIAL),
	GA_MUTATE_PROBABILITY_GAMMA_MSK = (1u << GA_MUTATE_PROBABILITY_GAMMA),
	GA_MUTATE_PROBABILITY_WEIBULL_MSK = (1u << GA_MUTATE_PROBABILITY_WEIBULL),
	GA_MUTATE_PROBABILITY_CHI_SQUARED_MSK = (1u << GA_MUTATE_PROBABILITY_CHI_SQUARED)
};

enum GA_METHOD_MutateType
{
	GA_MUTATE_BIT_FLIP,
	GA_MUTATE_ADD_SUB,
	GA_MUTATE_TOTAL_OPTIONS,
	GA_MUTATE_DEFAULT = GA_MUTATE_BIT_FLIP,

	GA_MUTATE_BIT_FLIP_MSK = (1u << GA_MUTATE_BIT_FLIP),
	GA_MUTATE_ADD_SUB_MSK = (1u << GA_MUTATE_ADD_SUB)
};

enum GA_METHOD_FitnessEvaluation
{
	GA_FITNESS_WEIGHTED_SUM,
	GA_FITNESS_MEAN_SQUARE_ERROR,
	GA_FITNESS_TOTAL_OPTIONS,
	GA_FITNESS_DEFAULT = GA_FITNESS_WEIGHTED_SUM,

	GA_FITNESS_WEIGHTED_SUM_MSK = (1u << GA_FITNESS_WEIGHTED_SUM),
	GA_FITNESS_MEAN_SQUARE_ERROR_MSK = (1u << GA_FITNESS_MEAN_SQUARE_ERROR)
};

enum GA_METHOD_ModelEvaluation
{
	GA_MODEL_STATE_SPACE,
	GA_MODEL_NEURAL_NETWORK,
	GA_MODEL_TOTAL_OPTIONS,

	GA_MODEL_STATE_SPACE_MSK = (1u << GA_MODEL_STATE_SPACE),
	GA_MODEL_NEURAL_NETWORK_MSK = (1u << GA_MODEL_NEURAL_NETWORK)
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
	GA_RESOLUTION_10DP,
	GA_RESOLUTION_DEFAULT = GA_RESOLUTION_3DP
};


/*-----------------------------------------------
* Model Data Types
*-----------------------------------------------*/

enum ModelSimulationType
{
	IMPULSE,
	STEP,
	RAMP,
	QUADRATIC,
	CUSTOM,
};

class SS_ModelBase
{
public:
	virtual Eigen::MatrixXd getA() = 0;
	virtual Eigen::MatrixXd getB() = 0;
	virtual Eigen::MatrixXd getC() = 0;
	virtual Eigen::MatrixXd getD() = 0;
	virtual Eigen::MatrixXd getX0() = 0;

	virtual int getNumInputs() = 0;
	virtual int getNumOutputs() = 0;
	virtual int getNumStates() = 0;

	virtual ~SS_ModelBase() = default;
private:
};
typedef boost::shared_ptr<SS_ModelBase> SS_ModelBase_sPtr;


class SS_NLTIVModel : public SS_ModelBase
{
public:
	Eigen::MatrixXd A, B, C, D, U, X0;
	int inputs, outputs, states;

	SS_NLTIVModel() = default;
	~SS_NLTIVModel() = default;

	/* Copy constructor */
	SS_NLTIVModel(const SS_ModelBase_sPtr& base)
	{
		//Note: Am I assuming too much by taking directly from the base class?
		inputs = base->getNumInputs();
		outputs = base->getNumOutputs();
		states = base->getNumStates();

		A = base->getA();
		B = base->getB();
		C = base->getC();
		D = base->getD();
		X0 = base->getX0();

		U.resize(inputs, 1); U.setZero(inputs, 1);
	}

	SS_NLTIVModel(const int Inputs, const int Outputs, const int States) : \
		inputs(Inputs), outputs(Outputs), states(States)
	{
		A.resize(states, states); A.setZero(states, states);
		B.resize(states, inputs); B.setZero(states, inputs);
		C.resize(outputs, states); C.setZero(outputs, states);
		D.resize(outputs, inputs); D.setZero(outputs, inputs);
		X0.resize(states, 1); X0.setZero(states, 1);
		U.resize(inputs, 1); U.setZero(inputs, 1);
	}

	/* System Function: compatible with Odeint solvers */
	void operator()(const Eigen::MatrixXd x, Eigen::MatrixXd &dxdt, double t)
	{
		dxdt = A*x + B*U;
	}


	Eigen::MatrixXd getA() override { return A; };
	Eigen::MatrixXd getB() override { return B; };
	Eigen::MatrixXd getC() override { return C; };
	Eigen::MatrixXd getD() override { return D; };
	Eigen::MatrixXd getX0() override { return X0; };
	
	int getNumInputs() override { return inputs; };
	int getNumOutputs() override { return outputs; };
	int getNumStates() override { return states; };

private:
};
typedef boost::shared_ptr<SS_NLTIVModel> SS_NLTIVModel_sPtr;


class NN_ModelBase
{
public:
	virtual int initialize() = 0;

	virtual ~NN_ModelBase() = default;
};
typedef boost::shared_ptr<NN_ModelBase> NN_ModelBase_sPtr;

/*-----------------------------------------------
* OUTPUT DATA
*-----------------------------------------------*/
enum PerformanceType
{
	PT_MEAN_SQUARED_ERROR,
	PT_NORMALIZED_MEAN_SQUARED_ERROR,
	PT_STEP,
	PT_RAMP,
	PT_QUADRATIC,
	PT_CUSTOM,
	PT_TOTAL_OPTIONS
};


struct MSEPerformance
{
	double POS = -1.0;		/* Percent Overshoot */
	double SSER = -1.0;		/* Steady State Error */
	double TS = -1.0;		/* Settling Time */
	double TR = -1.0;		/* Rise Time */
	double SSV = -1.0;		/* Steady State Value */
};
typedef boost::shared_ptr<MSEPerformance> MSEPerformance_sPtr;

struct NMSEPerformance
{
	double POS = -1.0;		/* Percent Overshoot */
	double SSER = -1.0;		/* Steady State Error */
	double dOS = -1.0;		/* delta Overshoot */
	double TS = -1.0;		/* Settling Time */
	double TR = -1.0;		/* Rise Time */
	double SSV = -1.0;		/* Steady State Value */
};
typedef boost::shared_ptr<NMSEPerformance> NMSEPerformance_sPtr;

struct StepPerformance
{
	GA_PIDChromosome<double> pidValues = { -1.0, -1.0, -1.0 };	/* PID Values that gave below performance metrics */
	Eigen::MatrixXd data;										/* Raw data that gave below performance metrics */

	double percentOvershoot_performance = -1.0;					/* Units: none, percentage */
	double steadyStateError_performance = -1.0;					/* Units: none, absolute */
	double delta_overshoot_performance = -1.0;					/* Units: none, absolute */
	double settlingTime_performance = -1.0;						/* Units: seconds */
	double riseTime_performance = -1.0;							/* Units: seconds */
	double steadyStateValue_performance = -1.0;					/* Units: user defined by problem */
			
			
	double settlingPcntRange = 0.0;								/* Percent range around final value that determines if a signal has settled */
	int settlingTime_Idx = 0;									/* Index in given data series */
	int riseTime_Idx[2] = { 0, 0 };								/* Start/Stop index in given data series */
};
typedef boost::shared_ptr<StepPerformance> StepPerformance_sPtr;

struct RampPerformance
{
	//currently just a placeholder
};
typedef boost::shared_ptr<RampPerformance> RampPerformance_sPtr;

struct QuadraticPerformance
{
	//currently just a placeholder
};
typedef boost::shared_ptr<QuadraticPerformance> QuadraticPerformance_sPtr;

struct CustomPerformance
{
	//currently just a placeholder
};
typedef boost::shared_ptr<CustomPerformance> CustomPerformance_sPtr;


struct PID_PerformanceData
{
	PerformanceType performanceDataType;

	StepPerformance_sPtr stepPerformanceData;	
	RampPerformance_sPtr rampPerformanceData;		
	QuadraticPerformance_sPtr quadPerformanceData;	
	CustomPerformance_sPtr customPerformanceData;	
	MSEPerformance_sPtr meanSquaredError;
	NMSEPerformance_sPtr normMeanSquaredError;
};

struct PID_FitnessScores
{
	/* Overall fitness score */
	double fitness_total = 0.0;				

	/* Sub fitness score */
	double fitness_POS = 0.0;
	double fitness_TS = 0.0;
	double fitness_TR = 0.0;
	double fitness_SSER = 0.0;
};


/*-----------------------------------------------
* FCSOptimizer Data Types
*-----------------------------------------------*/

/**
* \brief Describes how the optimizer will behave when limiting PID values to a user range
*/
enum FCSOptimizer_TunerLimiterBehavior
{
	FCS_LIMITER_FORCE_TO_BOUNDARY,			/* Upon exceeding user limits, it caps at the max/min values */
	FCS_LIMITER_REGENERATE_CHROMOSOME,		/* Upon exceeding user limits, a new random chromosome will be generated */
	FCS_LIMITER_DEFAULT = FCS_LIMITER_FORCE_TO_BOUNDARY
};

/**
* \brief Specifies basic solver convergence constraints
*/
struct FCSOptimizer_BasicConstraints
{
	struct RunTime
	{
		int minutes = 3;
		int seconds = 0;
	};

	RunTime timeConstraints;			/* Time based limitation given in minutes and seconds */

	double overallPerformance = 0.95;	/* Exits the solver if a solution is found that meets "x" percent of the specified goals.
										Essentially this is setting how good is "good enough" of a solution. */

	int maxNumberOfThreads = 16;		/* Forcefully limit how many threads can be spawned for a given optimizer. It's useful to
										specify this in multiples of 2. Currently the maximum value is hard limited to 16 threads. */

	int systemPrecision = 3;			/* Limit how many decimal places are used in calculations. This can speed up result generation. */

										//Add more as needed 
};

/**
* \brief Specifies advanced solver convergence constraints
*/
struct FCSOptimizer_AdvConstraints
{
	uint32_t populationSize = 20;			/* GA Population Size */
	uint32_t generationLimit = 100;			/* GA Generation Limit before exit */

	int iterations_before_refresh = 10;		/* Max number of iterations without progress using current GA Methods */
	int refresh_parameter_number = 2;		/* Specifies how many methods to replace when refreshing */
	int mutation_severity = 8;				/* Specifies how strong of a mutation can occur. Min: 1; Max: 16, increasing severity */
	double mutation_threshold = 0.5;		/* Specifies probability threshold for mutation. Range: 0-1, lower #s == higher probability */

											//Add more as needed
	FCSOptimizer_TunerLimiterBehavior limitingBehavior = FCS_LIMITER_DEFAULT;
};

/**
* \brief Descriptive information about a population member
*/
struct FCSOptimizer_PopulationMember
{
	GA_METHOD_ModelEvaluation modelType;		/* Specifies what kind of simulation model was used to generate member data */

	GA_PIDChromosome<double> dChrom;			
	GA_PIDChromosome<uint16_t> u16Chrom;		

	PID_FitnessScores fitnessScores;			/* Fitness score for the specified chromosome */
	PID_PerformanceData evaluationPerformance;	/* Evaluation performance for the specified chromosome */
};

/**
* \brief Conversion constants between uint16_t and double
* Holds pre-calculated constants used to map doubles into uint16_t
* types and back. This is useful for when chromosomes are represented
* as some form of bits instead of real numbers
*/
struct FCSOptimizer_MappingCoeff
{
	double bytes_precision;
	double x_offset;
	double x_lo;
	double x_hi;
	double x_dPow;
	double x_bPow;
	double x_sF;
	double x_sR;
};



#endif /* TYPES_H_ */