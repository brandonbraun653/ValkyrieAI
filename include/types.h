#pragma once
#ifndef TYPES_H_
#define TYPES_H_

/* C/C++ Includes */
#include <stdlib.h>
#include <stdint.h>

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
	double steadyStateError_goal = 0.0;	/* Units: none, absolute */
	double settlingTime_goal = 0.0;		/* Units: seconds */
	double riseTime_goal = 0.0;			/* Units: seconds */
};

/** 
* \brief Specified tolerance (% +/- of goal) for each parameter
*/
struct PID_PerformanceTolerance
{
	double percentOvershoot_pcntTol = 0.0;
	double steadyStateError_pcntTol = 0.0;
	double settlingTime_pcntTol = 0.0;
	double riseTime_pcntTol = 0.0;
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
		double upper;
		double lower;
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

/**
* \brief Chromosome representation for a set of PID values
*
* A highly generic way of representing PID data so that many 
* different approaches can be used for manipulation
*/
template<typename ChromType>
struct GA_PIDChromosome
{
	ChromType Kp;
	ChromType Ki;
	ChromType Kd;
};



/*-----------------------------------------------
* FCSOptimizer Data Types
*-----------------------------------------------*/

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
	uint32_t populationSize = 50;			/* GA Population Size */
	uint32_t generationLimit = 50;			/* GA Generation Limit before exit */

	int iterations_before_refresh = 10;		/* Max number of iterations without progress using current GA Methods */
	int refresh_parameter_number = 2;		/* Specifies how many methods to replace when refreshing */
	int mutation_severity = 8;				/* Specifies how strong of a mutation can occur. Min: 1; Max: 16, increasing severity */
	double mutation_threshold = 0.5;		/* Specifies probability threshold for mutation. Range: 0-1, lower #s == higher probability */

	//Add more as needed 
};

/**
* \brief Descriptive information about a population member
*/
struct FCSOptimizer_PopulationMember
{
	GA_PIDChromosome<double> realPID;			/* Real valued representation of PID tuning parameters */

	GA_PIDChromosome<uint16_t> mappedPID;		/* Mapped representation (double -> uint16_t) of realPID over the min/max
												range specified in the PID_ControlSettings.tuningLimits struct */

	double fitnessScore;						/* A generic score that is independent of fitness metrics. More specific
												performance indices will need to be stored elsewhere */
	
	int age;									/* Records how many iterations the particular instance has been alive */
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


/*-----------------------------------------------
* Model Data Types
*-----------------------------------------------*/

class SS_NLTIV_Dynamics
{
public:
	double integrator_dt, integrator_time_start, integrator_time_end;
	Eigen::MatrixXd A, B, C, D, U, X0;
	int inputs, outputs, states;

	SS_NLTIV_Dynamics(const int Inputs, const int Outputs, const int States)
	{
		A.resize(States, States); A.setZero(States, States);
		B.resize(States, Inputs); B.setZero(States, Inputs);
		C.resize(Outputs, States); C.setZero(Outputs, States);
		D.resize(Outputs, Inputs); D.setZero(Outputs, Inputs);
		X0.resize(States, 1); X0.setZero(States, 1);
		U.resize(Inputs, 1); U.setZero(Inputs, 1);
		integrator_dt = 0.0;
		integrator_time_start = 0.0;
		integrator_time_end = 0.0;

		inputs = Inputs;
		outputs = Outputs;
		states = States;
	}

	/* System Function */
	void operator()(const Eigen::MatrixXd x, Eigen::MatrixXd &dxdt, double t)
	{
		dxdt = A*x + B*U;
	}

private:
};

class SS_NLTIV_Observer
{
public:
	uint32_t currentStep;
	Eigen::MatrixXd C, D, U, &Y;

	/* Constructor */
	SS_NLTIV_Observer(Eigen::MatrixXd ssC, Eigen::MatrixXd ssD, Eigen::MatrixXd ssU, Eigen::MatrixXd& ssY)
		: C(ssC), D(ssD), U(ssU), Y(ssY)
	{
		currentStep = 0;
		Y.resizeLike(U);
		Y.setZero(U.rows(), U.cols());
	}

	/* Observer Function */
	void operator()(const Eigen::MatrixXd &x, double t)
	{
		Y.col(currentStep) = C*x + D*U.col(currentStep);
		currentStep += 1;
	}

private:
};

class SS_NLTIV_ModelBase
{
public:
	virtual void assignData(int argc, double* argv) = 0;

	virtual Eigen::MatrixXd getA() = 0;
	virtual Eigen::MatrixXd getB() = 0;
	virtual Eigen::MatrixXd getC() = 0;
	virtual Eigen::MatrixXd getD() = 0;
	virtual Eigen::MatrixXd getX0() = 0;

	virtual int getNumInputs() = 0;
	virtual int getNumOutputs() = 0;
	virtual int getNumStates() = 0;

private:
};
typedef std::shared_ptr<SS_NLTIV_ModelBase> SS_NLTIV_ModelBase_sPtr;

/*-----------------------------------------------
* OUTPUT DATA
*-----------------------------------------------*/

/**
*	\brief Contains full performance analysis of a system's step response
*/
struct StepPerformance
{
	double Kp = -1.0;
	double Ki = -1.0;
	double Kd = -1.0;

	double percentOvershoot_performance = -1.0;	/* Units: none, percentage */
	double steadyStateError_performance = -1.0;	/* Units: none, absolute */
	double delta_overshoot_performance = -1.0;	/* Units: none, absolute */
	double settlingTime_performance = -1.0;		/* Units: seconds */
	double riseTime_performance = -1.0;			/* Units: seconds */
	double finalValue_performance = -1.0;		/* Units: user defined by problem */

	double settlingTime_window = 0.0;	/* Percent range around final value */
	int settlingTime_Idx = 0;			/* Index in given data series */
	int riseTime_Idx[2] = { 0, 0 };		/* Start/Stop index in given data series */

										/* Stores the data that resulted in performance criteria above */
	bool performance_simulation_data_is_valid = false;
	Eigen::MatrixXd performance_simulation_data;
};
typedef boost::shared_ptr<StepPerformance> StepPerformance_sPtr;
typedef boost::container::vector<StepPerformance> StepPerformance_Vec;

/**
*	\brief Fitness score for a given chromosome
*/
struct PID_FitnessScores
{
	//Overall fitness score
	double global_fitness = -1.0;

	//Sub fitness score
	double percentOvershoot_fitness = -1.0;
	double settlingTime_fitness = -1.0;
	double riseTime_fitness = -1.0;
	double steadyStateError_fitness = -1.0;

	StepPerformance stepPerformanceData;
};
typedef boost::shared_ptr<PID_FitnessScores> PIDFitnessScores_sPtr;
typedef boost::container::vector<PID_FitnessScores> PIDFitness_Vec;



#endif /* TYPES_H_ */