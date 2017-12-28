#pragma once
#ifndef DATA_H_
#define DATA_H_

/* C/C++ Includes */
#include <stdlib.h>
#include <stdint.h>
#include <memory>

/* Boost Includes */
#include <boost/chrono.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/container/vector.hpp>

/* Eigen Includes */
#include <eigen/Eigen>

/*-----------------------------------------------
* CONTROL FLOW
*-----------------------------------------------*/
/* Specify which device performs calculations */
#define USE_GPU 0x677075	/*<"gpu" in hex>*/
#define USE_CPU	0x637075	/*<"cpu" in hex>*/

/*-----------------------------------------------
* GENERIC DATA TYPES
*-----------------------------------------------*/
typedef boost::container::vector<int> iVec;
typedef boost::container::vector<double> dVec;
typedef boost::container::vector<uint16_t> u16Vec;

/*-----------------------------------------------
* INPUT DATA
*-----------------------------------------------*/
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

/* High level descriptors for drone tuning on each axis
 * without having to specify "engineering" performance metrics.
 */
struct Pilot_HighLevelResponse_TypeDef
{
	ControlResponseJargon roll;
	ControlResponseJargon pitch;
	ControlResponseJargon yaw;
};

/* Engineering performance characteristics for
 * one PID control loop. These are only the tuner goals.
 */
struct PID_PerformanceGoals
{
	double percentOvershoot_goal = 0.0;	/* Units: none, absolute percentage -> 0.02 == 2% */
	double steadyStateError_goal = 0.0;	/* Units: none, absolute */
	double settlingTime_goal = 0.0;		/* Units: seconds */
	double riseTime_goal = 0.0;			/* Units: seconds */
};

/* Specified tolerance (% +/- of goal) for each parameter
 */
struct PID_PerformanceTolerance
{
	double percentOvershoot_pcntTol = 0.0;
	double steadyStateError_pcntTol = 0.0;
	double settlingTime_pcntTol = 0.0;
	double riseTime_pcntTol = 0.0;
};

/* Specifies how important the parameter is to the user.
 * Each parameter can have a value in the range 0.0 - 1.0, with
 * 0.0 as unimportant and 1.0 very important.
 */
struct PID_PerformanceWeights
{
	double percentOvershoot_weight = 0.0;
	double steadyStateError_weight = 0.0;
	double settlingTime_weight = 0.0;
	double riseTime_weight = 0.0;
};

/* Specifies an initial guess for a starting point in the
 * tuner algorithm. This is helpful if a known configuration
 * is close to the desired response.
 */
struct PID_Value
{
	double Kp = 0.0;
	double Ki = 0.0;
	double Kd = 0.0;
};

/* Specifies the lower/upper limit the tuner is allowed to
 * apply in the algorithm.
 */
struct PID_TuningLimits
{
	double Kp_limits_upper = 100.0;
	double Kp_limits_lower = 0.001;
	double Ki_limits_upper = 100.0;
	double Ki_limits_lower = 0.001;
	double Kd_limits_upper = 100.0;
	double Kd_limits_lower = 0.001;
};

/* A conglomeration of various structs that fully characterize
 * performance metrics for a single PID control loop.
 */
struct PID_ControlGoals
{
	PID_PerformanceGoals		performance_goals;
	PID_PerformanceTolerance	performance_tolerance;
	PID_PerformanceWeights		performance_weights;
	PID_TuningLimits			pid_limits;
	PID_Value					pid_initial_values;
};
typedef std::shared_ptr<PID_ControlGoals> PID_ControlGoals_sPtr;

/* Specifies what kind of convergence constraints are placed
 * upon the GA solver. Not all of these may be used.
 */
struct GA_ConverganceCriteria
{
	uint32_t populationSize = 0;
	uint32_t generationLimit = 0;
	int mutation_severity = 13;
	double mutation_threshold = 0.0;
};
typedef std::shared_ptr<GA_ConverganceCriteria> GA_ConverganceCriteria_sPtr;

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

struct hPID_DataVector
{
	dVec Kp;
	dVec Ki;
	dVec Kd;
};

struct hPID_Chromosomes
{
	u16Vec Kp;
	u16Vec Ki;
	u16Vec Kd;
};

struct hPID_SimVector
{
	// Basic Simulation Parameters
	size_t population_size;
	unsigned int simulation_steps;
	double simulation_dt_sec;
};

struct GAMOP_hPID_Data
{
	hPID_DataVector pid_data;
	hPID_SimVector sim_data;
};

typedef struct
{
	double bytes_precision;
	double x_offset;
	double x_lo;
	double x_hi;
	double x_dPow;
	double x_bPow;
	double x_sF;
	double x_sR;
} mapCoeff_t;

/* Non-Linear Time Invariant State Space Model:*/
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
struct PID_Fitness
{
	//Overall fitness score
	double global_fitness = -1.0;

	//Sub fitness score
	double percentOvershoot_fitness = -1.0;
	double settlingTime_fitness = -1.0;
	double riseTime_fitness = -1.0;
	double steadyStateError_fitness = -1.0;

	StepPerformance fitness_performance;
};
typedef boost::shared_ptr<PID_Fitness> PIDFitness_sPtr;
typedef boost::container::vector<PID_Fitness> PIDFitness_Vec;

struct PIDElitist
{
	/* Some data to keep track of things */
	int worst_performer_index;
	double worst_performer_value;

	PIDFitness_Vec EliteFitSolutions;
};

struct GA_EngineStatistics
{
	//add things here like total run time, performance, etc.
	boost::chrono::duration<double> totalRunTime;

	PID_Fitness topPerformer;
};
typedef boost::container::vector<GA_EngineStatistics> GAEngineStatistics_Vec;
#endif /* !DATA_H_ */