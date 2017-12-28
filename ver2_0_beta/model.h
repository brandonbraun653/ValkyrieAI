#pragma once
#ifndef MODEL_H_
#define MODEL_H_

/* C/C++ Includes */
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <memory>

/* Eigen Includes */
#include <eigen/Eigen>
#include <eigen/StdVector>

/* Boost Includes */
#include <boost/thread.hpp>
#include <boost/chrono.hpp>
#include <boost/container/vector.hpp>
#include <boost/timer/timer.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen.hpp>

/* Thrust Includes */
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

/* Local Includes */
#include "ga_config.h"
#include "debugger.h"
#include "data.h"
#include "model_simulation.h"
#include "signal_analysis.h"

namespace odeint = boost::numeric::odeint;
void model_error_exit(std::string error_msg, int line = __LINE__, std::string file = __FILE__);
bool dim_assert(size_t row_act, size_t row_exp, size_t col_act, size_t col_exp);


////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
/*-----------------------------------------------
* Basic Genetic Algorithm Model
* Description:
*	This base is intended to serve as a template for what
*	kinds of functions are to be made available in inheriting
*	classes. In addition, the polymorphic structure allows for
*	easy assignment of a diverse range of model types to a single
*	GA Engine instance.
*-----------------------------------------------*/
class GAModel
{
public:
	virtual void init() = 0;
	virtual void destroy() = 0;

	GAModel();
	~GAModel();

private:
};
typedef std::shared_ptr<GAModel> GAModel_sPtr;

/*-----------------------------------------------
* State Space Model
* Description:
*
*-----------------------------------------------*/
class SSModel : public GAModel
{
public:
	void init() override;
	void destroy() override;
	
	void assignSimulationTimeConstraints(double dt, double start_time, double end_time);
	

	StepPerformance stepResponseSingleThreaded(int member_num, SS_NLTIV_Dynamics &model, double Kp, double Ki, double Kd);


	void stepResponseMultiThreaded(int member_num, SS_NLTIV_Dynamics model, StepPerformance_Vec& output_data, boost::mutex& output_data_mutex,
		double Kp, double Ki, double Kd);

	StepPerformance step_performance;

	//Simulation specific data here
	double sim_dt;
	double sim_start_time;
	double sim_end_time;
	int total_time_steps;

	SSModel();
	~SSModel();

private:
	boost::mutex csv_mutex;
	boost::mutex cout_mutex;
	boost::mutex test_mutex;

	Eigen::MatrixXd raw_data;
	StepResponseSimulator step_simulator;
	StepResponseAnalyzer step_analyzer;
};
extern std::string SSModelName;
typedef std::shared_ptr<SSModel> SSModel_sPtr;

/*-----------------------------------------------
* Dynamic Model
* Description:
*	This model implements a state-space version of a linear or non-linear
*	control system.
*
* Notes: Use Matlab to generate the system matrices.
*-----------------------------------------------*/
class DynamicModel : public GAModel
{
public:
	void init() override;
	void destroy() override;
	void evaluate(int processor, GAMOP_hPID_Data *data);

	DynamicModel();
	~DynamicModel();

private:
};
extern std::string DYModelName;
typedef std::shared_ptr<DynamicModel> DynamicModel_sPtr;

/*-----------------------------------------------
* Live Model
* Description:
*
*-----------------------------------------------*/
class LiveModel : public GAModel
{
public:
	void init() override;
	void destroy() override;
	void evaluate(int processor, GAMOP_hPID_Data *data);

	LiveModel();
	~LiveModel();

private:
};
extern std::string LVModelName;
typedef std::shared_ptr<LiveModel> LiveModel_sPtr;

#endif /* !MODEL_H_ */