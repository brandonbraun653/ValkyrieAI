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

/* Local Includes */
#include "ga_config.h"
#include "debugger.h"
#include "data.h"
#include "model_simulation.h"
#include "signal_analysis.h"

namespace odeint = boost::numeric::odeint;
void model_error_exit(std::string error_msg, int line = __LINE__, std::string file = __FILE__);
bool dim_assert(size_t row_act, size_t row_exp, size_t col_act, size_t col_exp);


/*-----------------------------------------------
* State Space Model
* Description:
*
*-----------------------------------------------*/
class StateSpaceModel
{
public:
	void init();
	void destroy();
	
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

	StateSpaceModel();
	~StateSpaceModel();

private:
	boost::mutex csv_mutex;
	boost::mutex cout_mutex;
	boost::mutex test_mutex;

	Eigen::MatrixXd raw_data;
	StepResponseSimulator step_simulator;
	StepResponseAnalyzer step_analyzer;
};
typedef boost::shared_ptr<StateSpaceModel> SSModel_sPtr;


/*-----------------------------------------------
* Neural Network Model
* Description:
*
*-----------------------------------------------*/
class NeuralNetworkModel
{
public:

	//Will need to provide stubs to some external TensorFlow code somehow.

private:
};
typedef boost::shared_ptr<NeuralNetworkModel> NNModel_sPtr;

#endif /* !MODEL_H_ */