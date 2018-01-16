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
#include "config.h"
#include "debugger.h"
#include "types.h"
#include "model_simulation.h"
#include "signal_analysis.h"

namespace odeint = boost::numeric::odeint;
void model_error_exit(std::string error_msg, int line = __LINE__, std::string file = __FILE__);
bool dim_assert(size_t row_act, size_t row_exp, size_t col_act, size_t col_exp);





// 
// class StateSpaceModel
// {
// public:
// 	
// 	
// 
// 	StepPerformance stepResponseSingleThreaded(int member_num, SS_NLTIVModel &model, double Kp, double Ki, double Kd);
// 
// 
// 	void stepResponseMultiThreaded(int member_num, SS_NLTIVModel model, StepPerformance_Vec& output_data, boost::mutex& output_data_mutex,
// 		double Kp, double Ki, double Kd);
// 
// 	StepPerformance step_performance;
// 
// 	//Simulation specific data here
// 	int total_time_steps;
// 
// 	StateSpaceModel();
// 	~StateSpaceModel();
// 
// private:
// 	boost::mutex csv_mutex;
// 	boost::mutex cout_mutex;
// 	boost::mutex test_mutex;
// 
// 	Eigen::MatrixXd raw_data;
// 	StateSpaceSimulator step_simulator;
// 	StepResponseAnalyzer step_analyzer;
// };


#endif /* !MODEL_H_ */