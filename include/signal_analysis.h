#pragma once
#ifndef SIGNAL_ANALYSIS_H_
#define SIGNAL_ANALYSIS_H_

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
#include <boost/container/vector.hpp>
#include <boost/timer/timer.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen.hpp>

/* Local Includes */
#include "ga_config.h"
#include "data.h"
#include "debugger.h"

/*-----------------------------------------------
* Structs/Enums
*-----------------------------------------------*/
struct InflectionPoints
{
	boost::container::vector<double> values;    /* Data Value */
	boost::container::vector<double> time;		/* Time Stamp */
	boost::container::vector<int> index;		/* Step Index */

	boost::container::vector<double> differences;
};

enum SystemDamping
{
	NOT_CONVERGED = -1,			/* Distance between inflection points too large or system hasn't converged */
	SUFFICIENTLY_DAMPED = 1,	/* Enough inflection points to run calculations */
	UNDER_DAMPED,				/* Fewer inflection points or framing may cause analysis issues */
	OVER_DAMPED					/* No inflection points */
};

class StepResponseAnalyzer
{
public:

	StepPerformance analyze(Eigen::MatrixXd data);

	StepResponseAnalyzer();
	~StepResponseAnalyzer();

private:

	bool searchDirection;

	/* Gathered run-time data */
	Eigen::MatrixXd sim_data;
	SystemDamping settling_state;
	InflectionPoints extrema;
	StepPerformance performance;

	bool prefilter(double abs_threshold);
	void solveSettlingValue();
	void solveOvershoot();
	void solveSettlingTime();
	void solveRiseTime();
	void solveSteadyStateError();
};

#endif