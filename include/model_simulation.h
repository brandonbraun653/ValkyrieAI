#pragma once
#ifndef MODEL_SIMULATION_H_
#define MODEL_SIMULATION_H_

/* C/C++ Includes */
#include <stdlib.h>

/* Boost Includes */
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen.hpp>

/* Local Includes */
#include "ga_config.h"
#include "data.h"

class StepResponseSimulator
{
public:
	Eigen::MatrixXd simulate(SS_NLTIV_Dynamics model, double Kp, double Ki, double Kd);

	StepResponseSimulator();
	~StepResponseSimulator();
private:
};

#endif