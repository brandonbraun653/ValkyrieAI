#pragma once
#ifndef MODEL_SIMULATION_H_
#define MODEL_SIMULATION_H_

/* C/C++ Includes */
#include <stdlib.h>

/* Boost Includes */
#include <boost/shared_ptr.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen.hpp>

/* Local Includes */
#include "config.h"
#include "types.h"

class StateSpaceSimulator
{
public:

	Eigen::MatrixXd stepResponse(const double start, const double stop, const double dt, 
		SS_ModelBase_sPtr model, PID_Values pid);

	Eigen::MatrixXd rampResponse(const double start, const double stop, const double dt,
		SS_NLTIVModel model, PID_Values pid);

	Eigen::MatrixXd quadResponse(const double start, const double stop, const double dt,
		SS_NLTIVModel model, PID_Values pid);

	Eigen::MatrixXd customResponse(const double start, const double stop, const double dt,
		SS_NLTIVModel model, PID_Values pid, Eigen::MatrixXd customInput);

	StateSpaceSimulator();
	~StateSpaceSimulator();
private:

};
typedef boost::shared_ptr<StateSpaceSimulator> StateSpaceSimulator_sPtr;

#endif