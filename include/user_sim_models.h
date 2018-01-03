#pragma once
#ifndef USER_SIM_MODELS_H_
#define USER_SIM_MODELS_H_
/* C/C++ Includes */
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

/* Eigen Includes */
#include <eigen/Eigen>
#include <eigen/StdVector>

/* Local Includes */
#include "model.h"


class Example_SS_PID : public SS_NLTIV_ModelBase
{
public:
	void assignData(int argc, double* argv) override;

	Eigen::MatrixXd getA() override;
	Eigen::MatrixXd getB() override;
	Eigen::MatrixXd getC() override;
	Eigen::MatrixXd getD() override;
	Eigen::MatrixXd getX0() override;

	int getNumInputs();
	int getNumOutputs();
	int getNumStates();

	Example_SS_PID();
	~Example_SS_PID();
private:
	int inputs, outputs, states;
	Eigen::MatrixXd A, B, C, D, U, X0;
};

#endif