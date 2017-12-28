#include "user_sim_models.h"

Example_SS_PID::Example_SS_PID()
{
	/*  For R inputs and M outputs and N states,
		A: NxN
		B: NxR
		C: MxN
		D: MxR
		Y: Mx1
	 */
	inputs = 1;
	outputs = 1;
	states = 2;

	A.resize(states, states);
	B.resize(states, inputs);
	C.resize(outputs, states);
	D.resize(outputs, inputs);
	X0.resize(states, 1);

	/* Hard-code the SS Matrices */
	A << -8.202, -2.029, -0.149, - 3.25;
	B << 1.14, -1.23;
	C << 1.0, 0.0;
	D << 0.0;

	/* Define the initial conditions */
	X0 << 0.0, 0.0;
}

Example_SS_PID::~Example_SS_PID()
{

}

void Example_SS_PID::assignData(int argc, double* argv)
{
}

Eigen::MatrixXd Example_SS_PID::getA()
{
	return A;
}

Eigen::MatrixXd Example_SS_PID::getB()
{
	return B;
}

Eigen::MatrixXd Example_SS_PID::getC()
{
	return C;
}

Eigen::MatrixXd Example_SS_PID::getD()
{
	return D;
}

Eigen::MatrixXd Example_SS_PID::getX0()
{
	return X0;
}

int Example_SS_PID::getNumInputs()
{
	return inputs;
}

int Example_SS_PID::getNumOutputs()
{
	return outputs;
}

int Example_SS_PID::getNumStates()
{
	return states;
}