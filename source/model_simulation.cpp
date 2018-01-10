#include "model_simulation.h"

namespace odeint = boost::numeric::odeint;

StateSpaceSimulator::StateSpaceSimulator()
{
}

StateSpaceSimulator::~StateSpaceSimulator()
{
}
//TODO: CHANGE THIS TO PASS BY REFERENCE!!!!
Eigen::MatrixXd StateSpaceSimulator::stepResponse(const double start, const double stop, const double dt, 
	SS_NLTIVModel model, PID_Values pid)
{
	/* Copy the relevant data from the model */
	Eigen::MatrixXd x = model.getX0();
	Eigen::MatrixXd A = model.getA();
	Eigen::MatrixXd B = model.getB();
	Eigen::MatrixXd C = model.getC();
	Eigen::MatrixXd D = model.getD();

	int outputs = model.getNumOutputs();
	int inputs = model.getNumInputs();
	int states = model.getNumStates();


	/* Create the stepping engine that does all the simulation */
	odeint::runge_kutta4<Eigen::MatrixXd, double, Eigen::MatrixXd, double, odeint::vector_space_algebra> stepper;


	int numSteps = (int)ceil((stop - start) / dt);
	double currentTime = start;


	Eigen::MatrixXd dxdt = x, y, yOut(outputs + 1, numSteps);

	/* Setup the inputs and PID terms so they match output Y dimensions */
	Eigen::MatrixXd reference_input(inputs, 1); 
	Eigen::MatrixXd U(inputs, 1);				U.setZero();
	Eigen::MatrixXd error(inputs, 1);			error.setZero();
	Eigen::MatrixXd integral(inputs, 1);		integral.setZero();
	Eigen::MatrixXd derivative(inputs, 1);		derivative.setZero();
	Eigen::MatrixXd previousError(inputs, 1);	previousError.setZero();

	/* Create the simulation input */
	reference_input.setConstant(1.0);


	/* Run the whole simulation using PID feedback control */
	for (int i = 0; i < numSteps; i++)
	{
		// Take one step with the given engine
		stepper.do_step(model, x, currentTime, dxdt, dt);

		// Calculate the output
		y = C*x + D*U;

		// Log the output
		yOut.col(i) << y, currentTime;

		// Calculate the PID terms
		error = reference_input - y;
		integral += error*dt;
		derivative = (error - previousError) / dt;

		// Calculate the control input for next step
		U = pid.Kp*error + pid.Ki*integral + pid.Kd*derivative;

		// Update variables
		x = dxdt;
		currentTime += dt;
		previousError = error;
	}

	return yOut;
}


Eigen::MatrixXd StateSpaceSimulator::rampResponse(const double start, const double stop, const double dt,
	SS_NLTIVModel model, PID_Values pid)
{
	Eigen::MatrixXd temp;
	return temp;
}

Eigen::MatrixXd StateSpaceSimulator::quadResponse(const double start, const double stop, const double dt,
	SS_NLTIVModel model, PID_Values pid)
{
	Eigen::MatrixXd temp;
	return temp;
}

Eigen::MatrixXd StateSpaceSimulator::customResponse(const double start, const double stop, const double dt,
	SS_NLTIVModel model, PID_Values pid, Eigen::MatrixXd customInput)
{
	Eigen::MatrixXd temp;
	return temp;
}