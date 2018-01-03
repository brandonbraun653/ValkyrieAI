#include "model_simulation.h"

namespace odeint = boost::numeric::odeint;

StepResponseSimulator::StepResponseSimulator()
{
}

StepResponseSimulator::~StepResponseSimulator()
{
}

Eigen::MatrixXd StepResponseSimulator::simulate(SS_NLTIV_Dynamics model, double Kp, double Ki, double Kd)
{
	odeint::runge_kutta4<Eigen::MatrixXd, double, Eigen::MatrixXd, double, odeint::vector_space_algebra> stepper;

	int numSteps = floor((model.integrator_time_end - model.integrator_time_start) / model.integrator_dt);
	double t = model.integrator_time_start;
	double dt = model.integrator_dt;

	Eigen::MatrixXd x = model.X0, dxdt = x, y, yOut(model.outputs + 1, numSteps);

	/* Setup the input so that it matches the output Y dimensions */
	Eigen::MatrixXd
		reference_input(model.inputs, 1),
		error(model.inputs, 1),
		integral(model.inputs, 1),
		derivative(model.inputs, 1),
		previousError(model.inputs, 1);

	/* 1.0 because this is a step response simulator...*/
	reference_input.setConstant(1.0);

	error.setZero();
	integral.setZero();
	derivative.setZero();
	previousError.setZero();

	for (int i = 0; i < numSteps; i++)
	{
		stepper.do_step(model, x, t, dxdt, dt);

		// Calculate the output
		y = model.C*x + model.D*model.U;

		// Log the output
		yOut.col(i) << y, t;

		// Calculate the PID terms
		error = reference_input - y;
		integral += error*dt;
		derivative = (error - previousError) / dt;

		// Calculate the control input for next step
		model.U = Kp*error + Ki*integral + Kd*derivative;

		// Update variables
		x = dxdt;
		t += dt;
		previousError = error;
	}

	return yOut;
}