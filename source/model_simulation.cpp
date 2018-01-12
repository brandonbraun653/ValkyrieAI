#include "model_simulation.h"

namespace odeint = boost::numeric::odeint;

StateSpaceSimulator::StateSpaceSimulator()
{
}

StateSpaceSimulator::~StateSpaceSimulator()
{
}

Eigen::MatrixXd StateSpaceSimulator::stepResponse(const double start, const double stop, const double dt, 
	SS_ModelBase_sPtr model, PID_Values pid)
{
	/* Due to how Odeint works internally, it's not feasible to pass in a full state space model object 
	and run the algorithm directly on that. (Odeint REALLY dislikes pointers) Instead, a local copy of the 
	model will be made. For now, assume that only NLTIV models are handled. */
	SS_NLTIVModel localModel(model);

	/* Copy the relevant data from the model for faster access/cleaner code in main loop */
	Eigen::MatrixXd x = localModel.getX0();
	Eigen::MatrixXd A = localModel.getA();
	Eigen::MatrixXd B = localModel.getB();
	Eigen::MatrixXd C = localModel.getC();
	Eigen::MatrixXd D = localModel.getD();

	
	/*-----------------------------------------------
	* Create the stepping engine that does all the simulation
	* (Optionally, later include methods to change how this is done dynamically)
	*-----------------------------------------------*/
	odeint::runge_kutta4<Eigen::MatrixXd, double, Eigen::MatrixXd, double, odeint::vector_space_algebra> stepper;


	/*-----------------------------------------------
	* Calculate any boundary conditions 
	*-----------------------------------------------*/
	int numSteps = (int)ceil((stop - start) / dt);	/* Total iterations the Odeint solver will run. Rounded up for non integral values */
	int outputs = localModel.getNumOutputs();		/* Output vector column dimension, (Nx1) */
	int inputs = localModel.getNumInputs();			/* Input vector column dimension, (Nx1) */
	int states = localModel.getNumStates();			/* State vector column dimension, (Nx1) */


	/*-----------------------------------------------
	* Initialize all matrices to correct dimensions and values 
	*-----------------------------------------------*/
	Eigen::MatrixXd reference_input(inputs, 1);		/* User desired set point for the PID controller to reach. Ultimately this translates
													into the specific evaluation type, which is in this case a step function. */
	reference_input.setConstant(1.0);				/* Simple step function input beginning at t0. */

	Eigen::MatrixXd U(inputs, 1);					/* Actual plant input after modification by PID controller */
	U.setZero();

	Eigen::MatrixXd error(inputs, 1);				/* PID runtime error value */
	error.setZero();

	Eigen::MatrixXd integral(inputs, 1);			/* PID runtime integral value */
	integral.setZero();

	Eigen::MatrixXd derivative(inputs, 1);			/* PID runtime derivative value */
	derivative.setZero();

	Eigen::MatrixXd previousError(inputs, 1);		/* PID runtime previousError value (simple memory for derivative calculation)*/
	previousError.setZero();

	Eigen::MatrixXd dxdt = x;						/* Initialize dxdt assuming starting from resting conditions */
	
	Eigen::MatrixXd y;								/* Intermediate output variable for model (raw calculation output only) */
	
	Eigen::MatrixXd yOut(outputs + 1, numSteps);	/* Actual output variable reported to user. To be concatenated with timestamp. */

	

	/*-----------------------------------------------
	* Run the whole simulation using PID feedback control
	*-----------------------------------------------*/
	double currentTime = start;
	for (int i = 0; i < numSteps; i++)
	{
		/* Take one step with the given engine */
		stepper.do_step(localModel, x, currentTime, dxdt, dt);

		/* Calculate the output */
		y = C*x + D*U;

		/* Log the output to the user buffer */
		yOut.col(i) << y, currentTime;

		/* Calculate the runtime PID terms */
		error = reference_input - y;
		integral += error*dt;
		derivative = (error - previousError) / dt;

		/* Calculate the control input for next step */
		U = pid.Kp*error + pid.Ki*integral + pid.Kd*derivative;

		/* Update variables for next round */
		x = dxdt;
		currentTime += dt;
		previousError = error;
		localModel.U = U;			/* The model operator() needs to be made aware of the new input value. This form
									is necessary due to how Odeint expects ODE models to be structured. */
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