#include "ga_steps_evaluateModel.h"

//////////////////////////////////////////////////////////////////
/* CLASS:  StateSpaceModel */
//////////////////////////////////////////////////////////////////
/*-----------------------------------------------
* Public Functions
*-----------------------------------------------*/
void StateSpaceEvaluator::evaluate(const StateSpaceModelInput input, StateSpaceModelOutput& output)
{
	#if (SS_TRACE_ENABLE == 1)
		#if (SS_TRACE_EXECUTION_TIME == 1)
		auto start = boost::chrono::high_resolution_clock::now();
		#endif
	#endif

		output.stepPerformance = boost::make_shared<StepPerformance>();


		/*-----------------------------
		* Simulate the state space model using the given parameters and then perform 
		* analysis on the raw data to calculate various metrics. 
		*----------------------------*/
		Eigen::MatrixXd raw_simulation_data;
		switch (input.simulationType)
		{
		case STEP:
			raw_simulation_data = simulator.stepResponse(input.startTime, input.endTime, input.dt, input.model, input.pid);
			output.stepPerformance = stepAnalyzer.analyze(raw_simulation_data);

			//Need to set these AFTER the simulation completes otherwise the output of analyze() resets them
			output.stepPerformance->data = raw_simulation_data;
			output.stepPerformance->pidValues = input.pid; 
			
			break;

		case RAMP:

			break;

		case QUADRATIC:

			break;

		case CUSTOM:

			break;

		default: 
			output.errorCode = -1;
			break;
		};

	#if (SS_TRACE_ENABLE == 1)
		#if (SS_TRACE_EXECUTION_TIME == 1)
		auto stop = boost::chrono::high_resolution_clock::now();
		output.executionTime = boost::chrono::duration_cast<boost::chrono::microseconds>(stop - start);
		#endif
	#endif
}



//////////////////////////////////////////////////////////////////
/* CLASS:  NeuralNetworkModel */
//////////////////////////////////////////////////////////////////
/*-----------------------------------------------
* Public Functions
*-----------------------------------------------*/
void NeuralNetworkEvaluator::evaluate(const NeuralNetworkModelInput input, NeuralNetworkModelOutput& output)
{
	/* For now assume that we are only using the TCP model type */
	//TODO: Version3, change this to an actual input parameter
	NN_TCPModel_sPtr nn_model = boost::dynamic_pointer_cast<NN_TCPModel, NN_ModelBase>(input.model);
	

	/* Fill in the input data */
	SimCommand cmd;

	switch (input.simulationType)
	{
	case IMPULSE:
		cmd.simType = "impulse";
	case STEP:
		cmd.simType = "step";
		break;
	case RAMP:
		cmd.simType = "ramp";
		break;
	case QUADRATIC:
		cmd.simType = "quadratic";
		break;
	case CUSTOM:
		cmd.simType = "custom";
		break;
	default:
		throw std::runtime_error("No known simulation input type. Gonna crash now.");
	}

	cmd.axis = "pitch"; //TODO: Need to create as a parameter in init_struct perhaps 
	cmd.numTimeSteps = (int)floor((input.endTime - input.startTime) / input.dt);
	cmd.dt = (float)input.dt;
	cmd.start_time = (float)input.startTime;
	cmd.end_time = (float)input.endTime;
	cmd.stepMagnitude = (float)input.step_magnitude;
	cmd.angle_kp = (float)input.pid.Kp;
	cmd.angle_ki = (float)input.pid.Ki;
	cmd.angle_kd = (float)input.pid.Kd;

	std::string cmd_str = parseStruct(cmd);

	/* Send it off to the python code via TCP */
	nn_model->send_data(cmd_str);

	/* Block until we get some results back */
	SimResults results = parseResults(nn_model->recv_data());

	/* Fill in the output struct with the gathered data */
	StepPerformance_sPtr stepData = boost::make_shared<StepPerformance>();

	stepData->pidValues = input.pid;
	stepData->percentOvershoot_performance	= (isnan(results.overshoot) ? -1.0f : results.overshoot);
	stepData->riseTime_performance			= (isnan(results.riseTime) ? -1.0f : results.riseTime);
	stepData->settlingTime_performance		= (isnan(results.settlingTime) ? -1.0f : results.settlingTime);
	
	//Return the full data output and look at the last value?
	//stepData->steadyStateError_performance = results.

	//TODO: I need the full data!!!! both time and magnitude!! yikes that is a lot.

	output.stepPerformance = stepData;
}
