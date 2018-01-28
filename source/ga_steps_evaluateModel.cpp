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

}
