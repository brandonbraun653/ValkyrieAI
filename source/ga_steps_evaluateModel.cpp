#include "ga_steps_evaluateModel.h"


//////////////////////////////////////////////////////////////////
/* CLASS:  StateSpaceModel */
//////////////////////////////////////////////////////////////////
/*-----------------------------------------------
* Constructors/Destructor
*-----------------------------------------------*/
StateSpaceEvaluator::StateSpaceEvaluator()
{
	simulator = boost::make_shared<StateSpaceSimulator>();
}

StateSpaceEvaluator::~StateSpaceEvaluator()
{
}

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

		output.errorCode = 0; //currently just default to zero

		switch (input.simulationType)
		{
		case STEP:
			output.data = simulator->stepResponse(input.startTime, input.endTime, input.dt, input.model, input.pid);
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

/*-----------------------------------------------
* Private Functions
*-----------------------------------------------*/




//////////////////////////////////////////////////////////////////
/* CLASS:  NeuralNetworkModel */
//////////////////////////////////////////////////////////////////
/*-----------------------------------------------
* Constructors/Destructor
*-----------------------------------------------*/
NeuralNetworkEvaluator::NeuralNetworkEvaluator()
{
}

NeuralNetworkEvaluator::~NeuralNetworkEvaluator()
{
}

/*-----------------------------------------------
* Public Functions
*-----------------------------------------------*/
void NeuralNetworkEvaluator::evaluate(const NeuralNetworkModelInput input, NeuralNetworkModelOutput& output)
{

}

/*-----------------------------------------------
* Private Functions
*-----------------------------------------------*/
