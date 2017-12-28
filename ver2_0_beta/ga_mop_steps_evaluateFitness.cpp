#include "ga_mop_steps_evaluateFitness.h"

bool filterData(double input, double abs_threshold)
{
	if ((abs(input) < abs_threshold) && (input != -1.0))
		return true;
	else
		return false;
}

//////////////////////////////////////////////////////////////////
/* CLASS:  WeightedSum */
//////////////////////////////////////////////////////////////////
/*-----------------------------------------------
* Constructors/Destructor
*-----------------------------------------------*/
WeightedSum::WeightedSum(GA_RunMode execution_type)
{
	executionType = execution_type;
}

WeightedSum::~WeightedSum()
{
}

/*-----------------------------------------------
* Public Functions
*-----------------------------------------------*/
void WeightedSum::calculateFitness(StepPerformance_Vec input_data, PID_ControlGoals_sPtr input_goals, PIDFitness_Vec* output_fitness)
{
	ws_data = input_data;			//Explicit copy
	ws_goals = input_goals;			//Explicit copy
	ws_fitness = output_fitness;	//Remember to pass this in via boost::ref()

	if (executionType == SINGLE_THREADED)
		calculate_cpu_single_threaded();

	if (executionType == MULTI_THREADED)
		calculate_cpu_multi_threaded();

	if (executionType == SINGLE_THREADED_WITH_CUDA)
		calculate_gpu_single_threaded();

	if (executionType == MULTI_THREADED_WITH_CUDA)
		calculate_gpu_multi_threaded();
}

/*-----------------------------------------------
* Private Functions
*-----------------------------------------------*/
PID_Fitness WeightedSum::calculateMemberFit(int memberNum, bool dataValid, double POS, double TS, double TR, double SSERR)
{
	PID_Fitness localFit;
	localFit.fitness_performance = ws_data.data()[memberNum];
	
	if (dataValid)
	{
		/* Make a copy of all the user goals */
		PID_ControlGoals_sPtr lws_goals = std::make_shared<PID_ControlGoals>(*ws_goals);

		/* Create some additional computational variables */
		localFit.global_fitness = 0.0;
		double perfect_fit_val = 1.0;
		double numValidParameters = 4.0;
		double error = 0.0;
		double abs_pct_error = 0.0;
		double abs_error = 0.0;

		/*-----------------------------------------------
		* Percent Overshoot
		*-----------------------------------------------*/
		if (POS != -1.0 && lws_goals->performance_goals.percentOvershoot_goal != -1.0)
		{
			if (POS == 0.0)
			{
				/* Give partial credit for achieving this. It means that the simulation is more
				   or less an overdamped system that may or may not have properly converged. However,
				   due to the likelihood of a good solution, some credit must be given. This value
				   is arbitrarily chosen.*/
				localFit.percentOvershoot_fitness = 0.25;
			}
			else
			{
				#ifdef DEBUGGING_ENABLED
				double maxVal = ws_data.data()[memberNum].performance_simulation_data.maxCoeff();
				#endif

				// I guess this can represent percent undershoot as well...
				error = abs(POS) - 100.0*lws_goals->performance_goals.percentOvershoot_goal;
				abs_pct_error = abs(error) - (100.0*lws_goals->performance_goals.percentOvershoot_goal);

				/* If we got negative error or are within tolerance, this is good and the solution
				is given a perfect score. If not, figure out how to derate accordingly. */
				if (error <= 0.0 || abs_pct_error < lws_goals->performance_tolerance.percentOvershoot_pcntTol)
					localFit.percentOvershoot_fitness = perfect_fit_val;

				/* Derate the solution because it's baaaad. */
				else
				{
					abs_error = abs_pct_error*lws_goals->performance_goals.percentOvershoot_goal;
					localFit.percentOvershoot_fitness = perfect_fit_val*exp(-1.0*abs_error);
				}
			}
			
			localFit.global_fitness += localFit.percentOvershoot_fitness;
		}

		/*-----------------------------------------------
		* Settling Time
		*-----------------------------------------------*/
		if (TS != -1.0 && lws_goals->performance_goals.settlingTime_goal != -1.0)
		{
			if (TS == 0.0)
			{
				/* Same reasoning as POS is used here. */
				localFit.settlingTime_fitness = 0.25;
			}
			else
			{
				error = TS - lws_goals->performance_goals.settlingTime_goal;
				abs_pct_error = abs(error) / lws_goals->performance_goals.settlingTime_goal;

				if (error <= 0.0 || abs_pct_error < lws_goals->performance_tolerance.settlingTime_pcntTol)
					localFit.settlingTime_fitness = perfect_fit_val;

				else
				{
					abs_error = abs_pct_error*lws_goals->performance_goals.settlingTime_goal;
					localFit.settlingTime_fitness = perfect_fit_val*exp(-1.0*abs_error);
				}
			}
			
			localFit.global_fitness += localFit.settlingTime_fitness;
		}

		/*-----------------------------------------------
		* Rise Time
		*-----------------------------------------------*/
		if (TR != -1.0 && lws_goals->performance_goals.riseTime_goal != -1.0)
		{
			error = TR - lws_goals->performance_goals.riseTime_goal;
			abs_pct_error = abs(error) / lws_goals->performance_goals.riseTime_goal;

			if (error <= 0.0 || abs_pct_error < lws_goals->performance_tolerance.riseTime_pcntTol)
				localFit.riseTime_fitness = perfect_fit_val;

			else
			{
				abs_error = abs_pct_error*lws_goals->performance_goals.riseTime_goal;
				localFit.riseTime_fitness = perfect_fit_val*exp(-1.0*abs_error);
			}
			localFit.global_fitness += localFit.riseTime_fitness;
		}

		/*-----------------------------------------------
		* Steady State Error
		*-----------------------------------------------*/
		if (SSERR != -1.0 && lws_goals->performance_goals.steadyStateError_goal != -1.0)
		{
			error = abs(SSERR) - lws_goals->performance_goals.steadyStateError_goal;
			abs_pct_error = abs(error) / lws_goals->performance_goals.steadyStateError_goal;

			if (error <= 0.0 || abs_pct_error < lws_goals->performance_tolerance.steadyStateError_pcntTol)
				localFit.steadyStateError_fitness = perfect_fit_val;

			else
			{
				abs_error = abs_pct_error*lws_goals->performance_goals.steadyStateError_goal;

				/* The value is shifted here to enforce a harsher penalty on a solution that 
				   does not quite hit the steady state error goal & tolerance. It was found 
				   that without shifting, the converged solution would routinely exhibit too
				   large of an error. */
				localFit.steadyStateError_fitness = perfect_fit_val*exp(-1.0*(abs_error+0.65));
			}
			localFit.global_fitness += localFit.steadyStateError_fitness;
		}

		/*-----------------------------------------------
		* Global Fit
		*-----------------------------------------------*/
		localFit.global_fitness /= numValidParameters;
	
		#ifdef DEBUGGING_ENABLED
		double posFit = localFit.percentOvershoot_fitness;
		double tSettleFit = localFit.settlingTime_fitness;
		double tRiseFit = localFit.riseTime_fitness;
		double ssErrorFit = localFit.steadyStateError_fitness;
		double globalFit = localFit.global_fitness;
		bool bp_flag = false;

		if (posFit == 1.0 || tSettleFit == 1.0 || tRiseFit == 1.0 ||\
			ssErrorFit == 1.0 || globalFit == 1.0)
		{
			bp_flag = true;
			//std::cout << "Got a perfect fit on a performance criteria" << std::endl;
		}
		#endif
	}
			
	if (executionType == MULTI_THREADED)
	{
		ws_fitness_mutex.lock();
		ws_fitness->data()[memberNum] = localFit;
		ws_fitness_mutex.unlock();
	}

	return localFit;
}

void WeightedSum::calculate_cpu_single_threaded()
{
	size_t popSize = ws_fitness->size();

	for (int member = 0; member < popSize; member++)
	{

		#ifdef DEBUGGING_ENABLED
		bool validData = false;
		validData = ws_data.data()[member].performance_simulation_data_is_valid;
		#endif

		ws_fitness->data()[member] = calculateMemberFit(
			member,
			ws_data.data()[member].performance_simulation_data_is_valid,
			ws_data.data()[member].percentOvershoot_performance,
			ws_data.data()[member].settlingTime_performance,
			ws_data.data()[member].riseTime_performance,
			ws_data.data()[member].steadyStateError_performance);
	}
		
}

void WeightedSum::calculate_cpu_multi_threaded()
{
	size_t popSize = ws_fitness->size();

	boost::thread_group tgroup;

	for (int member = 0; member < popSize; member++)
	{
		tgroup.create_thread(boost::bind(&WeightedSum::calculateMemberFit, this,
			member,
			ws_data.data()[member].performance_simulation_data_is_valid,
			ws_data.data()[member].percentOvershoot_performance,
			ws_data.data()[member].settlingTime_performance,
			ws_data.data()[member].riseTime_performance,
			ws_data.data()[member].steadyStateError_performance
			));
	}

	tgroup.join_all();
}

void WeightedSum::calculate_gpu_single_threaded()
{
	std::cout << "GPU Single Threaded Mode Currently Not Supported." << std::endl;
}

void WeightedSum::calculate_gpu_multi_threaded()
{
	std::cout << "GPU Multi Threaded Mode Currently Not Supported." << std::endl;
}

//////////////////////////////////////////////////////////////////
/* CLASS:  NonDominatedSort */
//////////////////////////////////////////////////////////////////
/*-----------------------------------------------
* Constructors/Destructor
*-----------------------------------------------*/

/*-----------------------------------------------
* Public Functions
*-----------------------------------------------*/

/*-----------------------------------------------
* Private Functions
*-----------------------------------------------*/