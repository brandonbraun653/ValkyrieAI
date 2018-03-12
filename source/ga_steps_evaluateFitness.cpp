#include "ga_steps_evaluateFitness.h"

//////////////////////////////////////////////////////////////////
/* CLASS:  WeightedSum */
//////////////////////////////////////////////////////////////////
/*-----------------------------------------------
* Public Functions
*-----------------------------------------------*/
void WeightedSum::evaluateFitness(const GA_EvaluateFitnessDataInput input, GA_EvaluateFitnessDataOutput& output)
{
	//Note: If getting nan errors from this function, make sure -ffast-math isn't enabled
	//https://stackoverflow.com/questions/3596622/negative-nan-is-not-a-nan 

	PID_FitnessScores localFit;

	/* Copy in the performance metrics for easier reading */
	double POS = input.POS;
	double TS = input.TS;
	double TR = input.TR;
	double SSERR = input.SSER;


	/* Create some additional computational variables */
	localFit.fitness_total = 0.0;
	double perfect_fit_val = 1.0;
	double numValidParameters = 4.0;
	double error = 0.0;
	double abs_pct_error = 0.0;
	double abs_error = 0.0;

	/*-----------------------------------------------
	* Percent Overshoot
	*-----------------------------------------------*/
	if (POS != -1.0 && input.goals.percentOvershoot_goal != -1.0)
	{
		#ifdef _DEBUG
		if (input.goals.percentOvershoot_goal == 0.0)
			std::cout << "You forgot to define a percent overshoot goal." << std::endl;

		if (input.tolerance.percentOvershoot_pcntTol == 0.0)
			std::cout << "You forgot to define a percent overshoot goal tolerance." << std::endl;
		#endif

		if (POS == 0.0)
		{
			/* Give partial credit for achieving this. It means that the simulation is more
			or less an over damped system that may or may not have properly converged. This value
			is arbitrarily chosen.*/
			localFit.fitness_POS = 0.1;
		}
		else if (isnan(POS))
		{
			localFit.fitness_POS = 0.0;
		}
		else
		{
			#ifdef DEBUGGING_ENABLED
			//double maxVal = input.simulationData.maxCoeff();
			#endif

			// I guess this can represent percent undershoot as well...
			error = abs(POS) - 100.0*input.goals.percentOvershoot_goal;
			abs_pct_error = abs(error) - (100.0*input.goals.percentOvershoot_goal);

			/* If we got negative error or are within tolerance, this is good and the solution
			is given a perfect score. If not, figure out how to derate accordingly. */
			if (error <= 0.0 || abs_pct_error < input.tolerance.percentOvershoot_pcntTol)
				localFit.fitness_POS = perfect_fit_val;

			/* Derate the solution because it's baaaad. */
			else
			{
				abs_error = abs_pct_error*input.goals.percentOvershoot_goal;
				localFit.fitness_POS = perfect_fit_val*exp(-1.0*abs_error);
			}
		}

		localFit.fitness_total += localFit.fitness_POS;
	}

	/*-----------------------------------------------
	* Settling Time
	*-----------------------------------------------*/
	if (TS != -1.0 && input.goals.settlingTime_goal != -1.0)
	{
		#ifdef _DEBUG
		if (input.goals.settlingTime_goal == 0.0)
			std::cout << "You forgot to define a settling time goal." << std::endl;

		if (input.tolerance.settlingTime_pcntTol == 0.0)
			std::cout << "You forgot to define a settling time goal tolerance." << std::endl;
		#endif

		if (TS == 0.0)
		{
			/* Same reasoning as POS is used here. */
			localFit.fitness_TS = 0.1;
		}
		else if (isnan(TS))
		{
			localFit.fitness_TS = 0.0;
		}
		else
		{
			error = TS - input.goals.settlingTime_goal;
			abs_pct_error = abs(error) / input.goals.settlingTime_goal;

			if (error <= 0.0 || abs_pct_error < input.tolerance.settlingTime_pcntTol)
				localFit.fitness_TS = perfect_fit_val;

			else
			{
				abs_error = abs_pct_error*input.goals.settlingTime_goal;
				localFit.fitness_TS = perfect_fit_val*exp(-1.0*abs_error);
			}
		}

		localFit.fitness_total += localFit.fitness_TS;
	}

	/*-----------------------------------------------
	* Rise Time
	*-----------------------------------------------*/
	if (TR != -1.0 && input.goals.riseTime_goal != -1.0)
	{
		#ifdef _DEBUG
		if (input.goals.riseTime_goal == 0.0)
			std::cout << "You forgot to define a rise time goal." << std::endl;

		if (input.tolerance.riseTime_pcntTol == 0.0)
			std::cout << "You forgot to define a rise time goal tolerance." << std::endl;
		#endif

		if (isnan(TR))
		{
			localFit.fitness_TR = 0.0;
		}	
		else
		{
			error = TR - input.goals.riseTime_goal;
			abs_pct_error = abs(error) / input.goals.riseTime_goal;

			if (error <= 0.0 || abs_pct_error < input.tolerance.riseTime_pcntTol)
				localFit.fitness_TR = perfect_fit_val;

			else
			{
				abs_error = abs_pct_error * input.goals.riseTime_goal;
				localFit.fitness_TR = perfect_fit_val * exp(-1.0*abs_error);
			}
		}
		
		localFit.fitness_total += localFit.fitness_TR;
	}

	/*-----------------------------------------------
	* Steady State Error
	*-----------------------------------------------*/
	if (SSERR != -1.0 && input.goals.steadyStateError_goal != -1.0)
	{
		#ifdef _DEBUG
		if (input.goals.steadyStateError_goal == 0.0)
			std::cout << "You forgot to define a steady state error goal." << std::endl;

		if (input.tolerance.steadyStateError_pcntTol == 0.0)
			std::cout << "You forgot to define a steady stater error goal tolerance." << std::endl;
		#endif 

		if (isnan(SSERR))
		{
			localFit.fitness_SSER = 0.0;
		}
		else
		{
			error = abs(SSERR) - input.goals.steadyStateError_goal;
			abs_pct_error = abs(error) / input.goals.steadyStateError_goal;

			if (error <= 0.0 || abs_pct_error < input.tolerance.steadyStateError_pcntTol)
				localFit.fitness_SSER = perfect_fit_val;

			else
			{
				abs_error = abs_pct_error * input.goals.steadyStateError_goal;

				/* The value is shifted here to enforce a harsher penalty on a solution that
				does not quite hit the steady state error goal & tolerance. It was found
				that without shifting, the converged solution would routinely exhibit too
				large of an error. */
				localFit.fitness_SSER = perfect_fit_val * exp(-1.0*(abs_error + 0.65));
			}
		}

		localFit.fitness_total += localFit.fitness_SSER;
	}

	/*-----------------------------------------------
	* Global Fit
	*-----------------------------------------------*/
	localFit.fitness_total /= numValidParameters;

	output.fitness = localFit;
}



//////////////////////////////////////////////////////////////////
/* CLASS:  MeanSquareError */
//////////////////////////////////////////////////////////////////
/*-----------------------------------------------
* Public Functions
*-----------------------------------------------*/
void MeanSquareError::evaluateFitness(const GA_EvaluateFitnessDataInput input, GA_EvaluateFitnessDataOutput& output)
{

}

