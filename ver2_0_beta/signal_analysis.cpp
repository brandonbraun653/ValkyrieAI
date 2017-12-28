#include "signal_analysis.h"

//////////////////////////////////////////////////////////////////
/* CLASS: StepResponseAnalyzer */
//////////////////////////////////////////////////////////////////

/*-----------------------------------------------
* Constructors/Destructor
*-----------------------------------------------*/
StepResponseAnalyzer::StepResponseAnalyzer()
{
}

StepResponseAnalyzer::~StepResponseAnalyzer()
{
}

/*-----------------------------------------------
* Public Functions
*-----------------------------------------------*/
/* Function: Analyze
 *
 * Description:
 * The idea behind this function is that it's a one stop shop for completely
 * analyzing a single dimensional transient step response of some input data.
 * The function does not actually perform the stepping and simulating, but
 * only evaluates given transient data for performance characteristics.
 *
 * Input Type: Dynamic Matrix (double)
 *
 * Input Format: [2xN] Array with the first row corresponding to the output and
 * the second corresponding to the time steps
 */
StepPerformance StepResponseAnalyzer::analyze(Eigen::MatrixXd data)
{
	//Enforce formatting
	if (data.rows() != 2 || data.cols() == 0)
	{
		std::cout << "Incorrect matrix input for step response analyzer." << std::endl;
		return performance;
	}

	/*-----------------------------------------------
	* Initialization steps
	*-----------------------------------------------*/
	sim_data = data;
	performance.performance_simulation_data = data; //Probably should merge this

	/*-----------------------------------------------
	* Run through each step successively
	*-----------------------------------------------*/
	solveSettlingValue();
	solveOvershoot();
	solveSettlingTime();
	solveRiseTime();
	solveSteadyStateError();

	return performance;
}

/*-----------------------------------------------
* Private Functions
*-----------------------------------------------*/
bool StepResponseAnalyzer::prefilter(double abs_threshold)
{
	/* This function exists primarily to sort out any extremes in the 
	   returned simulation data. It was discovered through testing that
	   the other performance functions in this class were unable to discern 
	   between data that was well behaved and data that was not. */

	double maxVal = sim_data.maxCoeff();
	double minVal = sim_data.minCoeff();

	if ((maxVal > abs_threshold) || (minVal < -abs_threshold))
		return false;
	else
		return true;
}

void StepResponseAnalyzer::solveSettlingValue()
{
	/*-----------------------------------------------
	* Control variable setup
	*-----------------------------------------------*/
	settling_state = NOT_CONVERGED;
	double minSettlingVariance = 0.4;
	double maxSlopeCritDamp = 0.01;
	double maxSlopeDamp = 0.01;

	/*-----------------------------------------------
	* Figure out the initial search direction
	*-----------------------------------------------*/
	double avg = 0.0;
	for (int i = 0; i < 10; i++)
		avg += sim_data(0, i);

	avg /= 10.0;

	if (avg > sim_data(0, 0))
		searchDirection = true;		//Increasing search first
	else
		searchDirection = false;	//Decreasing search first

	/*-----------------------------------------------
	* Find all the inflection points
	*-----------------------------------------------*/
	double lastVal = sim_data(0, 0);
	for (int i = 0; i < sim_data.cols(); i++)
	{
		/* Searching for a maximum */
		if (searchDirection)
		{
			// Just switched from positive/zero slope to negative
			if (sim_data(0, i) < lastVal)
			{
				extrema.values.push_back(lastVal);
				extrema.time.push_back(sim_data(1, i - 1));
				extrema.index.push_back(i - 1);

				searchDirection = false;
			}

			// Still increasing, save for next iteration
			lastVal = sim_data(0, i);
		}

		/* Searching for a minimum */
		if (!searchDirection)
		{
			// Just switched from negative/zero slope to positive
			if (sim_data(0, i) > lastVal)
			{
				extrema.values.push_back(lastVal);
				extrema.time.push_back(sim_data(1, i - 1));
				extrema.index.push_back(i - 1);

				searchDirection = true;
			}

			// Still decreasing, save for next iteration
			lastVal = sim_data(0, i);
		}
	}

	/*-----------------------------------------------
	* Figure out what kind of system is present
	*-----------------------------------------------*/
	/* At least 1 inflection point is found */
	if (!extrema.values.empty())
	{
		/* If at least 1 pair of inflection points, calculate their difference and check for convergence */
		if (extrema.values.size() >= 2)
		{
			//Calculate distance pairs...
			for (int i = 0; i < extrema.values.size() - 1; i++)
				extrema.differences.push_back(abs(extrema.values[i] - extrema.values[i + 1]));

			//Find the minimum of all pairs
			double minDiffValue = *std::min_element(extrema.differences.begin(), extrema.differences.end());

			//Check for convergence
			if (minDiffValue < minSettlingVariance)
				settling_state = SUFFICIENTLY_DAMPED;
		}

		/* 1 or more inflection points found, but not converged yet. Still has a shot though! */
		if (settling_state == NOT_CONVERGED)
		{
			int simSize = sim_data.cols();			//Num time steps
			int infPtSize = extrema.values.size();	//Num inflection points

			//Check the slope between the last data point and the last inflection point.
			int nPts = simSize - extrema.index[infPtSize - 1];

			//testing
			double y1 = sim_data(0, simSize - 1);
			double y2 = extrema.values[infPtSize - 1];

			double ydiff = abs(y1 - y2);
			double slope = ydiff / (double)nPts;

			//If slope too large or there are too few points,
			//we reject the whole system as it isn't really settling out.
			if ((slope < maxSlopeDamp) && (nPts > 10) && (ydiff < 1.5*minSettlingVariance))
				settling_state = UNDER_DAMPED;
		}
	}

	/* No inflection points found, thus it may be critically damped or maybe exponential/ramp */
	else
	{
		//Figure out how many data points correspond to the last 15% of sim
		int lastIdx = sim_data.cols() - 1;
		int nPts = (int)floor(0.15*lastIdx);

		//Calculate the slope over the last 15% of sim
		double slope = (sim_data(0, lastIdx) - sim_data(0, lastIdx - nPts)) / nPts;

		//If the slope is small enough, we assume it's an over/critically damped system
		if (slope < maxSlopeCritDamp)
			settling_state = OVER_DAMPED;
	}

	/*-----------------------------------------------
	* Calculate the final value
	*-----------------------------------------------*/
	if (settling_state != NOT_CONVERGED)
	{
		if (settling_state == SUFFICIENTLY_DAMPED)
		{
			int eLastVal = extrema.values.size() - 1;
			int eDiffLastVal = extrema.differences.size() - 1;

			if (searchDirection)
				performance.finalValue_performance = extrema.values[eLastVal] + extrema.differences[eDiffLastVal] / 2.0;
			else
				performance.finalValue_performance = extrema.values[eLastVal] - extrema.differences[eDiffLastVal] / 2.0;
		}

		if (settling_state == UNDER_DAMPED)
		{
			int eLastVal = extrema.values.size() - 1;
			int simLastVal = sim_data.cols() - 1;

			double offset = (sim_data(0, simLastVal) - extrema.values[eLastVal]) / 2.0;
			if (searchDirection)
				performance.finalValue_performance = extrema.values[eLastVal] + offset;
			else
				performance.finalValue_performance = extrema.values[eLastVal] - offset;
		}

		if (settling_state == OVER_DAMPED)
		{
			performance.finalValue_performance = sim_data(0, sim_data.cols() - 1);
		}
	}
}

void StepResponseAnalyzer::solveOvershoot()
{	
	if (settling_state == SUFFICIENTLY_DAMPED || settling_state == UNDER_DAMPED)
	{
		double peak = *std::max_element(extrema.values.begin(), extrema.values.end());
		performance.delta_overshoot_performance = peak - performance.finalValue_performance;
		performance.percentOvershoot_performance = 100.0*(performance.delta_overshoot_performance) / performance.finalValue_performance;
	}

	/* The OVERDAMPED case does not have any overshoot */
	if (settling_state == OVER_DAMPED)
		performance.percentOvershoot_performance = 0.0;
}

void StepResponseAnalyzer::solveSettlingTime()
{
	/* Percent range +/- around final value that indicates settling */
	double settlingBoundary = 0.05;
	performance.settlingTime_window = settlingBoundary;

	if (settling_state == SUFFICIENTLY_DAMPED || settling_state == UNDER_DAMPED)
	{
		double delta = abs(settlingBoundary * performance.finalValue_performance);
		bool found = false;
		int time_index = 0;
		boost::container::vector<double> extreme_deltas;

		/* First, find the extrema that is closest to but just below the delta threshold */
		for (int i = 0; i < extrema.values.size(); i++)
		{
			extreme_deltas.push_back(abs(extrema.values[i] - performance.finalValue_performance));

			//Found a possible value
			if (extreme_deltas[i] < delta && found == false)
			{
				time_index = extrema.index[i];
				found = true;
			}

			//Oops. Went back outside the delta boundary at some point
			if (extreme_deltas[i] > delta && found == true)
			{
				time_index = 0;
				found = false;
			}
		}

		/* Now work backwards until you find the time step that corresponds
		 * to when you first enter the delta range
		 */
		for (int i = time_index; 0 < i; --i)
		{
			if (abs(sim_data(0, i) - performance.finalValue_performance) > delta)
			{
				time_index = i + 1;
				break;
			}
		}

		double dt = sim_data(1, (sim_data.cols() - 1)) / sim_data.cols();
		performance.settlingTime_performance = time_index*dt;
		performance.settlingTime_Idx = time_index;
	}

	/* The over damped case does not have settling time. Might need to change
	this value in practice though. */
	if (settling_state == OVER_DAMPED)
		performance.settlingTime_performance = 0.0;
}

void StepResponseAnalyzer::solveRiseTime()
{
	if (settling_state != NOT_CONVERGED)
	{
		double rise_start = 0.1*performance.finalValue_performance;
		double rise_end = 0.9*performance.finalValue_performance;

		double dt = sim_data(1, (sim_data.cols() - 1)) / sim_data.cols();

		double tStart = 0.0, tEnd = 0.0, target = rise_start;

		for (int i = 1; i < sim_data.cols(); i++)
		{
			if (abs(sim_data(0, i)) > abs(target))
			{
				if (tStart == 0)
				{
					tStart = sim_data(1, i - 1);
					target = abs(rise_end);
					performance.riseTime_Idx[0] = i - 1;
				}
				else if (tStart != 0 && tEnd == 0)
				{
					tEnd = sim_data(1, i - 1);
					performance.riseTime_Idx[1] = i - 1;
					break;
				}
			}
		}

		performance.riseTime_performance = tEnd - tStart;
	}
}

void StepResponseAnalyzer::solveSteadyStateError()
{
	if (settling_state != NOT_CONVERGED)
	{
		/* Because this is a step response analyzer, the output of this system 
		   SHOULD eventually converge to 1.0. This is our reference for error.*/

		performance.steadyStateError_performance = performance.finalValue_performance - 1.0;
	}
}