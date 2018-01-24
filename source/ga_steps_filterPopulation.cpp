#include "ga_steps_filterPopulation.h"

// void FCSOptimizer::filterPopulation()
// {
// 	#ifdef GA_TRACE_FILTER_POPULATION
// 	auto start = boost::chrono::high_resolution_clock::now();
// 	#endif
// 
// 	if (ga_instance_model_name == SSModelName)
// 	{
// 		#if defined(GA_CPU_SINGLE_THREADED) || defined(GA_CPU_MULTI_THREADED)
// 		/*-----------------------------------
// 		* Elite Performance Set
// 		*-----------------------------------*/
// 		if (GA_ElitistSolutions.EliteFitSolutions.size() < 10)
// 		{
// 			/* First time around, just update the variables directly */
// 			if (GA_BestFitnessValues.empty())
// 			{
// 				GA_ElitistSolutions.worst_performer_index = 0;
// 				GA_ElitistSolutions.worst_performer_value = 1.0;
// 			}
// 			/* Otherwise, push solution to the end */
// 			else
// 				GA_ElitistSolutions.EliteFitSolutions.push_back(GA_BestFitnessValues.back());
// 		}
// 
// 		/* Assuming the buffer is full, decide whether or not to add a solution in. */
// 		else
// 		{
// 			/* IF the latest "BestFitnessValue" is higher than the current worst fitness, this means 
// 			the Elitist solution set must be updated to include the new value. The goal is for the elite 
// 			set to update to higher and higher performing members each time, assuming such a member exists. */
// 			if (GA_ElitistSolutions.worst_performer_value < GA_BestFitnessValues.back().global_fitness)
// 			{
// 				GA_ElitistSolutions.EliteFitSolutions.data()[GA_ElitistSolutions.worst_performer_index] =
// 					GA_BestFitnessValues.back();
// 			}
// 		}
// 
// 		/* Find the new lowest performing individual out of the full set */
// 		double worst_performer = 2.0;
// 		for (int i = 0; i < GA_ElitistSolutions.EliteFitSolutions.size(); i++)
// 		{
// 			if (GA_ElitistSolutions.EliteFitSolutions.data()[i].global_fitness < worst_performer)
// 			{
// 				GA_ElitistSolutions.worst_performer_value = GA_ElitistSolutions.EliteFitSolutions.data()[i].global_fitness;
// 				GA_ElitistSolutions.worst_performer_index = i;
// 
// 				worst_performer = GA_ElitistSolutions.worst_performer_value;
// 			}
// 		}
// 
// 		/*-----------------------------------
// 		* Rejection of Population Members:
// 		* In this case, it is a simulation of natural disaster/survival of fittest
// 		*-----------------------------------*/
// 		boost::container::vector<int> rejectionIdxs;
// 
// 		/* Create a random number generator that is used to stepResponse how 
// 		   many population members are killed off at each generation. */
// 		int maxReplacements = (int)floor(hData.sim_data.population_size*0.4);
// 		std::random_device rd1;
// 		std::mt19937 rng1(rd1());
// 		std::uniform_int_distribution<uint32_t> uniform_int(0, maxReplacements);
// 
// 		/* Use a static threshold to filter through member performance */
// 		if (ga_instance_step_methods.filterType == GA_POPULATION_STATIC_FILTER)
// 		{
// 			int totalReplacements = uniform_int(rng1);
// 			double filter_threshold = 0.1;
// 
// 			for (int i = 0; i < totalReplacements; i++)
// 			{
// 				if (SS_FitnessValues.data()[i].global_fitness < filter_threshold)
// 					rejectionIdxs.push_back(i);
// 			}
// 		}
// 
// 		/* Dynamically choose a threshold to filter through member performance */
// 		if (ga_instance_step_methods.filterType == GA_POPULATION_DYNAMIC_FILTER)
// 		{
// 			std::random_device rd2;
// 			std::mt19937 rng2(rd2());
// 			std::uniform_real_distribution<double> uniform_dbl(0.1, 0.7);
// 
// 			int totalReplacements = uniform_int(rng1);
// 			double filter_threshold = uniform_dbl(rng2);
// 
// 			for (int i = 0; i < totalReplacements; i++)
// 			{
// 				if (SS_FitnessValues.data()[i].global_fitness < filter_threshold)
// 					rejectionIdxs.push_back(i);
// 			}
// 		}
// 
// 		/*-----------------------------------
// 		* Creation of New Members
// 		*-----------------------------------*/
// 		//Random, elite, copy of existing solutions?
// 
// 		double kpl = ga_instance_pid_config_data->tuningLimits.Kp_limits_lower;
// 		double kpu = ga_instance_pid_config_data->tuningLimits.Kp_limits_upper;
// 
// 		double kil = ga_instance_pid_config_data->tuningLimits.Ki_limits_lower;
// 		double kiu = ga_instance_pid_config_data->tuningLimits.Ki_limits_upper;
// 
// 		double kdl = ga_instance_pid_config_data->tuningLimits.Kd_limits_lower;
// 		double kdu = ga_instance_pid_config_data->tuningLimits.Kd_limits_upper;
// 
// 		/* Pure random generation of replacement members */
// 		for (int i = 0; i < rejectionIdxs.size(); i++)
// 		{
// 			hData.pid_data.Kp[rejectionIdxs[i]] = uniformRandomNumber(kpl, kpu);
// 			hData.pid_data.Ki[rejectionIdxs[i]] = uniformRandomNumber(kil, kiu);
// 			hData.pid_data.Kd[rejectionIdxs[i]] = uniformRandomNumber(kdl, kdu);
// 		}
// 		#endif
// 	}
// 
// 	#ifdef GA_TRACE_FILTER_POPULATION
// 	auto stop = boost::chrono::high_resolution_clock::now();
// 	#endif
// }


///////////////////////////////////////////////////
/* CLASS:  StaticFilter */
///////////////////////////////////////////////////
/*-----------------------------------------------
* Constructors/Destructor
*-----------------------------------------------*/
StaticFilter::StaticFilter(const double KpMax, const double KpMin, const double KiMax, const double KiMin, const double KdMax, const double KdMin)
{
	auto distKp = boost::random::uniform_real_distribution<>(KpMin, KpMax);
	rngKp = boost::make_shared <RNGInstance<boost::mt19937, boost::random::uniform_real_distribution<>>>(distKp);

	auto distKi = boost::random::uniform_real_distribution<>(KiMin, KiMax);
	rngKi = boost::make_shared <RNGInstance<boost::mt19937, boost::random::uniform_real_distribution<>>>(distKi);

	auto distKd = boost::random::uniform_real_distribution<>(KdMin, KdMax);
	rngKd = boost::make_shared <RNGInstance<boost::mt19937, boost::random::uniform_real_distribution<>>>(distKd);
}

StaticFilter::~StaticFilter()
{

}

/*-----------------------------------------------
* Public Functions
*-----------------------------------------------*/
void StaticFilter::filter(const GA_PopulationFilterDataInput input, GA_PopulationFilterDataOutput& output)
{
	double performanceThreshold = 0.0;
	if (input.static_performanceThreshold == 0.0)
		performanceThreshold = 0.1;
	else
		performanceThreshold = input.static_performanceThreshold;


	/* Generate new genetic material if a population member doesn't meet spec */
	rngKp->acquireEngine();
	rngKi->acquireEngine();
	rngKd->acquireEngine();

	for (int member = 0; member < input.currentGlobalFitScores.size(); member++)
	{
		if (input.currentGlobalFitScores[member] < performanceThreshold)
		{
			output.replacedMemberIndexes.push_back(member);
			output.replacementPIDValues.push_back({
				rngKp->getDouble(),
				rngKi->getDouble(),
				rngKd->getDouble()
			});
		}
	}

	rngKp->releaseEngine();
	rngKi->releaseEngine();
	rngKd->releaseEngine();
}

/*-----------------------------------------------
* Private Functions
*-----------------------------------------------*/



///////////////////////////////////////////////////
/* CLASS:  DynamicFilter */
///////////////////////////////////////////////////
/*-----------------------------------------------
* Constructors/Destructor
*-----------------------------------------------*/
DynamicFilter::DynamicFilter()
{

}

DynamicFilter::~DynamicFilter()
{

}

/*-----------------------------------------------
* Public Functions
*-----------------------------------------------*/
void DynamicFilter::filter(const GA_PopulationFilterDataInput input, GA_PopulationFilterDataOutput& output)
{

}

/*-----------------------------------------------
* Private Functions
*-----------------------------------------------*/


