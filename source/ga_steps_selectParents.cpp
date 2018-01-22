#include "ga_steps_selectParents.h"

///////////////////////////////////////////////////
/* CLASS:  RankedSelection */
///////////////////////////////////////////////////
/*-----------------------------------------------
* Constructors/Destructor
*-----------------------------------------------*/
RankedSelection::RankedSelection()
{
}

RankedSelection::~RankedSelection()
{
}
/*-----------------------------------------------
* Public Functions
*-----------------------------------------------*/
void RankedSelection::selectParent(const GA_SelectParentDataInput input, GA_SelectParentDataOutput& output)
{

}

/*-----------------------------------------------
* Private Functions
*-----------------------------------------------*/
void RankedSelection::selectParentKp()
{
	
}

void RankedSelection::selectParentKi()
{
	
}

void RankedSelection::selectParentKd()
{
	
}


///////////////////////////////////////////////////
/* CLASS:  RandomSelection */
///////////////////////////////////////////////////
/*-----------------------------------------------
* Constructors/Destructor
*-----------------------------------------------*/
RandomSelection::RandomSelection(const int populationSize)
{
	auto distribution = boost::random::uniform_int_distribution<>(0, populationSize - 1);
	rng_engine = boost::make_shared<RNGInstance<boost::mt19937, boost::random::uniform_int_distribution<>>>(distribution);
}

RandomSelection::~RandomSelection()
{
}
/*-----------------------------------------------
* Public Functions
*-----------------------------------------------*/
void RandomSelection::selectParent(const GA_SelectParentDataInput input, GA_SelectParentDataOutput& output)
{
	rng_engine->acquireEngine();
	for (uint32_t parent = 0; parent < input.populationSize; parent++)
		output.parentSelections[parent] = rng_engine->getInt();
	rng_engine->releaseEngine();
}

/*-----------------------------------------------
* Private Functions
*-----------------------------------------------*/
void RandomSelection::selectParentKp()
{

}

void RandomSelection::selectParentKi()
{

}

void RandomSelection::selectParentKd()
{

}


///////////////////////////////////////////////////
/* CLASS:  RouletteSelection */
///////////////////////////////////////////////////
/*-----------------------------------------------
* Constructors/Destructor
*-----------------------------------------------*/
RouletteSelection::RouletteSelection()
{
}

RouletteSelection::~RouletteSelection()
{
}
/*-----------------------------------------------
* Public Functions
*-----------------------------------------------*/
void RouletteSelection::selectParent(const GA_SelectParentDataInput input, GA_SelectParentDataOutput& output)
{

}

/*-----------------------------------------------
* Private Functions
*-----------------------------------------------*/
void RouletteSelection::selectParentKp()
{

}

void RouletteSelection::selectParentKi()
{

}

void RouletteSelection::selectParentKd()
{

}


///////////////////////////////////////////////////
/* CLASS:  StochasticSelection */
///////////////////////////////////////////////////
/*-----------------------------------------------
* Constructors/Destructor
*-----------------------------------------------*/
StochasticSelection::StochasticSelection()
{
}

StochasticSelection::~StochasticSelection()
{
}
/*-----------------------------------------------
* Public Functions
*-----------------------------------------------*/
void StochasticSelection::selectParent(const GA_SelectParentDataInput input, GA_SelectParentDataOutput& output)
{

}

/*-----------------------------------------------
* Private Functions
*-----------------------------------------------*/
void StochasticSelection::selectParentKp()
{

}

void StochasticSelection::selectParentKi()
{

}

void StochasticSelection::selectParentKd()
{

}


///////////////////////////////////////////////////
/* CLASS:  TournamentSelection */
///////////////////////////////////////////////////
/*-----------------------------------------------
* Constructors/Destructor
*-----------------------------------------------*/
TournamentSelection::TournamentSelection()
{
}

TournamentSelection::~TournamentSelection()
{
}
/*-----------------------------------------------
* Public Functions
*-----------------------------------------------*/
void TournamentSelection::selectParent(const GA_SelectParentDataInput input, GA_SelectParentDataOutput& output)
{

}

/*-----------------------------------------------
* Private Functions
*-----------------------------------------------*/
void TournamentSelection::selectParentKp()
{

}

void TournamentSelection::selectParentKi()
{

}

void TournamentSelection::selectParentKd()
{

}


///////////////////////////////////////////////////
/* CLASS:  ElitistSelection */
///////////////////////////////////////////////////
/*-----------------------------------------------
* Constructors/Destructor
*-----------------------------------------------*/
ElitistSelection::ElitistSelection()
{
}

ElitistSelection::~ElitistSelection()
{
}
/*-----------------------------------------------
* Public Functions
*-----------------------------------------------*/
void ElitistSelection::selectParent(const GA_SelectParentDataInput input, GA_SelectParentDataOutput& output)
{

}

/*-----------------------------------------------
* Private Functions
*-----------------------------------------------*/
void ElitistSelection::selectParentKp()
{

}

void ElitistSelection::selectParentKi()
{

}

void ElitistSelection::selectParentKd()
{

}



// void RandomSelection::selectParents()
// {
// // 	rs_selections = out_selections;
// // 
// // 	/* Unless there is a REALLY large population, splitting this into threads probably
// // 	* is not very beneficial */
// // 	if (executionType == SINGLE_THREADED || executionType == MULTI_THREADED)
// // 		calculate_cpu_single_threaded();
// // 
// // 	if (executionType == SINGLE_THREADED_WITH_CUDA || executionType == MULTI_THREADED_WITH_CUDA)
// // 		calculate_gpu_single_threaded();
// }
// 
// /*-----------------------------------------------
// * Private Functions
// *-----------------------------------------------*/
// void RandomSelection::calculate_cpu_single_threaded()
// {
// 	std::random_device rd;
// 	std::mt19937 rng(rd());
// 	std::uniform_int_distribution<uint32_t> uniform_int(0, rs_selections->size() - 1);
// 
// 	/* Randomly select a mate out of the current gene pool */
// 	for (uint32_t parent = 0; parent < rs_selections->size(); parent++)
// 		rs_selections->data()[parent] = uniform_int(rng);
// }
// 
// void RandomSelection::calculate_gpu_single_threaded()
// {
// 	std::cout << "GPU Single/Multi Threaded Mode Currently Not Supported." << std::endl;
// }

// /*-----------------------------------------------
// * Public Functions
// *-----------------------------------------------*/
// void TournamentSelection::selectParents(PIDFitness_Vec in_fitnessValues, iVec* out_selections)
// {
// 	ts_data = in_fitnessValues;
// 	ts_selections = out_selections;
// 
// // 	if (executionType == SINGLE_THREADED)
// // 		calculate_cpu_single_threaded();
// // 
// // 	if (executionType == MULTI_THREADED)
// // 		calculate_cpu_multi_threaded();
// // 
// // 	if (executionType == SINGLE_THREADED_WITH_CUDA)
// // 		calculate_gpu_single_threaded();
// // 
// // 	if (executionType == MULTI_THREADED_WITH_CUDA)
// // 		calculate_gpu_multi_threaded();
// }
// 
// /*-----------------------------------------------
// * Private Functions
// *-----------------------------------------------*/
// void TournamentSelection::calculate_cpu_single_threaded()
// {
// 	size_t popSize = ts_data.size();
// 
// 	/*------------------------------
// 	* Create RNG
// 	*------------------------------*/
// 	std::random_device rd;
// 	std::mt19937 rng(rd());
// 
// 	/* Ensure there are at least two competitors per tourney*/
// 	std::uniform_int_distribution<uint32_t> uniform_num_competitors(2, popSize - 1);
// 
// 	/* Allow any members in the current population to be chosen without bias to compete */
// 	std::uniform_int_distribution<uint32_t> uniform_competitor(0, popSize - 1);
// 
// 	/*------------------------------
// 	* Run many tournaments to determine breedSelection parents
// 	*------------------------------*/
// 	for (int parent = 0; parent < popSize; parent++)
// 	{
// 		/*------------------------------
// 		* Select competitors
// 		*------------------------------*/
// 		size_t num_competitors = uniform_num_competitors(rng);
// 		iVec competitors(num_competitors);
// 
// 		for (int x = 0; x < num_competitors; x++)
// 			competitors[x] = uniform_competitor(rng);
// 
// 		/*------------------------------
// 		* Evaluate the tournament
// 		*------------------------------*/
// 		double bestFitVal = 0.0;
// 
// 		for (int x = 0; x < num_competitors; x++)
// 		{
// 			uint32_t competitor_idx = competitors[x];
// 			double currentFitVal = ts_data.data()[competitor_idx].global_fitness;
// 
// 			//Maximization problem
// 			if (currentFitVal > bestFitVal)
// 				ts_selections->data()[parent] = competitor_idx;
// 		}
// 
// 		/* It's possible to get all -1.0 fitness values (i.e. nothing met criteria). If so, just
// 		* choose the last idx value in the competitors[] selections as the winner for now.
// 		*/
// 		if (bestFitVal == -1.0)
// 		{
// 			ts_selections->data()[parent] = competitors[num_competitors - 1];
// 
// 			//TODO: Select based on the best fitness in each sub category
// 		}
// 	}
// }