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
TournamentSelection::TournamentSelection(const int populationSize)
{
	/* Ensure there are at least two possible competitor per tourney */
	auto dist1 = boost::random::uniform_int_distribution<>(2, populationSize - 1);
	tourneySizeSelectorRNG = boost::make_shared <RNGInstance<boost::mt19937, boost::random::uniform_int_distribution<>>>(dist1);

	/* Allow any members in the current population to be chosen without bias to compete */
	auto dist2 = boost::random::uniform_int_distribution<>(0, populationSize - 1);
	tourneyCompetitorSelectorRNG = boost::make_shared <RNGInstance<boost::mt19937, boost::random::uniform_int_distribution<>>>(dist2);
}

/*-----------------------------------------------
* Public Functions
*-----------------------------------------------*/
void TournamentSelection::selectParent(const GA_SelectParentDataInput input, GA_SelectParentDataOutput& output)
{
	tourneySizeSelectorRNG->acquireEngine();
	tourneyCompetitorSelectorRNG->acquireEngine();

	/*------------------------------
	* Run many tournaments to determine breedSelection parents
	*------------------------------*/
	for (int parent = 0; parent < input.populationSize; parent++)
	{
		/*------------------------------
		* Select competitor
		*------------------------------*/
		const int competitionSize = tourneySizeSelectorRNG->getInt();

		int* competitor = new int[competitionSize];

		/* Select the population members that will go head to head in competition this round! */
		for (int x = 0; x < competitionSize; x++)
			competitor[x] = tourneyCompetitorSelectorRNG->getInt();

		/*------------------------------
		* Evaluate the tournament
		*------------------------------*/
		double bestFitVal = 0.0;
		int bestFitIdx = 0;

		for (int x = 0; x < competitionSize; x++)
		{
			if (input.popGlobalFitScores[competitor[x]] > bestFitVal)
			{
				bestFitIdx = competitor[x];
				bestFitVal = input.popGlobalFitScores[competitor[x]];
			}
		}

		/* It's possible to get all -1.0 fitness values (i.e. nothing met criteria). If so, just
		* choose the last idx value in the competitor[] selections as the winner for now.
		*/
		if (bestFitVal == 0.0)
			bestFitIdx = competitor[competitionSize - 1];

		output.parentSelections[parent] = bestFitIdx;
		delete[] competitor;
	}

	tourneySizeSelectorRNG->releaseEngine();
	tourneyCompetitorSelectorRNG->releaseEngine();
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