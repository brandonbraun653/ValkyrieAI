#include "ga_mop_steps_selectParents.h"

///////////////////////////////////////////////////
/* CLASS:  RankedSelection */
///////////////////////////////////////////////////
/*-----------------------------------------------
* Constructors/Destructor
*-----------------------------------------------*/
RankedSelection::RankedSelection(GA_RunMode execution_type)
{
	executionType = execution_type;
}

RankedSelection::~RankedSelection()
{
}
/*-----------------------------------------------
* Public Functions
*-----------------------------------------------*/

/*-----------------------------------------------
* Private Functions
*-----------------------------------------------*/

///////////////////////////////////////////////////
/* CLASS:  RandomSelection */
///////////////////////////////////////////////////
/*-----------------------------------------------
* Constructors/Destructor
*-----------------------------------------------*/
RandomSelection::RandomSelection(GA_RunMode execution_type)
{
	executionType = execution_type;
}

RandomSelection::~RandomSelection()
{
}

/*-----------------------------------------------
* Public Functions
*-----------------------------------------------*/
void RandomSelection::selectParents(iVec* out_selections)
{
	rs_selections = out_selections;

	/* Unless there is a REALLY large population, splitting this into threads probably
	* is not very beneficial */
	if (executionType == SINGLE_THREADED || executionType == MULTI_THREADED)
		calculate_cpu_single_threaded();

	if (executionType == SINGLE_THREADED_WITH_CUDA || executionType == MULTI_THREADED_WITH_CUDA)
		calculate_gpu_single_threaded();
}

/*-----------------------------------------------
* Private Functions
*-----------------------------------------------*/
void RandomSelection::calculate_cpu_single_threaded()
{
	std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_int_distribution<uint32_t> uniform_int(0, rs_selections->size() - 1);

	/* Randomly select a mate out of the current gene pool */
	for (uint32_t parent = 0; parent < rs_selections->size(); parent++)
		rs_selections->data()[parent] = uniform_int(rng);
}

void RandomSelection::calculate_gpu_single_threaded()
{
	std::cout << "GPU Single/Multi Threaded Mode Currently Not Supported." << std::endl;
}

///////////////////////////////////////////////////
/* CLASS:  TournamentSelection */
///////////////////////////////////////////////////
/*-----------------------------------------------
* Constructors/Destructor
*-----------------------------------------------*/
TournamentSelection::TournamentSelection(GA_RunMode execution_type)
{
	executionType = execution_type;
}

TournamentSelection::~TournamentSelection()
{
}

/*-----------------------------------------------
* Public Functions
*-----------------------------------------------*/
void TournamentSelection::selectParents(PIDFitness_Vec in_fitnessValues, iVec* out_selections)
{
	ts_data = in_fitnessValues;
	ts_selections = out_selections;

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
void TournamentSelection::calculate_cpu_single_threaded()
{
	size_t popSize = ts_data.size();

	/*------------------------------
	* Create RNG
	*------------------------------*/
	std::random_device rd;
	std::mt19937 rng(rd());

	/* Ensure there are at least two competitors per tourney*/
	std::uniform_int_distribution<uint32_t> uniform_num_competitors(2, popSize - 1);

	/* Allow any members in the current population to be chosen without bias to compete */
	std::uniform_int_distribution<uint32_t> uniform_competitor(0, popSize - 1);

	/*------------------------------
	* Run many tournaments to determine breedSelection parents
	*------------------------------*/
	for (int parent = 0; parent < popSize; parent++)
	{
		/*------------------------------
		* Select competitors
		*------------------------------*/
		size_t num_competitors = uniform_num_competitors(rng);
		iVec competitors(num_competitors);

		for (int x = 0; x < num_competitors; x++)
			competitors[x] = uniform_competitor(rng);

		/*------------------------------
		* Evaluate the tournament
		*------------------------------*/
		double bestFitVal = 0.0;

		for (int x = 0; x < num_competitors; x++)
		{
			uint32_t competitor_idx = competitors[x];
			double currentFitVal = ts_data.data()[competitor_idx].global_fitness;

			//Maximization problem
			if (currentFitVal > bestFitVal)
				ts_selections->data()[parent] = competitor_idx;
		}

		/* It's possible to get all -1.0 fitness values (i.e. nothing met criteria). If so, just
		* choose the last idx value in the competitors[] selections as the winner for now.
		*/
		if (bestFitVal == -1.0)
		{
			ts_selections->data()[parent] = competitors[num_competitors - 1];

			//TODO: Select based on the best fitness in each sub category
		}
	}
}

void TournamentSelection::calculate_cpu_multi_threaded()
{
	std::cout << "CPU Multi Threaded Mode Currently Not Supported." << std::endl;
}

void TournamentSelection::calculate_gpu_single_threaded()
{
	std::cout << "GPU Single Threaded Mode Currently Not Supported." << std::endl;
}

void TournamentSelection::calculate_gpu_multi_threaded()
{
	std::cout << "GPU Multi Threaded Mode Currently Not Supported." << std::endl;
}