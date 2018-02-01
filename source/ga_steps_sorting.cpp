#include "ga_steps_sorting.h"

struct DominationData
{
	int dominationCount = 0;
	boost::container::vector<int> dominationSet;
};

///////////////////////////////////////////////////
/* CLASS:  FastNonDominatedSort */
///////////////////////////////////////////////////
void FastNonDominatedSort::sort(const GA_SortingInput input, GA_SortingOutput& output)
{
	/* Combine the two populations into one. The parent input is forced to exist, so
	there only needs to be a check for the child input. */
	boost::container::vector<double> R = input.parentChildFitScores;
	boost::container::vector<boost::container::vector<int>> F;
	
	F.resize(R.size() + 1);
	int frontNum = 0;


	/* Each population member needs to externally keep track of how many times 
	it has been dominated AND which members it dominates */
	boost::container::vector<int> memberDominationCount;
	boost::container::vector<boost::container::vector<int>> memberDominationSet;

	memberDominationSet.resize(R.size());
	memberDominationCount.resize(R.size());

	/*-------------------------------------
	* Find all the population members that belong to the first front
	*------------------------------------*/
	for (int p = 0; p < R.size(); p++)
	{
		/* Reset variables for next iteration */
		memberDominationCount[p] = 0;
		memberDominationSet[p].clear();


		for (int q = 0; q < R.size(); q++)
		{
			//Assuming maximization...

			/* If i dominates q, add q to the dominated solution set for population member i*/
			if (R[p] > R[q])
				memberDominationSet[p].push_back(q);

			/* If q dominates i, i cannot be the best solution, so increase i's domination counter */
			else if (R[p] < R[q])
				memberDominationCount[p]++;
		}

		/* Assuming Np is 0, this means we found a solution that was
		not dominated by other population members */
		if (memberDominationCount[p] == 0)
		{
			//Add p to the first front of non-dominated solutions 
			F[frontNum].push_back(p);
		}
	}

	/*-------------------------------------
	* Find all the population members that belong to the remaining fronts
	*------------------------------------*/
	while (F[frontNum].size() != 0)
	{
		boost::container::vector<int> Q;

		for (int p = 0; p < F[frontNum].size(); p++)
		{
			int member = F[frontNum][p];

			for (int q = 0; q < memberDominationSet[member].size(); q++)
			{
				//Grab the dominated population member
				int memberIdx = memberDominationSet[member][q];

				//Decrement its domination count 
				memberDominationCount[memberIdx]--;

				//If it becomes zero, it belongs to then next front
				if (memberDominationCount[memberIdx] == 0)
					Q.push_back(memberIdx);
			}
		}

		frontNum++;
		F[frontNum] = Q;
	}


	
	for (int i = 0; i < F.size(); i++)
	{
		for (int j = 0; j < F[i].size(); j++)
		{
			output.sortedPopulation.push_back(F[i][j]);
		}
	}
}

