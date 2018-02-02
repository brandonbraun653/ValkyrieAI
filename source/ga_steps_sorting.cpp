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
	boost::container::vector<double> R;							/* Combined population */
	boost::container::vector<boost::container::vector<int>> F;	/* Ranking of population into fronts */
	
	/* Force 1 decimal point resolution on the fitness scores */
	for (int i = 0; i < input.parentChildFitScores.size(); i++)
	{
		R.push_back(input.parentChildFitScores[i].fitness_total);

		double fracPart = 0.0;
		double intPart = 0.0;

		/* Decompose the data into integral and fractional parts */
		fracPart = std::modf(R[i], &intPart);

		/* Shift up, truncate, shift down data */
		fracPart *= std::pow(10.0, (int)GA_RESOLUTION_1DP);
		fracPart = floor(fracPart);
		fracPart /= std::pow(10.0, (int)GA_RESOLUTION_1DP);

		/* Reassign truncated value */
		R[i] = intPart + fracPart;
	}
		


	F.resize(R.size() + 1);
	int frontNum = 0;


	/* Each population member needs to externally keep track of how many times 
	it has been dominated AND which members it dominates */
	boost::container::vector<int> memberRank;
	boost::container::vector<int> memberDominationCount;
	boost::container::vector<double> memberDistance;
	boost::container::vector<boost::container::vector<int>> memberDominationSet;


	/* Oh gosh this is disgusting */
	boost::container::vector<boost::container::vector<double>> memberObjFuncFit;
	memberObjFuncFit.resize(4);

	for (int i = 0; i < input.parentChildFitScores.size(); i++)
		memberObjFuncFit[0].push_back(input.parentChildFitScores[i].fitness_POS);

	for (int i = 0; i < input.parentChildFitScores.size(); i++)
		memberObjFuncFit[1].push_back(input.parentChildFitScores[i].fitness_SSER);

	for (int i = 0; i < input.parentChildFitScores.size(); i++)
		memberObjFuncFit[2].push_back(input.parentChildFitScores[i].fitness_TR);

	for (int i = 0; i < input.parentChildFitScores.size(); i++)
		memberObjFuncFit[3].push_back(input.parentChildFitScores[i].fitness_TS);



	memberDominationSet.resize(R.size());
	memberDominationCount.resize(R.size());
	memberDistance.resize(R.size());
	memberRank.resize(R.size());

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
			memberRank[p] = 0;
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
				{
					Q.push_back(memberIdx);
					memberRank[memberIdx] = frontNum + 1;
				}
					
			}
		}

		frontNum++;
		F[frontNum] = Q;
	}

	/*-------------------------------------
	* Crowding Distance Assignment
	*------------------------------------*/
	for (int front = 0; front < F.size(); front++)
	{
		int frontMemberLength = F[front].size();

		//TODO: Need to specify this somehow in a global config file
		int numObjectives = 4;
		boost::container::vector<int> I;
		for (int objective = 0; objective < numObjectives; objective++)
		{
			if (F[front].size() == 0)
				break;

			//re-sort the front based on objective function fitness 
			I = sortObjectiveFunc(F[front], memberObjFuncFit[objective]);

			//Ensure boundary points are selected
			memberDistance[I[0]] = DBL_MAX;
			memberDistance[I[I.size()-1]] = DBL_MAX;

			//Calculate the distance for all other points 
			for (int i = 1; i < frontMemberLength - 1; i++)
			{
				memberDistance[I[i]] = memberDistance[I[i]] +
					(memberObjFuncFit[objective][I[i + 1]] - memberObjFuncFit[objective][I[i - 1]]);
			}
		}
	}

	
	/*-------------------------------------
	* New Order by Crowding Distance
	*------------------------------------*/
	for (int i = 0; i < R.size(); i++)
	{
		if (memberRank[i] == 0)
		{
			output.sortedPopulation.push_back(i);
			memberRank[i] = -1; //To say we already selected it
		}
		else
		{
			int lowestRankIndex = 0;
			for (int j = 0; j < R.size(); j++)
			{
				if (memberRank[i] < memberRank[j])
					lowestRankIndex = i;

				else if ((memberRank[i] == memberRank[j]) && (memberDistance[i] > memberDistance[j]))
					lowestRankIndex = i;
			}

			output.sortedPopulation.push_back(lowestRankIndex);
			memberRank[lowestRankIndex] = -1;
		}
	}

	
}

boost::container::vector<int> FastNonDominatedSort::sortObjectiveFunc(
	boost::container::vector<int> memberSet,
	boost::container::vector<double> objFuncScores)
{
	boost::container::vector<int> output;

	//First, build the vector containing all the relevant scores
	boost::container::vector<double> scores;
	for (int i = 0; i < memberSet.size(); i++)
		scores.push_back(objFuncScores[memberSet[i]]);


	//Now sort those scores in ascending order
	std::sort(scores.begin(), scores.end());

	//Rebuild the output vector by matching sorted scores to input indexes
	for (int j = 0; j < scores.size(); j++)
	{
		for (int i = 0; i < memberSet.size(); i++)
		{
			if (objFuncScores[memberSet[i]] == scores[j])
			{
				output.push_back(memberSet[i]);

				break; //should only push out of the inner loop?
			}
		}
	}

	return output;
}

