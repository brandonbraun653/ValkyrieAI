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

		//No need to process because there is nothing to sort!
		if (frontMemberLength == 1)
			continue;

		
		int numObjectives = 4; //TODO: Need to specify this somehow in a global config file
		
		//Intermediary container for holding the sorted 
		boost::container::vector<int> sortedIdxs; 
		for (int objective = 0; objective < numObjectives; objective++)
		{
			if (F[front].size() == 0)
				break;

			/* Re-sorts the population members in a given front in ascending order
			according to how well they performed for an objective function. */ 
			sortedIdxs = sortObjectiveFunc(F[front], memberObjFuncFit[objective]);

			/* Ensure boundary points are selected so there is a relative reference
			to compare crowding against. Because the indexes are already sorted, it
			is only a matter of selecting the first and last indices. */
			int first = *sortedIdxs.begin();
			int last = *(sortedIdxs.end()-1);

			memberDistance[first] = DBL_MAX;	//Technically should be infinity, but ya know...
			memberDistance[last] = DBL_MAX;

			/* Calculate the distance between all other points in the front, using 
			the first and last points chosen above as a reference. The official 
			algorithm is also normalized, but the max/min values are 1.0 and 0.0, so the 
			division by 1.0 is pointless to perform.*/
			for (int i = 1; i < frontMemberLength - 1; i++)
			{
				int previous = sortedIdxs[i - 1];
				int current = sortedIdxs[i];
				int next = sortedIdxs[i + 1];

				memberDistance[current] = memberDistance[current] + 
					(memberObjFuncFit[objective][next] - memberObjFuncFit[objective][previous]); //  / (fmax - fmin);
			}
		}
	}

	
	/*-------------------------------------
	* New Order by Crowding Distance
	*------------------------------------*/
	for (int front = 0; front < F.size(); front++)
	{
		int frontMemberLength = F[front].size();

		//Nothing to compete against. Add directly to the output
		if (frontMemberLength == 1)
		{
			output.sortedPopulation.push_back(F[front][0]);
			continue;
		}

		/* Sorts all the population members in this front in descending order
		according their crowding distance. The idea is that to maintain diversity, solutions
		that are farther apart from each other are better than those closer together. */
		auto sortedIdxs = sortCrowdingDistance(F[front], memberDistance);
		
		/* Add the sorted data to the end of the output buffer*/
		output.sortedPopulation.insert(output.sortedPopulation.end(), sortedIdxs.begin(), sortedIdxs.end());
	}
	
}

//TODO: build a custom sorter functor that gets the index along with the score
boost::container::vector<int> FastNonDominatedSort::sortObjectiveFunc(
	boost::container::vector<int> memberSet,
	boost::container::vector<double> objFuncScores)
{
	boost::container::vector<int> output;

	//First, build the vector containing all the relevant scores
	boost::container::vector<double> scores;
	for (int i = 0; i < memberSet.size(); i++)
		scores.push_back(objFuncScores[memberSet[i]]);

	//Sort the scores in ASCENDING order
	std::sort(scores.begin(), scores.end());

	//Rebuild the output vector by matching sorted scores to input indexes
	for (int j = 0; j < scores.size(); j++)
	{
		for (int i = 0; i < memberSet.size(); i++)
		{
			if (objFuncScores[memberSet[i]] == scores[j])
			{
				output.push_back(memberSet[i]);
				break; //should only push out of the inner loop
			}
		}
	}

	return output;
}


//TODO: build a custom sorter functor that gets the index along with the score
boost::container::vector<int> FastNonDominatedSort::sortCrowdingDistance(
	boost::container::vector<int> memberSet,
	boost::container::vector<double> crowdingDistances)
{
	boost::container::vector<int> output;

	//First, build the vector containing all the relevant scores
	boost::container::vector<double> distances;
	for (int i = 0; i < memberSet.size(); i++)
		distances.push_back(crowdingDistances[memberSet[i]]);

	//Sort the scores in DESCENDING order
	std::sort(distances.rbegin(), distances.rend());

	//Rebuild the output vector by matching sorted scores to input indexes
	for (int j = 0; j < distances.size(); j++)
	{
		for (int i = 0; i < memberSet.size(); i++)
		{
			if (crowdingDistances[memberSet[i]] == distances[j])
			{
				output.push_back(memberSet[i]);

				/* Prevents accidental inclusion of the same population member when 
				several members share the same crowding distance. This occurs if 
				a set of PID values really don't work and meet none of the user's goals.*/
				memberSet.erase(memberSet.begin() + i);
			}
		}
	}

	return output;
}