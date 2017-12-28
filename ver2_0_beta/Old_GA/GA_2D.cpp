//////////////////////////////////////////////////////////////////
/* CLASS: GACLASS_2D */
//////////////////////////////////////////////////////////////////

/*-----------------------------------------------
* Constructors/Destructors
*-----------------------------------------------*/
GAClass_2D::GAClass_2D()
{
	algStatus = GA_OK;
	algRuntimeData.bestFit = DBL_MAX;
	algRuntimeData.best_x = DBL_MAX;
	algRuntimeData.best_y = DBL_MAX;
	algRuntimeData.iterationNumber = 0;
}

GAClass_2D::~GAClass_2D()
{
	/* Free up dynamically allocated arrays */
	for (uint32_t i = 0; i < algControl.populationSize; i++)
	{
		delete[] currentPopulation[i];
		delete[] breedSelections[i];
	}

	delete[] currentPopulation;
	delete[] breedSelections;
}

/*-----------------------------------------------
* Public Functions
*-----------------------------------------------*/
void GAClass_2D::assignAlgorithm(GA2D_FunctionBenchMark benchMark)
{
	assignFitnessFunction(benchMark.function);
	assignSearchConstraint(benchMark.limits);

	algBreedStyle = benchMark.breedStyle;
	algMateSelectStyle = benchMark.mateSelectStyle;
	algMutationStyle = benchMark.mutationStyle;
}

void GAClass_2D::localSearch(double x_center, double y_center, uint32_t axis_ticks,
	double *results, double percentage)
{
	double dX = abs(percentage * x_center);
	double dY = abs(percentage * y_center);

	if (dX == 0.0)
		dX = 0.15;
	if (dY == 0.0)
		dY = 0.15;

	std::vector<double> xPts = linspace((x_center - dX), (x_center + dX), axis_ticks);
	std::vector<double> yPts = linspace((y_center - dY), (y_center + dY), axis_ticks);

	/* Linearly search along the coordinate axis */
	double lowestVal = DBL_MAX;
	double currentVal = 0.0;
	for (uint32_t x = 0; x < xPts.size(); x++)
	{
		for (uint32_t y = 0; y < yPts.size(); y++)
		{
			currentVal = fitnessFunction(xPts[x], yPts[y]);

			if (currentVal < lowestVal)
			{
				results[X] = xPts[x];
				results[Y] = yPts[y];
				results[F] = currentVal;

				lowestVal = currentVal;
			}
		}
	}
}

GA_StatusTypeDef GAClass_2D::begin(uint32_t populationSize, uint32_t generationLimit, double mutationProbabilityThreshold)
{
	/* Copy in the general control information and initialize the algorithm */
	algControl.populationSize = populationSize;
	algControl.generationLimit = generationLimit;
	algControl.mutation_rate = mutationProbabilityThreshold;

	return initGeneration();
}

GA_StatusTypeDef GAClass_2D::reset()
{
	algStatus = GA_READY;
	generation_counter = 0;

	/*-------------------------------------------
	* Clear out the runtime logs
	*-------------------------------------------*/
	algRuntimeData.bestFit = DBL_MAX;
	algRuntimeData.best_x = DBL_MAX;
	algRuntimeData.best_y = DBL_MAX;

	/*-------------------------------------------
	* Completely randomize the chromosomes
	*-------------------------------------------*/
	for (uint32_t i = 0; i < algControl.populationSize; i++)
	{
		currentPopulation[i][X] = uniformRandomNumber(algConstraints.x_lower, algConstraints.x_upper);
		currentPopulation[i][Y] = uniformRandomNumber(algConstraints.y_lower, algConstraints.y_upper);
	}

	return algStatus;
}

GA_StatusTypeDef GAClass_2D::iterateGeneration()
{
	algStatus = GA_OK;
	generation_counter += 1;

	/* Execute generational simulation */
	evaluateFitness();
	selectMates(algMateSelectStyle);
	breedGeneration(algBreedStyle);
	mutateChromosomes(algMutationStyle);

	/* Stopping Condition Check*/
	if (generation_counter == algControl.generationLimit)
		algStatus = GA_COMPLETE;

	return algStatus;
}

GA_StatusTypeDef GAClass_2D::peekStatus(GA2D_RuntimeData *container)
{
	if (algRuntimeData.locked)
		return GA_BUSY;

	//TODO: Can this be copied easier than this explicit non-sense?
	//		I don't want a pointer to the real data...I want a copy.
	container->iterationNumber = algRuntimeData.iterationNumber;
	container->bestFit = algRuntimeData.bestFit;
	container->best_x = algRuntimeData.best_x;
	container->best_y = algRuntimeData.best_y;

	return algStatus;
}

GA_StatusTypeDef GAClass_2D::peekStatus()
{
	return algStatus;
}

/*-----------------------------------------------
* Private Functions
*-----------------------------------------------*/
void GAClass_2D::assignFitnessFunction(fitfunc_2d_t fitfunc)
{
	fitnessFunction = fitfunc;
}

void GAClass_2D::assignSearchConstraint(GA2D_Constraints problemLimits)
{
	algConstraints = problemLimits;

	/*-------------------------------------------
	* Handle possible div/0 error
	*-------------------------------------------*/
	if (algConstraints.x_gridDiv == 0 || algConstraints.y_gridDiv == 0)
	{
		std::cout << "Divide by zero condition created from given SearchConstraint container!\nOverriding with GridSize of 100.";
		algConstraints.x_gridDiv = 100;
		algConstraints.y_gridDiv = 100;
	}

	/*-------------------------------------------
	* Calculate Mapping Coefficients (x,y)
	*-------------------------------------------*/
	bytes_precision = (double)algConstraints.numBytes_precision;

	x_lo = algConstraints.x_lower;	y_lo = algConstraints.y_lower;
	x_hi = algConstraints.x_upper;	y_hi = algConstraints.y_upper;

	x_dPow = 8 * bytes_precision - 1;	y_dPow = x_bPow;
	x_bPow = ceil(log2(x_hi - x_lo));	y_bPow = ceil(log2(y_hi - y_lo));

	x_sF = pow(2, (x_dPow - x_bPow));	y_sF = pow(2, (y_dPow - y_bPow));
	x_sR = pow(2, (x_bPow - x_dPow));	y_sR = pow(2, (y_bPow - y_dPow));

	x_offset = 0.0;	if (x_lo < 0.0 && x_hi > 0.0) x_offset = (x_sF*x_lo) - floor(x_sF*x_lo);
	y_offset = 0.0; if (y_lo < 0.0 && y_hi > 0.0) y_offset = (y_sF*y_lo) - floor(y_sF*y_lo);
}

GA_StatusTypeDef GAClass_2D::initGeneration()
{
	const uint32_t popSize = algControl.populationSize;

	/*-------------------------------------------
	* Create containers to hold population data
	*-------------------------------------------*/
	currentPopulation = new double *[popSize];
	breedSelections = new uint16_t *[popSize];
	currentFitness = new double[popSize];

	/* Expand the arrays to hold one set of (x,y) coordinates per population member */
	for (uint32_t popMember = 0; popMember < popSize; popMember++)
	{
		const uint8_t numerical_dimension = 2u;
		currentPopulation[popMember] = new double[numerical_dimension];
		breedSelections[popMember] = new uint16_t[numerical_dimension];
	}

	/*-------------------------------------------
	* Fill in the initial chromosomes with random data
	*-------------------------------------------*/
	for (uint32_t i = 0; i < popSize; i++)
	{
		currentPopulation[i][0] = uniformRandomNumber(algConstraints.x_lower, algConstraints.x_upper);
		currentPopulation[i][1] = uniformRandomNumber(algConstraints.y_lower, algConstraints.y_upper);
	}

	algStatus = GA_OK;
	return algStatus;
}

GA_StatusTypeDef GAClass_2D::evaluateFitness()
{
	double bestFit = DBL_MAX, bestX = DBL_MAX, bestY = DBL_MAX;

	//NOTE: Parallelization opportunity
	for (uint32_t member = 0; member < algControl.populationSize; member++)
	{
		/* Calculate fitness */
		double x = currentPopulation[member][X];
		double y = currentPopulation[member][Y];
		currentFitness[member] = fitnessFunction(x, y);

		/* Check for local best solution */
		if (currentFitness[member] < bestFit)
		{
			bestFit = currentFitness[member];
			bestX = x;
			bestY = y;
		}
	}

	/* Check for global best solution */
	if (bestFit < algRuntimeData.bestFit)
	{
		algRuntimeData.bestFit = bestFit;
		algRuntimeData.best_x = bestX;
		algRuntimeData.best_y = bestY;
	}

	return algStatus;
}

GA_StatusTypeDef GAClass_2D::selectMates(GA_SelectStyle_TypeDef selection_approach)
{
	/*NOTE:
	1. Regardless of the selection process chosen, it is always assumed that there are two
	parents and two children.

	2. The chosen parents are stored in sequential pairs in breedSelections. i.e. [1] & [2] mate,
	[3] & [4] mate, [5] & [6] mate, etc.

	3. The data format for breedSelections is a uint16_t array indexed by [parent][coordinate axis]
	*/

	if (selection_approach == GA_SELECT_RANDOM)
	{
		std::random_device rd;
		std::mt19937 rng(rd());
		std::uniform_int_distribution<uint32_t> uniform_int(0, algControl.populationSize - 1);

		/* Randomly select a mate out of the current gene pool */
		for (uint32_t i = 0; i < algControl.populationSize; i++)
		{
			uint32_t mate = uniform_int(rng);
			breedSelections[i][X] = coordinate2Chromosome(currentPopulation[mate][X], X);
			breedSelections[i][Y] = coordinate2Chromosome(currentPopulation[mate][Y], Y);
		}
	}

	if (selection_approach == GA_SELECT_RANKED)
	{
		//Don't necessarily implement this one yet. Requires weird sorting
	}

	if (selection_approach == GA_SELECT_ROULETTE)
	{
		/*------------------------------
		* Find the sum of current fitness values
		*------------------------------*/
		double fitSum = 0.0;
		for (uint32_t i = 0; i < algControl.populationSize; i++)
			fitSum += currentFitness[i];

		/*------------------------------
		* Spin the Roulette wheel
		*------------------------------*/
		std::random_device rd;
		std::mt19937 rng(rd());
		std::uniform_real_distribution<double> rouletteGenerator(0.0, fitSum);

		/* Run N number of times to generate N parents */
		for (uint32_t parent = 0; parent < algControl.populationSize; parent++)
		{
			/* Reset the roulette wheel */
			double rouletteSum = 0.0;
			double rouletteThreshold = rouletteGenerator(rng);

			/* Add up fitness values until the threshold value is exceeded */
			for (uint32_t i = 0; i < algControl.populationSize; i++)
			{
				rouletteSum += currentFitness[i];

				/* Assign winner to next parent slot */
				if (rouletteSum >= rouletteThreshold)
				{
					breedSelections[parent][X] = coordinate2Chromosome(currentPopulation[i][X], X);
					breedSelections[parent][Y] = coordinate2Chromosome(currentPopulation[i][Y], Y);
				}
			}
		}
	}

	if (selection_approach == GA_SELECT_TOURNAMENT)
	{
		/* NOTE:
		1. Choose N number of chromosomes at random. Have them "compete" for the highest rank.
		2. Parent1 is winner of (1)
		3. Repeat (1-2) for finding parent 2, 3, 4...etc.
		*/

		/*------------------------------
		* Create RNG
		*------------------------------*/
		std::random_device rd;
		std::mt19937 rng(rd());

		/* Ensure there are at least two competitors per tourney*/
		std::uniform_int_distribution<uint32_t> uniform_num_competitors(2, algControl.populationSize - 1);

		/* Allow any members in the current population to be chosen without bias to compete */
		std::uniform_int_distribution<uint32_t> uniform_competitor(0, algControl.populationSize - 1);

		/*------------------------------
		* Run many tournaments to determine breedSelection parents
		*------------------------------*/
		for (uint32_t i = 0; i < algControl.populationSize; i++)
		{
			/*------------------------------
			* Select competitors
			*------------------------------*/
			size_t num_competitors = uniform_num_competitors(rng);
			std::vector<uint32_t> competitors(num_competitors);

			for (uint32_t x = 0; x < num_competitors; x++)
				competitors[x] = uniform_competitor(rng);

			/*------------------------------
			* Evaluate the tournament
			*------------------------------*/
			double winner_fitness = DBL_MAX;
			uint32_t winner_idx;

			for (uint32_t x = 0; x < num_competitors; x++)
			{
				uint32_t fitIdx = competitors[x];
				double fitVal = currentFitness[fitIdx];

				if (fitVal < winner_fitness)
				{
					winner_idx = fitIdx;
					winner_fitness = fitVal;
				}
			}

			/*------------------------------
			* Assign the winner as the next parent
			*------------------------------*/
			breedSelections[i][X] = coordinate2Chromosome(currentPopulation[winner_idx][X], X);
			breedSelections[i][Y] = coordinate2Chromosome(currentPopulation[winner_idx][Y], Y);
		}
	}

	return algStatus;
}

GA_StatusTypeDef GAClass_2D::breedGeneration(GA_BreedStyle_TypeDef breeding_approach)
{
	/* Assumes two parents, simple crossover, two children */
	if (breeding_approach == GA_BREED_SIMPLE_CROSSOVER)
	{
		uint16_t upper_mask = 0xFF00;
		uint16_t lower_mask = 0x00FF;
		uint16_t *parent1, *parent2, *child1, *child2;

		//TODO: Check the ranging of this loop. Might accidentally exclude values.
		for (uint32_t i = 0; i < algControl.populationSize - 1; i += 2)
		{
			/* Grab the arrays to mate */
			parent1 = breedSelections[i];
			parent2 = breedSelections[i + 1];
			child1 = parent1;
			child2 = parent2;

			/* Handle the X-Coordinate breeding */
			child1[X] = (parent1[X] & upper_mask) | (parent2[X] & lower_mask);
			child2[X] = (parent1[X] & lower_mask) | (parent2[X] & upper_mask);

			/* Handle the Y-Coordinate breeding */
			child1[Y] = (parent1[Y] & upper_mask) | (parent2[Y] & lower_mask);
			child2[Y] = (parent1[Y] & lower_mask) | (parent2[Y] & upper_mask);

			/* Write over the original data */
			breedSelections[i] = child1;
			breedSelections[i + 1] = child2;
		}
	}

	return algStatus;
}

GA_StatusTypeDef GAClass_2D::mutateChromosomes(GA_MutateDistribution_TypeDef mutation_approach)
{
	/* Check each member in the population for mutations */
	for (uint32_t i = 0; i < algControl.populationSize; i++)
	{
		if (mutationProbability(mutation_approach) > algControl.mutation_rate)
		{
			uint16_t severity_mask = (1u << algControl.mutation_severity);

			/* Check current state of bit to be flipped:
			If lo -> flip hi.
			If hi -> flip lo. */
			if (breedSelections[i][X] & severity_mask)
				breedSelections[i][X] &= ~(severity_mask);
			else
				breedSelections[i][X] |= severity_mask;

			if (breedSelections[i][Y] & severity_mask)
				breedSelections[i][Y] &= ~(severity_mask);
			else
				breedSelections[i][Y] |= severity_mask;
		}

		/* Map the resulting data, mutated or not, back into a format usable for the fitness function */
		double newX = chromosome2Coordinate(breedSelections[i][X], X);
		double newY = chromosome2Coordinate(breedSelections[i][Y], Y);

		/* Slow check if overrunning bounds */
		if (newX < algConstraints.x_lower)
			newX = algConstraints.x_lower;

		if (newX > algConstraints.x_upper)
			newX = algConstraints.x_upper;

		if (newY < algConstraints.y_lower)
			newY = algConstraints.y_lower;

		if (newY > algConstraints.y_upper)
			newY = algConstraints.y_upper;

		currentPopulation[i][X] = newX;
		currentPopulation[i][Y] = newY;
	}

	return algStatus;
}

uint16_t GAClass_2D::coordinate2Chromosome(double coordinate, uint32_t axis)
{
	/* Algorithm from Motion Imagery Standards Board: MISB ST 1201.1 */
	switch (axis)
	{
	case X:
		return (uint16_t)(trunc(x_sF*(coordinate - x_lo) + x_offset));
		break;

	case Y:
		return (uint16_t)(trunc(y_sF*(coordinate - y_lo) + y_offset));
		break;

	default:
		return -0;
		break;
	}
}

double GAClass_2D::chromosome2Coordinate(uint16_t chromosome, uint32_t axis)
{
	/* Algorithm from Motion Imagery Standards Board: MISB ST 1201.1*/
	switch (axis)
	{
	case X:
		return (x_sR*(((float)chromosome) - x_offset) + x_lo);
		break;

	case Y:
		return (y_sR*(((float)chromosome) - y_offset) + y_lo);
		break;

	default:
		return -0.0;
		break;
	}
}

double GAClass_2D::uniformRandomNumber(double lower_bound, double upper_bound)
{
	/* Create a new RNG using the Mersenne Twister engine. A new object must
	be created each time to guarantee a good level of non-repeating randomness.*/
	std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_real_distribution<double> uniform_number(lower_bound, upper_bound);
	return uniform_number(rng);
}

double GAClass_2D::mutationProbability(GA_MutateDistribution_TypeDef distribution)
{
	switch (distribution)
	{
	case GA_MUTATE_EXPONENTIAL:
		/* Yields a probability equal to 1.0 when x is 0.0 and nearly zero when x is 5.0 */
		return exp(-1 * uniformRandomNumber(0.0, 10.0));
		break;

	default: return 0.0;
		break;
	}
}