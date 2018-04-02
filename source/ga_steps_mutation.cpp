#include "ga_steps_mutation.h"

///////////////////////////////////////////////////
/* CLASS:  MutateProbGenerator */
///////////////////////////////////////////////////
MutateProbGenerator::MutateProbGenerator()
{
	auto dist = boost::random::uniform_real_distribution<>(0.0, 10.0);
	rngUniform = boost::make_shared <RNGInstance<boost::mt19937, boost::random::uniform_real_distribution<>>>(dist);
}

MutateProbGenerator::~MutateProbGenerator()
{
}

double MutateProbGenerator::get(GA_METHOD_MutateProbability probType)
{
	//TODO: add boundary checking to ensure the return value will ALWAYS lie between 0.0 and 1.0

	if (probType == GA_MUTATE_PROBABILITY_EXPONENTIAL)
	{
		rngUniform->acquireEngine();
		double val = exp(-1.0 * rngUniform->getDouble());
		rngUniform->releaseEngine();

		return val;
	}

	if (probType == GA_MUTATE_PROBABILITY_CHI_SQUARED)
	{
		std::cout << "Selected mutation method currently not supported." << std::endl;
		return 0.0;
	}

	if (probType == GA_MUTATE_PROBABILITY_GAMMA)
	{
		std::cout << "Selected mutation method currently not supported." << std::endl;
		return 0.0;
	}

	if (probType == GA_MUTATE_PROBABILITY_POISSON)
	{
		std::cout << "Selected mutation method currently not supported." << std::endl;
		return 0.0;
	}

	if (probType == GA_MUTATE_PROBABILITY_WEIBULL)
	{
		std::cout << "Selected mutation method currently not supported." << std::endl;
		return 0.0;
	}
	else
		return 0.0;
}


///////////////////////////////////////////////////
/* CLASS:  BitFlipMutator */
///////////////////////////////////////////////////
BitFlipMutator::BitFlipMutator()
{
	auto dist = boost::random::uniform_int_distribution<>(0, maxSeverity);
	severityRNG = boost::make_shared <RNGInstance<boost::mt19937, boost::random::uniform_int_distribution<>>>(dist);
	mutateRNG = boost::make_shared<MutateProbGenerator>();
}

void BitFlipMutator::mutate(const GA_MutateDataInput input, GA_MutateDataOutput& output)
{
	boost::container::vector<GA_PIDChromosome<uint16_t>> localChrom;

	/*------------------------------
	* Convert to bit fields before mutation if needed
	*-------------------------------*/
	if (input.chromType == MAPPING_TYPE_REAL)
	{
		/* Ensure the output has been correctly allocated */
		localChrom.resize(input.d_chrom.size());

		for (int member = 0; member < input.d_chrom.size(); member++)
		{
			localChrom[member].Kp = data2Chromosome(input.mapCoeff_Kp, input.d_chrom[member].Kp);
			localChrom[member].Ki = data2Chromosome(input.mapCoeff_Ki, input.d_chrom[member].Ki);
			localChrom[member].Kd = data2Chromosome(input.mapCoeff_Kd, input.d_chrom[member].Kd);
		}
	}
	else if (input.chromType == MAPPING_TYPE_BIT_FIELD)
		localChrom = input.u16_chrom;

	/*------------------------------
	* Do the actual mutation 
	*-------------------------------*/
	severityRNG->acquireEngine();

	for (int member = 0; member < localChrom.size(); member++)
	{
		double mutateProb = mutateRNG->get(input.mutateProbType);

		if (mutateProb > input.mutationProbabilityThreshold)
		{
			//Toggle a random bit in the chromosome data 
			localChrom[member].Kp ^= (1 << severityRNG->getInt());
			localChrom[member].Ki ^= (1 << severityRNG->getInt());
			localChrom[member].Kd ^= (1 << severityRNG->getInt());
		}
	}

	severityRNG->releaseEngine();

	/*------------------------------
	* Convert back to the correct optimizer chromosome representation 
	*-------------------------------*/
	if (input.optimizerChromType == MAPPING_TYPE_REAL)
	{
		output.chromType = input.optimizerChromType;

		output.d_chrom.resize(localChrom.size());

		for (int member = 0; member < localChrom.size(); member++)
		{
			output.d_chrom[member].Kp = chromosome2Data(input.mapCoeff_Kp, localChrom[member].Kp);
			output.d_chrom[member].Ki = chromosome2Data(input.mapCoeff_Ki, localChrom[member].Ki);
			output.d_chrom[member].Kd = chromosome2Data(input.mapCoeff_Kd, localChrom[member].Kd);
		}
	}
	else if (input.optimizerChromType == MAPPING_TYPE_BIT_FIELD)
	{
		output.chromType = input.optimizerChromType;
		output.u16_chrom = localChrom;
	}
}


///////////////////////////////////////////////////
/* CLASS:  AddSubMutator */
///////////////////////////////////////////////////
AddSubMutator::AddSubMutator()
{
	auto dist = boost::random::uniform_int_distribution<>(-maxRange, maxRange);
	severityRNG = boost::make_shared <RNGInstance<boost::mt19937, boost::random::uniform_int_distribution<>>>(dist);
	mutateRNG = boost::make_shared<MutateProbGenerator>();
}

void AddSubMutator::mutate(const GA_MutateDataInput input, GA_MutateDataOutput& output)
{
	boost::container::vector<GA_PIDChromosome<double>> localChrom;

	/*------------------------------
	* Convert from bit fields before mutation if needed
	*-------------------------------*/
	if (input.chromType == MAPPING_TYPE_BIT_FIELD)
	{
		/* Ensure the output has been correctly allocated */
		localChrom.resize(input.u16_chrom.size());

		for (int member = 0; member < input.u16_chrom.size(); member++)
		{
			localChrom[member].Kp = chromosome2Data(input.mapCoeff_Kp, input.u16_chrom[member].Kp);
			localChrom[member].Ki = chromosome2Data(input.mapCoeff_Ki, input.u16_chrom[member].Ki);
			localChrom[member].Kd = chromosome2Data(input.mapCoeff_Kd, input.u16_chrom[member].Kd);
		}
	}
	else if (input.chromType == MAPPING_TYPE_REAL)
	{
		localChrom = input.d_chrom;
	}

	/*------------------------------
	* Do the actual mutation
	*-------------------------------*/
	severityRNG->acquireEngine();

	for (int member = 0; member < localChrom.size(); member++)
	{
		double mutateProb = mutateRNG->get(input.mutateProbType);

		if (mutateProb > input.mutationProbabilityThreshold)
		{
			localChrom[member].Kp += severityRNG->getInt();
			localChrom[member].Ki += severityRNG->getInt();
			localChrom[member].Kd += severityRNG->getInt();
		}
	}

	severityRNG->releaseEngine();

	/*------------------------------
	* Convert back to the correct optimizer chromosome representation
	*-------------------------------*/
	if (input.optimizerChromType == MAPPING_TYPE_BIT_FIELD)
	{
		output.chromType = input.optimizerChromType;

		output.u16_chrom.resize(localChrom.size());

		for (int member = 0; member < localChrom.size(); member++)
		{
			output.u16_chrom[member].Kp = data2Chromosome(input.mapCoeff_Kp, localChrom[member].Kp);
			output.u16_chrom[member].Ki = data2Chromosome(input.mapCoeff_Ki, localChrom[member].Ki);
			output.u16_chrom[member].Kd = data2Chromosome(input.mapCoeff_Kd, localChrom[member].Kd);
		}
	}
	else if (input.optimizerChromType == MAPPING_TYPE_REAL)
	{
		output.chromType = input.optimizerChromType;
		output.d_chrom = localChrom;
	}
}