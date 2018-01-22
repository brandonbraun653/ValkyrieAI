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
/*-----------------------------------------------
* Constructors/Destructor
*-----------------------------------------------*/
BitFlipMutator::BitFlipMutator()
{
	auto dist = boost::random::uniform_int_distribution<>(0, maxSeverity);
	severityRNG = boost::make_shared <RNGInstance<boost::mt19937, boost::random::uniform_int_distribution<>>>(dist);
	
	mutateRNG = boost::make_shared<MutateProbGenerator>();
}

BitFlipMutator::~BitFlipMutator()
{
}

/*-----------------------------------------------
* Public Functions
*-----------------------------------------------*/
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

/*-----------------------------------------------
* Private Functions
*-----------------------------------------------*/
void BitFlipMutator::mutateKp()
{

}

void BitFlipMutator::mutateKi()
{
	
}

void BitFlipMutator::mutateKd()
{

}


///////////////////////////////////////////////////
/* CLASS:  AddSubMutator */
///////////////////////////////////////////////////
/*-----------------------------------------------
* Constructors/Destructor
*-----------------------------------------------*/
AddSubMutator::AddSubMutator()
{
}

AddSubMutator::~AddSubMutator()
{
}

/*-----------------------------------------------
* Public Functions
*-----------------------------------------------*/
void AddSubMutator::mutate(const GA_MutateDataInput input, GA_MutateDataOutput& output)
{

}

/*-----------------------------------------------
* Private Functions
*-----------------------------------------------*/
void AddSubMutator::mutateKp()
{

}

void AddSubMutator::mutateKi()
{

}

void AddSubMutator::mutateKd()
{

}

// ///////////////////////////////////////////////////
// /* CLASS:  MutateProbGenerator */
// ///////////////////////////////////////////////////
// /*-----------------------------------------------
// * Constructors/Destructor
// *-----------------------------------------------*/
// MutateProbGenerator::MutateProbGenerator(GA_MutateProbabilityMethod prob_type)
// {
// 	probType = prob_type;
// }
// 
// MutateProbGenerator::~MutateProbGenerator()
// {
// }
// 
// /*-----------------------------------------------
// * Public Functions
// *-----------------------------------------------*/
// double MutateProbGenerator::get()
// {
// 	if (probType == GA_MUTATE_PROBABILITY_EXPONENTIAL)
// 	{
// 		return exp(-1 * uniformRandomNumber(0.0, 10.0));
// 	}
// 
// 	if (probType == GA_MUTATE_PROBABILITY_CHI_SQUARED)
// 	{
// 		std::cout << "Selected mutation method currently not supported." << std::endl;
// 		return 0.0;
// 	}
// 
// 	if (probType == GA_MUTATE_PROBABILITY_GAMMA)
// 	{
// 		std::cout << "Selected mutation method currently not supported." << std::endl;
// 		return 0.0;
// 	}
// 
// 	if (probType == GA_MUTATE_PROBABILITY_POISSON)
// 	{
// 		std::cout << "Selected mutation method currently not supported." << std::endl;
// 		return 0.0;
// 	}
// 
// 	if (probType == GA_MUTATE_PROBABILITY_WEIBULL)
// 	{
// 		std::cout << "Selected mutation method currently not supported." << std::endl;
// 		return 0.0;
// 	}
// 	else
// 		return 0.0;
// }
// 
// ///////////////////////////////////////////////////
// /* CLASS:  BitFlipMutator */
// ///////////////////////////////////////////////////
// /*-----------------------------------------------
// * Constructors/Destructor
// *-----------------------------------------------*/
// BitFlipMutator::BitFlipMutator(GA_RunMode execution_type, GA_MutateProbabilityMethod mutationProb_type)
// {
// 	executionType = execution_type;
// 	mutationProbType = mutationProb_type;
// }
// 
// BitFlipMutator::~BitFlipMutator()
// {
// }
// 
// /*-----------------------------------------------
// * Public Functions
// *-----------------------------------------------*/
// void BitFlipMutator::mutate(hPID_Chromosomes* in_bredChrom, FCSOptimizer_AdvConstraints_sPtr in_convgCriteria, PID_ControlSettings_sPtr in_config,
// 	hPID_DataVector* out_pidVals, FCSOptimizer_MappingCoeff* in_kp, FCSOptimizer_MappingCoeff* in_ki, FCSOptimizer_MappingCoeff* in_kd)
// {
// 	cm_data = out_pidVals;
// 	cm_chrom = in_bredChrom;
// 	cm_convergence = in_convgCriteria;
// 	cm_config = in_config;
// 
// 	mapCM_Kp = in_kp;
// 	mapCM_Ki = in_ki;
// 	mapCM_Kd = in_kd;
// 
// 	if (executionType == SINGLE_THREADED)
// 		calculate_cpu_single_threaded();
// 
// 	if (executionType == MULTI_THREADED)
// 		calculate_cpu_multi_threaded();
// 
// 	if (executionType == SINGLE_THREADED_WITH_CUDA)
// 		calculate_gpu_single_threaded();
// 
// 	if (executionType == MULTI_THREADED_WITH_CUDA)
// 		calculate_gpu_multi_threaded();
// }
// 
// /*-----------------------------------------------
// * Private Functions
// *-----------------------------------------------*/
// void BitFlipMutator::calculate_cpu_single_threaded()
// {
// 	mutateKp();
// 	mutateKi();
// 	mutateKd();
// }
// 
// void BitFlipMutator::calculate_cpu_multi_threaded()
// {
// 	boost::thread_group tgroup;
// 
// 	tgroup.create_thread(boost::bind(&BitFlipMutator::mutateKp, this));
// 	tgroup.create_thread(boost::bind(&BitFlipMutator::mutateKi, this));
// 	tgroup.create_thread(boost::bind(&BitFlipMutator::mutateKd, this));
// 
// 	tgroup.join_all();
// }
// 
// void BitFlipMutator::calculate_gpu_single_threaded()
// {
// 	std::cout << "GPU Single Threaded Mode Currently Not Supported." << std::endl;
// }
// 
// void BitFlipMutator::calculate_gpu_multi_threaded()
// {
// 	std::cout << "GPU Multi Threaded Mode Currently Not Supported." << std::endl;
// }
// 
// void BitFlipMutator::mutateKp()
// {
// 	MutateProbGenerator mutationProbability(mutationProbType);
// 	uint16_t severity_mask = (1u << cm_convergence->mutation_severity);
// 
// 	/* Check each population member for possibility to mutate */
// 	for (int member = 0; member < cm_chrom->Kp.size(); member++)
// 	{
// 		if (mutationProbability.get() > cm_convergence->mutation_threshold)
// 		{
// 			/* Check current state of bit to be flipped:
// 			If lo -> flip hi.
// 			If hi -> flip lo. */
// 
// 			/* Kp*/
// 			if (cm_chrom->Kp.data()[member] & severity_mask)
// 				cm_chrom->Kp.data()[member] &= ~(severity_mask);
// 			else
// 				cm_chrom->Kp.data()[member] |= severity_mask;
// 		}
// 
// 		/* Map the resulting data, mutated or not, back into a usable format */
// 		cm_data->Kp.data()[member] = chromosome2Data(mapCM_Kp, cm_chrom->Kp.data()[member]);
// 		
// 
// 		/* Force boundaries */
// 		if (cm_data->Kp.data()[member] > cm_config->tuningLimits.Kp_limits_upper)
// 			cm_data->Kp.data()[member] = cm_config->tuningLimits.Kp_limits_upper;
// 		if (cm_data->Kp.data()[member] < cm_config->tuningLimits.Kp_limits_lower) 
// 			cm_data->Kp.data()[member] = cm_config->tuningLimits.Kp_limits_lower; 
// 	}
// }
// 
// void BitFlipMutator::mutateKi()
// {
// 	MutateProbGenerator mutationProbability(mutationProbType);
// 	uint16_t severity_mask = (1u << cm_convergence->mutation_severity);
// 
// 	/* Check each population member for possibility to mutate */
// 	for (int member = 0; member < cm_chrom->Ki.size(); member++)
// 	{
// 		if (mutationProbability.get() > cm_convergence->mutation_threshold)
// 		{
// 			/* Check current state of bit to be flipped:
// 			If lo -> flip hi.
// 			If hi -> flip lo. */
// 
// 			/* Ki*/
// 			if (cm_chrom->Ki.data()[member] & severity_mask)
// 				cm_chrom->Ki.data()[member] &= ~(severity_mask);
// 			else
// 				cm_chrom->Ki.data()[member] |= severity_mask;
// 		}
// 
// 		/* Map the resulting data, mutated or not, back into a usable format */
// 		cm_data->Ki.data()[member] = chromosome2Data(mapCM_Ki, cm_chrom->Ki.data()[member]);
// 
// 
// 		/* Force boundaries */
// 		if (cm_data->Ki.data()[member] > cm_config->tuningLimits.Ki_limits_upper)
// 			cm_data->Ki.data()[member] = cm_config->tuningLimits.Ki_limits_upper;
// 		if (cm_data->Ki.data()[member] < cm_config->tuningLimits.Ki_limits_lower)
// 			cm_data->Ki.data()[member] = cm_config->tuningLimits.Ki_limits_lower;
// 	}
// }
// 
// void BitFlipMutator::mutateKd()
// {
// 	MutateProbGenerator mutationProbability(mutationProbType);
// 	uint16_t severity_mask = (1u << cm_convergence->mutation_severity);
// 
// 	/* Check each population member for possibility to mutate */
// 	for (int member = 0; member < cm_chrom->Kd.size(); member++)
// 	{
// 		if (mutationProbability.get() > cm_convergence->mutation_threshold)
// 		{
// 			/* Check current state of bit to be flipped:
// 			If lo -> flip hi.
// 			If hi -> flip lo. */
// 
// 			/* Kd*/
// 			if (cm_chrom->Kd.data()[member] & severity_mask)
// 				cm_chrom->Kd.data()[member] &= ~(severity_mask);
// 			else
// 				cm_chrom->Kd.data()[member] |= severity_mask;
// 		}
// 
// 		/* Map the resulting data, mutated or not, back into a usable format */
// 		cm_data->Kd.data()[member] = chromosome2Data(mapCM_Kd, cm_chrom->Kd.data()[member]);
// 
// 
// 		/* Force boundaries */
// 		if (cm_data->Kd.data()[member] > cm_config->tuningLimits.Kd_limits_upper)
// 			cm_data->Kd.data()[member] = cm_config->tuningLimits.Kd_limits_upper;
// 		if (cm_data->Kd.data()[member] < cm_config->tuningLimits.Kd_limits_lower)
// 			cm_data->Kd.data()[member] = cm_config->tuningLimits.Kd_limits_lower;
// 	}
// }
// 
// ///////////////////////////////////////////////////
// /* CLASS:  AddSubMutator */
// ///////////////////////////////////////////////////
// /*-----------------------------------------------
// * Constructors/Destructor
// *-----------------------------------------------*/
// AddSubMutator::AddSubMutator(GA_RunMode execution_type, GA_MutateProbabilityMethod mutationProb_type, GA_Resolution res_type)
// {
// 	executionType = execution_type;
// 	mutationProbType = mutationProb_type;
// 	resolutionType = res_type;
// }
// 
// AddSubMutator::~AddSubMutator()
// {
// }
// 
// /*-----------------------------------------------
// * Public Functions
// *-----------------------------------------------*/
// void AddSubMutator::mutate(hPID_Chromosomes* in_bredChrom, FCSOptimizer_AdvConstraints_sPtr in_convgCriteria, PID_ControlSettings_sPtr in_config,
// 	hPID_DataVector* out_pidVals, FCSOptimizer_MappingCoeff* in_kp, FCSOptimizer_MappingCoeff* in_ki, FCSOptimizer_MappingCoeff* in_kd)
// {
// 	as_data = out_pidVals;
// 	as_chrom = in_bredChrom;
// 	as_convergence = in_convgCriteria;
// 	as_config = in_config;
// 
// 	mapAS_Kp = in_kp;
// 	mapAS_Ki = in_ki;
// 	mapAS_Kd = in_kd;
// 
// 	if (executionType == SINGLE_THREADED)
// 		calculate_cpu_single_threaded();
// 
// 	if (executionType == MULTI_THREADED)
// 		calculate_cpu_multi_threaded();
// 
// 	if (executionType == SINGLE_THREADED_WITH_CUDA)
// 		calculate_gpu_single_threaded();
// 
// 	if (executionType == MULTI_THREADED_WITH_CUDA)
// 		calculate_gpu_multi_threaded();
// }
// 
// /*-----------------------------------------------
// * Private Functions
// *-----------------------------------------------*/
// void AddSubMutator::calculate_cpu_single_threaded()
// {
// 	mutateKp();
// 	mutateKi();
// 	mutateKd();
// }
// 
// void AddSubMutator::calculate_cpu_multi_threaded()
// {
// 	boost::thread_group tgroup;
// 
// 	tgroup.create_thread(boost::bind(&AddSubMutator::mutateKp, this));
// 	tgroup.create_thread(boost::bind(&AddSubMutator::mutateKi, this));
// 	tgroup.create_thread(boost::bind(&AddSubMutator::mutateKd, this));
// 
// 	tgroup.join_all();
// }
// 
// void AddSubMutator::calculate_gpu_single_threaded()
// {
// 	std::cout << "GPU Single Threaded Mode Currently Not Supported." << std::endl;
// }
// 
// void AddSubMutator::calculate_gpu_multi_threaded()
// {
// 	std::cout << "GPU Multi Threaded Mode Currently Not Supported." << std::endl;
// }
// 
// void AddSubMutator::mutateKp()
// {
// 	MutateProbGenerator mutationProbability(mutationProbType);
// 
// 	/* Instantiate the RNG device */
// 	std::random_device rd;
// 	std::mt19937 rng(rd());
// 
// 	double maxPct = 0.25;
// 
// 	/* Configure the range of mutation severity */
// 	std::uniform_real_distribution<double> kp_Adjustment(
// 		-maxPct*as_config->tuningLimits.Kp_limits_upper,
// 		maxPct*as_config->tuningLimits.Kp_limits_upper);
// 
// 
// 	for (int member = 0; member < as_chrom->Kp.size(); member++)
// 	{
// 		/* Map the resulting data, mutated or not, back into a usable format */
// 		as_data->Kp.data()[member] = chromosome2Data(mapAS_Kp, as_chrom->Kp.data()[member]);
// 
// 		/* Check for possibility to mutate */
// 		if (mutationProbability.get() > as_convergence->mutation_threshold)
// 			as_data->Kp.data()[member] += kp_Adjustment(rng);
// 
// 		/* Force boundaries */
// 		if (as_data->Kp.data()[member] > as_config->tuningLimits.Kp_limits_upper) //Greater than max
// 			as_data->Kp.data()[member] = as_config->tuningLimits.Kp_limits_upper; //	Set as max
// 		if (as_data->Kp.data()[member] < as_config->tuningLimits.Kp_limits_lower) //Less than min
// 			as_data->Kp.data()[member] = as_config->tuningLimits.Kp_limits_lower; //  Set as min
// 	}
// }
// 
// void AddSubMutator::mutateKi()
// {
// 	MutateProbGenerator mutationProbability(mutationProbType);
// 
// 	/* Instantiate the RNG device */
// 	std::random_device rd;
// 	std::mt19937 rng(rd());
// 
// 	double maxPct = 0.25;
// 
// 	/* Configure the range of mutation severity */
// 	std::uniform_real_distribution<double> ki_Adjustment(
// 		-maxPct*as_config->tuningLimits.Ki_limits_upper,
// 		maxPct*as_config->tuningLimits.Ki_limits_upper);
// 
// 	for (int member = 0; member < as_chrom->Ki.size(); member++)
// 	{
// 		/* Map the resulting data, mutated or not, back into a usable format */
// 		as_data->Ki.data()[member] = chromosome2Data(mapAS_Ki, as_chrom->Ki.data()[member]);
// 
// 		/* Check for possibility to mutate */
// 		if (mutationProbability.get() > as_convergence->mutation_threshold)
// 			as_data->Ki.data()[member] += ki_Adjustment(rng);
// 
// 		/* Force boundaries */
// 		if (as_data->Ki.data()[member] > as_config->tuningLimits.Ki_limits_upper)
// 			as_data->Ki.data()[member] = as_config->tuningLimits.Ki_limits_upper;
// 		if (as_data->Ki.data()[member] < as_config->tuningLimits.Ki_limits_lower)
// 			as_data->Ki.data()[member] = as_config->tuningLimits.Ki_limits_lower;
// 	}
// }
// 
// void AddSubMutator::mutateKd()
// {
// 	MutateProbGenerator mutationProbability(mutationProbType);
// 
// 	/* Instantiate the RNG device */
// 	std::random_device rd;
// 	std::mt19937 rng(rd());
// 
// 	double maxPct = 0.25;
// 
// 	/* Configure the range of mutation severity */
// 	std::uniform_real_distribution<double> kd_Adjustment(
// 		-maxPct*as_config->tuningLimits.Kd_limits_upper,
// 		maxPct*as_config->tuningLimits.Kd_limits_upper);
// 
// 
// 	for (int member = 0; member < as_chrom->Kd.size(); member++)
// 	{
// 		/* Map the resulting data, mutated or not, back into a usable format */
// 		as_data->Kd.data()[member] = chromosome2Data(mapAS_Kd, as_chrom->Kd.data()[member]);
// 
// 		/* Check for possibility to mutate */
// 		if (mutationProbability.get() > as_convergence->mutation_threshold)
// 			as_data->Kd.data()[member] += kd_Adjustment(rng);
// 
// 		/* Force boundaries */
// 		if (as_data->Kd.data()[member] > as_config->tuningLimits.Kd_limits_upper)
// 			as_data->Kd.data()[member] = as_config->tuningLimits.Kd_limits_upper;
// 		if (as_data->Kd.data()[member] < as_config->tuningLimits.Kd_limits_lower)
// 			as_data->Kd.data()[member] = as_config->tuningLimits.Kd_limits_lower;
// 
// 	}
// }