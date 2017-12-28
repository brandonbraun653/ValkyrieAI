#include "ga_mop_steps_mutation.h"

///////////////////////////////////////////////////
/* CLASS:  MutateProbGenerator */
///////////////////////////////////////////////////
/*-----------------------------------------------
* Constructors/Destructor
*-----------------------------------------------*/
MutateProbGenerator::MutateProbGenerator(GA_MutateProbabilityMethod prob_type)
{
	probType = prob_type;
}

MutateProbGenerator::~MutateProbGenerator()
{
}

/*-----------------------------------------------
* Public Functions
*-----------------------------------------------*/
double MutateProbGenerator::get()
{
	if (probType == GA_MUTATE_PROBABILITY_EXPONENTIAL)
	{
		return exp(-1 * uniformRandomNumber(0.0, 10.0));
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
BitFlipMutator::BitFlipMutator(GA_RunMode execution_type, GA_MutateProbabilityMethod mutationProb_type)
{
	executionType = execution_type;
	mutationProbType = mutationProb_type;
}

BitFlipMutator::~BitFlipMutator()
{
}

/*-----------------------------------------------
* Public Functions
*-----------------------------------------------*/
void BitFlipMutator::mutate(hPID_Chromosomes* in_bredChrom, GA_ConverganceCriteria_sPtr in_convgCriteria, PID_ControlGoals_sPtr in_config,
	hPID_DataVector* out_pidVals, mapCoeff_t* in_kp, mapCoeff_t* in_ki, mapCoeff_t* in_kd)
{
	cm_data = out_pidVals;
	cm_chrom = in_bredChrom;
	cm_convergence = in_convgCriteria;
	cm_config = in_config;

	mapCM_Kp = in_kp;
	mapCM_Ki = in_ki;
	mapCM_Kd = in_kd;

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
void BitFlipMutator::calculate_cpu_single_threaded()
{
	mutateKp();
	mutateKi();
	mutateKd();
}

void BitFlipMutator::calculate_cpu_multi_threaded()
{
	boost::thread_group tgroup;

	tgroup.create_thread(boost::bind(&BitFlipMutator::mutateKp, this));
	tgroup.create_thread(boost::bind(&BitFlipMutator::mutateKi, this));
	tgroup.create_thread(boost::bind(&BitFlipMutator::mutateKd, this));

	tgroup.join_all();
}

void BitFlipMutator::calculate_gpu_single_threaded()
{
	std::cout << "GPU Single Threaded Mode Currently Not Supported." << std::endl;
}

void BitFlipMutator::calculate_gpu_multi_threaded()
{
	std::cout << "GPU Multi Threaded Mode Currently Not Supported." << std::endl;
}

void BitFlipMutator::mutateKp()
{
	MutateProbGenerator mutationProbability(mutationProbType);
	uint16_t severity_mask = (1u << cm_convergence->mutation_severity);

	/* Check each population member for possibility to mutate */
	for (int member = 0; member < cm_chrom->Kp.size(); member++)
	{
		if (mutationProbability.get() > cm_convergence->mutation_threshold)
		{
			/* Check current state of bit to be flipped:
			If lo -> flip hi.
			If hi -> flip lo. */

			/* Kp*/
			if (cm_chrom->Kp.data()[member] & severity_mask)
				cm_chrom->Kp.data()[member] &= ~(severity_mask);
			else
				cm_chrom->Kp.data()[member] |= severity_mask;
		}

		/* Map the resulting data, mutated or not, back into a usable format */
		cm_data->Kp.data()[member] = chromosome2Data(mapCM_Kp, cm_chrom->Kp.data()[member]);
		

		/* Force boundaries */
		if (cm_data->Kp.data()[member] > cm_config->pid_limits.Kp_limits_upper)
			cm_data->Kp.data()[member] = cm_config->pid_limits.Kp_limits_upper;
		if (cm_data->Kp.data()[member] < cm_config->pid_limits.Kp_limits_lower) 
			cm_data->Kp.data()[member] = cm_config->pid_limits.Kp_limits_lower; 
	}
}

void BitFlipMutator::mutateKi()
{
	MutateProbGenerator mutationProbability(mutationProbType);
	uint16_t severity_mask = (1u << cm_convergence->mutation_severity);

	/* Check each population member for possibility to mutate */
	for (int member = 0; member < cm_chrom->Ki.size(); member++)
	{
		if (mutationProbability.get() > cm_convergence->mutation_threshold)
		{
			/* Check current state of bit to be flipped:
			If lo -> flip hi.
			If hi -> flip lo. */

			/* Ki*/
			if (cm_chrom->Ki.data()[member] & severity_mask)
				cm_chrom->Ki.data()[member] &= ~(severity_mask);
			else
				cm_chrom->Ki.data()[member] |= severity_mask;
		}

		/* Map the resulting data, mutated or not, back into a usable format */
		cm_data->Ki.data()[member] = chromosome2Data(mapCM_Ki, cm_chrom->Ki.data()[member]);


		/* Force boundaries */
		if (cm_data->Ki.data()[member] > cm_config->pid_limits.Ki_limits_upper)
			cm_data->Ki.data()[member] = cm_config->pid_limits.Ki_limits_upper;
		if (cm_data->Ki.data()[member] < cm_config->pid_limits.Ki_limits_lower)
			cm_data->Ki.data()[member] = cm_config->pid_limits.Ki_limits_lower;
	}
}

void BitFlipMutator::mutateKd()
{
	MutateProbGenerator mutationProbability(mutationProbType);
	uint16_t severity_mask = (1u << cm_convergence->mutation_severity);

	/* Check each population member for possibility to mutate */
	for (int member = 0; member < cm_chrom->Kd.size(); member++)
	{
		if (mutationProbability.get() > cm_convergence->mutation_threshold)
		{
			/* Check current state of bit to be flipped:
			If lo -> flip hi.
			If hi -> flip lo. */

			/* Kd*/
			if (cm_chrom->Kd.data()[member] & severity_mask)
				cm_chrom->Kd.data()[member] &= ~(severity_mask);
			else
				cm_chrom->Kd.data()[member] |= severity_mask;
		}

		/* Map the resulting data, mutated or not, back into a usable format */
		cm_data->Kd.data()[member] = chromosome2Data(mapCM_Kd, cm_chrom->Kd.data()[member]);


		/* Force boundaries */
		if (cm_data->Kd.data()[member] > cm_config->pid_limits.Kd_limits_upper)
			cm_data->Kd.data()[member] = cm_config->pid_limits.Kd_limits_upper;
		if (cm_data->Kd.data()[member] < cm_config->pid_limits.Kd_limits_lower)
			cm_data->Kd.data()[member] = cm_config->pid_limits.Kd_limits_lower;
	}
}

///////////////////////////////////////////////////
/* CLASS:  AddSubMutator */
///////////////////////////////////////////////////
/*-----------------------------------------------
* Constructors/Destructor
*-----------------------------------------------*/
AddSubMutator::AddSubMutator(GA_RunMode execution_type, GA_MutateProbabilityMethod mutationProb_type, GA_Resolution res_type)
{
	executionType = execution_type;
	mutationProbType = mutationProb_type;
	resolutionType = res_type;
}

AddSubMutator::~AddSubMutator()
{
}

/*-----------------------------------------------
* Public Functions
*-----------------------------------------------*/
void AddSubMutator::mutate(hPID_Chromosomes* in_bredChrom, GA_ConverganceCriteria_sPtr in_convgCriteria, PID_ControlGoals_sPtr in_config,
	hPID_DataVector* out_pidVals, mapCoeff_t* in_kp, mapCoeff_t* in_ki, mapCoeff_t* in_kd)
{
	as_data = out_pidVals;
	as_chrom = in_bredChrom;
	as_convergence = in_convgCriteria;
	as_config = in_config;

	mapAS_Kp = in_kp;
	mapAS_Ki = in_ki;
	mapAS_Kd = in_kd;

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
void AddSubMutator::calculate_cpu_single_threaded()
{
	mutateKp();
	mutateKi();
	mutateKd();
}

void AddSubMutator::calculate_cpu_multi_threaded()
{
	boost::thread_group tgroup;

	tgroup.create_thread(boost::bind(&AddSubMutator::mutateKp, this));
	tgroup.create_thread(boost::bind(&AddSubMutator::mutateKi, this));
	tgroup.create_thread(boost::bind(&AddSubMutator::mutateKd, this));

	tgroup.join_all();
}

void AddSubMutator::calculate_gpu_single_threaded()
{
	std::cout << "GPU Single Threaded Mode Currently Not Supported." << std::endl;
}

void AddSubMutator::calculate_gpu_multi_threaded()
{
	std::cout << "GPU Multi Threaded Mode Currently Not Supported." << std::endl;
}

void AddSubMutator::mutateKp()
{
	MutateProbGenerator mutationProbability(mutationProbType);

	/* Instantiate the RNG device */
	std::random_device rd;
	std::mt19937 rng(rd());

	double maxPct = 0.25;

	/* Configure the range of mutation severity */
	std::uniform_real_distribution<double> kp_Adjustment(
		-maxPct*as_config->pid_limits.Kp_limits_upper,
		maxPct*as_config->pid_limits.Kp_limits_upper);


	for (int member = 0; member < as_chrom->Kp.size(); member++)
	{
		/* Map the resulting data, mutated or not, back into a usable format */
		as_data->Kp.data()[member] = chromosome2Data(mapAS_Kp, as_chrom->Kp.data()[member]);

		/* Check for possibility to mutate */
		if (mutationProbability.get() > as_convergence->mutation_threshold)
			as_data->Kp.data()[member] += kp_Adjustment(rng);

		/* Force boundaries */
		if (as_data->Kp.data()[member] > as_config->pid_limits.Kp_limits_upper) //Greater than max
			as_data->Kp.data()[member] = as_config->pid_limits.Kp_limits_upper; //	Set as max
		if (as_data->Kp.data()[member] < as_config->pid_limits.Kp_limits_lower) //Less than min
			as_data->Kp.data()[member] = as_config->pid_limits.Kp_limits_lower; //  Set as min
	}
}

void AddSubMutator::mutateKi()
{
	MutateProbGenerator mutationProbability(mutationProbType);

	/* Instantiate the RNG device */
	std::random_device rd;
	std::mt19937 rng(rd());

	double maxPct = 0.25;

	/* Configure the range of mutation severity */
	std::uniform_real_distribution<double> ki_Adjustment(
		-maxPct*as_config->pid_limits.Ki_limits_upper,
		maxPct*as_config->pid_limits.Ki_limits_upper);

	for (int member = 0; member < as_chrom->Ki.size(); member++)
	{
		/* Map the resulting data, mutated or not, back into a usable format */
		as_data->Ki.data()[member] = chromosome2Data(mapAS_Ki, as_chrom->Ki.data()[member]);

		/* Check for possibility to mutate */
		if (mutationProbability.get() > as_convergence->mutation_threshold)
			as_data->Ki.data()[member] += ki_Adjustment(rng);

		/* Force boundaries */
		if (as_data->Ki.data()[member] > as_config->pid_limits.Ki_limits_upper)
			as_data->Ki.data()[member] = as_config->pid_limits.Ki_limits_upper;
		if (as_data->Ki.data()[member] < as_config->pid_limits.Ki_limits_lower)
			as_data->Ki.data()[member] = as_config->pid_limits.Ki_limits_lower;
	}
}

void AddSubMutator::mutateKd()
{
	MutateProbGenerator mutationProbability(mutationProbType);

	/* Instantiate the RNG device */
	std::random_device rd;
	std::mt19937 rng(rd());

	double maxPct = 0.25;

	/* Configure the range of mutation severity */
	std::uniform_real_distribution<double> kd_Adjustment(
		-maxPct*as_config->pid_limits.Kd_limits_upper,
		maxPct*as_config->pid_limits.Kd_limits_upper);


	for (int member = 0; member < as_chrom->Kd.size(); member++)
	{
		/* Map the resulting data, mutated or not, back into a usable format */
		as_data->Kd.data()[member] = chromosome2Data(mapAS_Kd, as_chrom->Kd.data()[member]);

		/* Check for possibility to mutate */
		if (mutationProbability.get() > as_convergence->mutation_threshold)
			as_data->Kd.data()[member] += kd_Adjustment(rng);

		/* Force boundaries */
		if (as_data->Kd.data()[member] > as_config->pid_limits.Kd_limits_upper)
			as_data->Kd.data()[member] = as_config->pid_limits.Kd_limits_upper;
		if (as_data->Kd.data()[member] < as_config->pid_limits.Kd_limits_lower)
			as_data->Kd.data()[member] = as_config->pid_limits.Kd_limits_lower;

	}
}