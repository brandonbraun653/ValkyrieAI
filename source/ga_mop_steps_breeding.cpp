#include "ga_mop_steps_breeding.h"

///////////////////////////////////////////////////
/* CLASS:  SimpleCrossover */
///////////////////////////////////////////////////
/*-----------------------------------------------
* Constructors/Destructor
*-----------------------------------------------*/
SimpleCrossover::SimpleCrossover(GA_RunMode execution_type)
{
	executionType = execution_type;
}

SimpleCrossover::~SimpleCrossover()
{
}

/*-----------------------------------------------
* Public Functions
*-----------------------------------------------*/
void SimpleCrossover::breed(iVec in_parents, hPID_DataVector* in_pidVals, hPID_Chromosomes* out_breed, mapCoeff_t* in_kp, mapCoeff_t* in_ki, mapCoeff_t* in_kd)
{
	sc_parents = in_parents;
	sc_data = in_pidVals;
	sc_chrom = out_breed;

	mapCF_Kp = in_kp;
	mapCF_Ki = in_ki;
	mapCF_Kd = in_kd;

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
void SimpleCrossover::calculate_cpu_single_threaded()
{
	breedKp();
	breedKi();
	breedKd();
}

void SimpleCrossover::calculate_cpu_multi_threaded()
{
	boost::thread_group tgroup;

	tgroup.create_thread(boost::bind(&SimpleCrossover::breedKp, this));
	tgroup.create_thread(boost::bind(&SimpleCrossover::breedKi, this));
	tgroup.create_thread(boost::bind(&SimpleCrossover::breedKd, this));

	tgroup.join_all();
}

void SimpleCrossover::calculate_gpu_single_threaded()
{
	std::cout << "GPU Single Threaded Mode Currently Not Supported." << std::endl;
}

void SimpleCrossover::calculate_gpu_multi_threaded()
{
	std::cout << "GPU Multi Threaded Mode Currently Not Supported." << std::endl;
}

void SimpleCrossover::breedKp()
{
	uint16_t upper_mask = 0xFF00;
	uint16_t lower_mask = 0x00FF;
	uint16_t parent1_chrom, parent2_chrom, child1_chrom, child2_chrom;

	int parent1_idx, parent2_idx;
	double parent1_data, parent2_data;

	for (int i = 0; i < sc_data->Kp.size() - 1; i += 2)
	{
		/* Grab the two breeding parent index identifiers */
		parent1_idx = sc_parents[i];
		parent2_idx = sc_parents[i + 1];

		//Convert the current parent data to a chromosome
		parent1_chrom = data2Chromosome(mapCF_Kp, sc_data->Kp.data()[parent1_idx]);
		parent2_chrom = data2Chromosome(mapCF_Kp, sc_data->Kp.data()[parent2_idx]);

		//Create new children
		sc_chrom->Kp.data()[parent1_idx] = (parent1_chrom & upper_mask) | (parent2_chrom & lower_mask);
		sc_chrom->Kp.data()[parent2_idx] = (parent1_chrom & lower_mask) | (parent2_chrom & upper_mask);
	}
}

void SimpleCrossover::breedKi()
{
	uint16_t upper_mask = 0xFF00;
	uint16_t lower_mask = 0x00FF;
	uint16_t parent1_chrom, parent2_chrom, child1_chrom, child2_chrom;

	int parent1_idx, parent2_idx;
	double parent1_data, parent2_data;

	for (int i = 0; i < sc_data->Ki.size() - 1; i += 2)
	{
		/* Grab the two breeding parent index identifiers */
		parent1_idx = sc_parents[i];
		parent2_idx = sc_parents[i + 1];

		//Convert the current parent data to a chromosome
		parent1_chrom = data2Chromosome(mapCF_Ki, sc_data->Ki.data()[parent1_idx]);
		parent2_chrom = data2Chromosome(mapCF_Ki, sc_data->Ki.data()[parent2_idx]);

		//Create new children
		sc_chrom->Ki.data()[parent1_idx] = (parent1_chrom & upper_mask) | (parent2_chrom & lower_mask);
		sc_chrom->Ki.data()[parent2_idx] = (parent1_chrom & lower_mask) | (parent2_chrom & upper_mask);
	}
}

void SimpleCrossover::breedKd()
{
	uint16_t upper_mask = 0xFF00;
	uint16_t lower_mask = 0x00FF;
	uint16_t parent1_chrom, parent2_chrom, child1_chrom, child2_chrom;

	int parent1_idx, parent2_idx;
	double parent1_data, parent2_data;

	for (int i = 0; i < sc_data->Kd.size() - 1; i += 2)
	{
		/* Grab the two breeding parent index identifiers */
		parent1_idx = sc_parents[i];
		parent2_idx = sc_parents[i + 1];

		//Convert the current parent data to a chromosome
		parent1_chrom = data2Chromosome(mapCF_Kd, sc_data->Kd.data()[parent1_idx]);
		parent2_chrom = data2Chromosome(mapCF_Kd, sc_data->Kd.data()[parent2_idx]);

		//Create new children
		sc_chrom->Kd.data()[parent1_idx] = (parent1_chrom & upper_mask) | (parent2_chrom & lower_mask);
		sc_chrom->Kd.data()[parent2_idx] = (parent1_chrom & lower_mask) | (parent2_chrom & upper_mask);
	}
}

///////////////////////////////////////////////////
/* CLASS:  DynamicCrossover */
///////////////////////////////////////////////////
/*-----------------------------------------------
* Constructors/Destructor
*-----------------------------------------------*/
DynamicCrossover::DynamicCrossover(GA_RunMode execution_type)
{
	executionType = execution_type;
}

DynamicCrossover::~DynamicCrossover()
{
}
/*-----------------------------------------------
* Public Functions
*-----------------------------------------------*/
void DynamicCrossover::breed(iVec in_parents, hPID_DataVector* in_pidVals, hPID_Chromosomes* out_breed, mapCoeff_t* in_kp, mapCoeff_t* in_ki, mapCoeff_t* in_kd)
{
	dc_parents = in_parents;
	dc_data = in_pidVals;
	dc_chrom = out_breed;

	mapCF_Kp = in_kp;
	mapCF_Ki = in_ki;
	mapCF_Kd = in_kd;

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
void DynamicCrossover::calculate_cpu_single_threaded()
{
	breedKp();
	breedKi();
	breedKd();
}

void DynamicCrossover::calculate_cpu_multi_threaded()
{
	boost::thread_group tgroup;

	tgroup.create_thread(boost::bind(&DynamicCrossover::breedKp, this));
	tgroup.create_thread(boost::bind(&DynamicCrossover::breedKi, this));
	tgroup.create_thread(boost::bind(&DynamicCrossover::breedKd, this));

	tgroup.join_all();
}

void DynamicCrossover::calculate_gpu_single_threaded()
{
	std::cout << "GPU Single Threaded Mode Currently Not Supported." << std::endl;
}

void DynamicCrossover::calculate_gpu_multi_threaded()
{
	std::cout << "GPU Multi Threaded Mode Currently Not Supported." << std::endl;
}

void DynamicCrossover::breedKp()
{
	/* Instantiate the RNG device */
	std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_real_distribution<double> ratioGen(0.0, 1.0);

	double dc_ratio = 0.0;
	uint16_t upper_mask = 0;
	uint16_t lower_mask = 0;
	uint16_t parent1_chrom, parent2_chrom, child1_chrom, child2_chrom;

	int parent1_idx, parent2_idx;
	double parent1_data, parent2_data;

	for (int i = 0; i < dc_data->Kp.size() - 1; i += 2)
	{
		/* Choose genetic exchange ratio */
		dc_ratio = ratioGen(rng);
		upper_mask = (uint16_t)floor(dc_ratio*((double)0xFFFF));
		lower_mask = ~upper_mask;

		/* Grab the two breeding parent index identifiers */
		parent1_idx = dc_parents[i];
		parent2_idx = dc_parents[i + 1];

		//Convert the current parent data to a chromosome
		parent1_chrom = data2Chromosome(mapCF_Kp, dc_data->Kp.data()[parent1_idx]);
		parent2_chrom = data2Chromosome(mapCF_Kp, dc_data->Kp.data()[parent2_idx]);

		//Create new children
		dc_chrom->Kp.data()[parent1_idx] = (parent1_chrom & upper_mask) | (parent2_chrom & lower_mask);
		dc_chrom->Kp.data()[parent2_idx] = (parent1_chrom & lower_mask) | (parent2_chrom & upper_mask);
	}
}

void DynamicCrossover::breedKi()
{
	/* Instantiate the RNG device */
	std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_real_distribution<double> ratioGen(0.0, 1.0);

	double dc_ratio = 0.0;
	uint16_t upper_mask = 0;
	uint16_t lower_mask = 0;
	uint16_t parent1_chrom, parent2_chrom, child1_chrom, child2_chrom;

	int parent1_idx, parent2_idx;
	double parent1_data, parent2_data;

	for (int i = 0; i < dc_data->Ki.size() - 1; i += 2)
	{
		/* Choose genetic exchange ratio */
		dc_ratio = ratioGen(rng);
		upper_mask = (uint16_t)floor(dc_ratio*((double)0xFFFF));
		lower_mask = ~upper_mask;

		/* Grab the two breeding parent index identifiers */
		parent1_idx = dc_parents[i];
		parent2_idx = dc_parents[i + 1];

		//Convert the current parent data to a chromosome
		parent1_chrom = data2Chromosome(mapCF_Ki, dc_data->Ki.data()[parent1_idx]);
		parent2_chrom = data2Chromosome(mapCF_Ki, dc_data->Ki.data()[parent2_idx]);

		//Create new children
		dc_chrom->Ki.data()[parent1_idx] = (parent1_chrom & upper_mask) | (parent2_chrom & lower_mask);
		dc_chrom->Ki.data()[parent2_idx] = (parent1_chrom & lower_mask) | (parent2_chrom & upper_mask);
	}
}

void DynamicCrossover::breedKd()
{
	/* Instantiate the RNG device */
	std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_real_distribution<double> ratioGen(0.0, 1.0);

	double dc_ratio = 0.0;
	uint16_t upper_mask = 0;
	uint16_t lower_mask = 0;
	uint16_t parent1_chrom, parent2_chrom, child1_chrom, child2_chrom;

	int parent1_idx, parent2_idx;
	double parent1_data, parent2_data;

	for (int i = 0; i < dc_data->Kd.size() - 1; i += 2)
	{
		/* Choose genetic exchange ratio */
		dc_ratio = ratioGen(rng);
		upper_mask = (uint16_t)floor(dc_ratio*((double)0xFFFF));
		lower_mask = ~upper_mask;

		/* Grab the two breeding parent index identifiers */
		parent1_idx = dc_parents[i];
		parent2_idx = dc_parents[i + 1];

		//Convert the current parent data to a chromosome
		parent1_chrom = data2Chromosome(mapCF_Kd, dc_data->Kd.data()[parent1_idx]);
		parent2_chrom = data2Chromosome(mapCF_Kd, dc_data->Kd.data()[parent2_idx]);

		//Create new children
		dc_chrom->Kd.data()[parent1_idx] = (parent1_chrom & upper_mask) | (parent2_chrom & lower_mask);
		dc_chrom->Kd.data()[parent2_idx] = (parent1_chrom & lower_mask) | (parent2_chrom & upper_mask);
	}
}


///////////////////////////////////////////////////
/* CLASS:  FixedRatioCrossover */
///////////////////////////////////////////////////
/*-----------------------------------------------
* Constructors/Destructor
*-----------------------------------------------*/
FixedRatioCrossover::FixedRatioCrossover(GA_RunMode execution_type, double ratio)
{
	executionType = execution_type;
	fr_ratio = ratio;
}

FixedRatioCrossover::~FixedRatioCrossover()
{
}
/*-----------------------------------------------
* Public Functions
*-----------------------------------------------*/
void FixedRatioCrossover::breed(iVec in_parents, hPID_DataVector* in_pidVals, hPID_Chromosomes* out_breed, mapCoeff_t* in_kp, mapCoeff_t* in_ki, mapCoeff_t* in_kd)
{
	fr_parents = in_parents;
	fr_data = in_pidVals;
	fr_chrom = out_breed;

	mapCF_Kp = in_kp;
	mapCF_Ki = in_ki;
	mapCF_Kd = in_kd;

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
void FixedRatioCrossover::calculate_cpu_single_threaded()
{
	breedKp();
	breedKi();
	breedKd();
}

void FixedRatioCrossover::calculate_cpu_multi_threaded()
{
	boost::thread_group tgroup;

	tgroup.create_thread(boost::bind(&FixedRatioCrossover::breedKp, this));
	tgroup.create_thread(boost::bind(&FixedRatioCrossover::breedKi, this));
	tgroup.create_thread(boost::bind(&FixedRatioCrossover::breedKd, this));

	tgroup.join_all();
}

void FixedRatioCrossover::calculate_gpu_single_threaded()
{
	std::cout << "GPU Single Threaded Mode Currently Not Supported." << std::endl;
}

void FixedRatioCrossover::calculate_gpu_multi_threaded()
{
	std::cout << "GPU Multi Threaded Mode Currently Not Supported." << std::endl;
}

void FixedRatioCrossover::breedKp()
{
	uint16_t upper_mask = (uint16_t)floor(fr_ratio*((double)0xFFFF));
	uint16_t lower_mask = ~upper_mask;

	uint16_t parent1_chrom, parent2_chrom, child1_chrom, child2_chrom;

	int parent1_idx, parent2_idx;
	double parent1_data, parent2_data;

	for (int i = 0; i < fr_data->Kp.size() - 1; i += 2)
	{
		/* Grab the two breeding parent index identifiers */
		parent1_idx = fr_parents[i];
		parent2_idx = fr_parents[i + 1];

		//Convert the current parent data to a chromosome
		parent1_chrom = data2Chromosome(mapCF_Kp, fr_data->Kp.data()[parent1_idx]);
		parent2_chrom = data2Chromosome(mapCF_Kp, fr_data->Kp.data()[parent2_idx]);

		//Create new children
		fr_chrom->Kp.data()[parent1_idx] = (parent1_chrom & upper_mask) | (parent2_chrom & lower_mask);
		fr_chrom->Kp.data()[parent2_idx] = (parent1_chrom & lower_mask) | (parent2_chrom & upper_mask);
	}
}

void FixedRatioCrossover::breedKi()
{
	uint16_t upper_mask = (uint16_t)floor(fr_ratio*((double)0xFFFF));
	uint16_t lower_mask = ~upper_mask;

	uint16_t parent1_chrom, parent2_chrom, child1_chrom, child2_chrom;

	int parent1_idx, parent2_idx;
	double parent1_data, parent2_data;

	for (int i = 0; i < fr_data->Ki.size() - 1; i += 2)
	{
		/* Grab the two breeding parent index identifiers */
		parent1_idx = fr_parents[i];
		parent2_idx = fr_parents[i + 1];

		//Convert the current parent data to a chromosome
		parent1_chrom = data2Chromosome(mapCF_Ki, fr_data->Ki.data()[parent1_idx]);
		parent2_chrom = data2Chromosome(mapCF_Ki, fr_data->Ki.data()[parent2_idx]);

		//Create new children
		fr_chrom->Ki.data()[parent1_idx] = (parent1_chrom & upper_mask) | (parent2_chrom & lower_mask);
		fr_chrom->Ki.data()[parent2_idx] = (parent1_chrom & lower_mask) | (parent2_chrom & upper_mask);
	}
}

void FixedRatioCrossover::breedKd()
{
	uint16_t upper_mask = (uint16_t)floor(fr_ratio*((double)0xFFFF));
	uint16_t lower_mask = ~upper_mask;

	uint16_t parent1_chrom, parent2_chrom, child1_chrom, child2_chrom;

	int parent1_idx, parent2_idx;
	double parent1_data, parent2_data;

	for (int i = 0; i < fr_data->Kd.size() - 1; i += 2)
	{
		/* Grab the two breeding parent index identifiers */
		parent1_idx = fr_parents[i];
		parent2_idx = fr_parents[i + 1];

		//Convert the current parent data to a chromosome
		parent1_chrom = data2Chromosome(mapCF_Kd, fr_data->Kd.data()[parent1_idx]);
		parent2_chrom = data2Chromosome(mapCF_Kd, fr_data->Kd.data()[parent2_idx]);

		//Create new children
		fr_chrom->Kd.data()[parent1_idx] = (parent1_chrom & upper_mask) | (parent2_chrom & lower_mask);
		fr_chrom->Kd.data()[parent2_idx] = (parent1_chrom & lower_mask) | (parent2_chrom & upper_mask);
	}
}

///////////////////////////////////////////////////
/* CLASS:  SimulatedBinaryCrossover */
///////////////////////////////////////////////////
/*-----------------------------------------------
* Constructors/Destructor
*-----------------------------------------------*/
SimulatedBinaryCrossover::SimulatedBinaryCrossover(GA_RunMode execution_type)
{
	executionType = execution_type;
}

SimulatedBinaryCrossover::~SimulatedBinaryCrossover()
{
}
/*-----------------------------------------------
* Public Functions
*-----------------------------------------------*/

/*-----------------------------------------------
* Private Functions
*-----------------------------------------------*/
void SimulatedBinaryCrossover::calculate_cpu_single_threaded()
{
}

void SimulatedBinaryCrossover::calculate_cpu_multi_threaded()
{
	std::cout << "CPU Multi Threaded Mode Currently Not Supported." << std::endl;
}

void SimulatedBinaryCrossover::calculate_gpu_single_threaded()
{
	std::cout << "GPU Single Threaded Mode Currently Not Supported." << std::endl;
}

void SimulatedBinaryCrossover::calculate_gpu_multi_threaded()
{
	std::cout << "GPU Multi Threaded Mode Currently Not Supported." << std::endl;
}