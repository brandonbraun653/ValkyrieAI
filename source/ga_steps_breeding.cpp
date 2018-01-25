#include "ga_steps_breeding.h"

/* Creates a bit mask by setting bits "a" through "b" to 1 */
#define BIT_MASK_RANGED(a, b) (((unsigned) -1 >> (31-(b))) & ~((1U << (a)) - 1))


///////////////////////////////////////////////////
/* CLASS:  SimpleCrossover */
///////////////////////////////////////////////////
/*-----------------------------------------------
* Constructors/Destructor
*-----------------------------------------------*/
SimpleCrossover::SimpleCrossover()
{
}

SimpleCrossover::~SimpleCrossover()
{
}

/*-----------------------------------------------
* Public Functions
*-----------------------------------------------*/
void SimpleCrossover::breed(GA_BreedingDataInput input, GA_BreedingDataOutput& output)
{
	uint16_t upper_mask = 0xFF00;
	uint16_t lower_mask = 0x00FF;
	uint16_t parent1_chrom, parent2_chrom;

	int parent1_idx, parent2_idx;

	/* Simple crossover outputs data as a bit field regardless of input type */
	output.chromType = MAPPING_TYPE_BIT_FIELD;

	if (input.chromType == MAPPING_TYPE_REAL)
	{
		/* Ensure the output has been correctly allocated */
		if (output.u16_chrom.size() != input.d_chrom.size())
			output.u16_chrom.resize(input.d_chrom.size());


		/* First, convert everything over to a bitfield. Not all population members will be selected 
		to breed. This ensures that every PID field has an accurate value whether it was mated or not. */
		for (int member = 0; member < input.d_chrom.size(); member++)
		{
			output.u16_chrom[member].Kp = data2Chromosome(input.mapCoeff_Kp, input.d_chrom[member].Kp);
			output.u16_chrom[member].Ki = data2Chromosome(input.mapCoeff_Ki, input.d_chrom[member].Ki);
			output.u16_chrom[member].Kd = data2Chromosome(input.mapCoeff_Kd, input.d_chrom[member].Kd);
		}


		/* Now perform the actual mating procedure. Indexing limits must be size()-1 to account for
		the +2 increment */

		if (input.swap_both_chrom_halves)
		{
			/*
			* EXAMPLE:
			* Chromosome 1: [a b c d | e f g h]
			* Chromosome 2: [i j k l | m n o p]
			* --CROSSOVER--
			* Chromosome 1: [a b c d | m n o p]
			* Chromosome 2: [e f g h | i j k l]
			**/

			for (int parent = 0; parent < input.parentSelections.size() - 1; parent += 2)
			{
				parent1_idx = input.parentSelections[parent];
				parent2_idx = input.parentSelections[parent + 1];

				//KP 
				parent1_chrom = output.u16_chrom[parent1_idx].Kp;
				parent2_chrom = output.u16_chrom[parent2_idx].Kp;

				output.u16_chrom[parent1_idx].Kp = (parent1_chrom & upper_mask) | (parent2_chrom & lower_mask);
				output.u16_chrom[parent2_idx].Kp = (parent1_chrom & lower_mask) | (parent2_chrom & upper_mask);


				//Ki 
				parent1_chrom = output.u16_chrom[parent1_idx].Ki;
				parent2_chrom = output.u16_chrom[parent2_idx].Ki;

				output.u16_chrom[parent1_idx].Ki = (parent1_chrom & upper_mask) | (parent2_chrom & lower_mask);
				output.u16_chrom[parent2_idx].Ki = (parent1_chrom & lower_mask) | (parent2_chrom & upper_mask);


				//Kd 
				parent1_chrom = output.u16_chrom[parent1_idx].Kd;
				parent2_chrom = output.u16_chrom[parent2_idx].Kd;

				output.u16_chrom[parent1_idx].Kd = (parent1_chrom & upper_mask) | (parent2_chrom & lower_mask);
				output.u16_chrom[parent2_idx].Kd = (parent1_chrom & lower_mask) | (parent2_chrom & upper_mask);
			}
		}
		else
		{
			/*
			* EXAMPLE: (swap_lower_chrom_half == true)
			* Chromosome 1: [a b c d | e f g h]
			* Chromosome 2: [i j k l | m n o p]
			* --CROSSOVER--
			* Chromosome 1: [a b c d | m n o p]
			* Chromosome 2: [i j k l | e f g h]
			**/

			/* Reverses the gene swapping directions if true. Default is to only swap the LSB masked bits. */
			if (!input.swap_lower_chrom_half)
			{
				upper_mask = ~upper_mask;
				lower_mask = ~lower_mask;
			}

			for (int parent = 0; parent < input.parentSelections.size() - 1; parent += 2)
			{
				parent1_idx = input.parentSelections[parent];
				parent2_idx = input.parentSelections[parent + 1];

				//KP 
				parent1_chrom = output.u16_chrom[parent1_idx].Kp;
				parent2_chrom = output.u16_chrom[parent2_idx].Kp;

				//										save upper bits					replace lower bits 
				output.u16_chrom[parent1_idx].Kp = (parent1_chrom & upper_mask) | (parent2_chrom & lower_mask);
				output.u16_chrom[parent2_idx].Kp = (parent2_chrom & upper_mask) | (parent1_chrom & lower_mask);


				//Ki 
				parent1_chrom = output.u16_chrom[parent1_idx].Ki;
				parent2_chrom = output.u16_chrom[parent2_idx].Ki;

				output.u16_chrom[parent1_idx].Ki |= (parent1_chrom & upper_mask) | (parent2_chrom & lower_mask);
				output.u16_chrom[parent2_idx].Ki |= (parent2_chrom & upper_mask) | (parent1_chrom & lower_mask);

				//Kd 
				parent1_chrom = output.u16_chrom[parent1_idx].Kd;
				parent2_chrom = output.u16_chrom[parent2_idx].Kd;

				output.u16_chrom[parent1_idx].Kd |= (parent1_chrom & upper_mask) | (parent2_chrom & lower_mask);
				output.u16_chrom[parent2_idx].Kd |= (parent2_chrom & upper_mask) | (parent1_chrom & lower_mask);
			}
		}
	}
	else if (input.chromType == MAPPING_TYPE_BIT_FIELD)
	{
		std::cout << "Using bit field map for breeding...not currently implemented!!" << std::endl;
	}
}

/*-----------------------------------------------
* Private Functions
*-----------------------------------------------*/
void SimpleCrossover::breedKp()
{

}

void SimpleCrossover::breedKi()
{

}

void SimpleCrossover::breedKd()
{

}


///////////////////////////////////////////////////
/* CLASS:  DynamicCrossover */
///////////////////////////////////////////////////
/*-----------------------------------------------
* Constructors/Destructor
*-----------------------------------------------*/
DynamicCrossover::DynamicCrossover()
{
}

DynamicCrossover::~DynamicCrossover()
{
}

/*-----------------------------------------------
* Public Functions
*-----------------------------------------------*/
void DynamicCrossover::breed(GA_BreedingDataInput input, GA_BreedingDataOutput& output)
{
	std::cout << "Hello from the dynamic crossover breeding function!" << std::endl;
}

/*-----------------------------------------------
* Private Functions
*-----------------------------------------------*/
void DynamicCrossover::breedKp()
{

}

void DynamicCrossover::breedKi()
{

}

void DynamicCrossover::breedKd()
{

}


///////////////////////////////////////////////////
/* CLASS:  FixedRatioCrossover */
///////////////////////////////////////////////////
/*-----------------------------------------------
* Constructors/Destructor
*-----------------------------------------------*/
FixedPointCrossover::FixedPointCrossover()
{
}

FixedPointCrossover::~FixedPointCrossover()
{
}

/*-----------------------------------------------
* Public Functions
*-----------------------------------------------*/
void FixedPointCrossover::breed(GA_BreedingDataInput input, GA_BreedingDataOutput& output)
{
	uint16_t lower_mask = BIT_MASK_RANGED(0, input.crossoverPoint-1);
	uint16_t upper_mask = ~lower_mask;
	uint16_t parent1_chrom, parent2_chrom;

	int parent1_idx, parent2_idx;

	/* Simple crossover outputs data as a bit field regardless of input type */
	output.chromType = MAPPING_TYPE_BIT_FIELD;

	if (input.chromType == MAPPING_TYPE_REAL)
	{
		/* Ensure the output has been correctly allocated */
		if (output.u16_chrom.size() != input.d_chrom.size())
			output.u16_chrom.resize(input.d_chrom.size());


		/* First, convert everything over to a bitfield. Not all population members will be selected
		to breed. This ensures that every PID field has an accurate value whether it was mated or not. */
		for (int member = 0; member < input.d_chrom.size(); member++)
		{
			output.u16_chrom[member].Kp = data2Chromosome(input.mapCoeff_Kp, input.d_chrom[member].Kp);
			output.u16_chrom[member].Ki = data2Chromosome(input.mapCoeff_Ki, input.d_chrom[member].Ki);
			output.u16_chrom[member].Kd = data2Chromosome(input.mapCoeff_Kd, input.d_chrom[member].Kd);
		}


		/* Now perform the actual mating procedure. The loops are split up to get around evaluating the
		   input.swap_both_chroms variable many times unnecessarily. Indexing limits must be size()-1 to account for
		   the +2 increment */
		if (input.swap_both_chrom_halves)
		{
			/*
			* EXAMPLE:
			* Chromosome 1: [a b c d | e f g h]
			* Chromosome 2: [i j k l | m n o p]
			* --CROSSOVER--
			* Chromosome 1: [a b c d | m n o p]
			* Chromosome 2: [e f g h | i j k l]
			**/

			for (int parent = 0; parent < input.parentSelections.size() - 1; parent += 2)
			{
				parent1_idx = input.parentSelections[parent];
				parent2_idx = input.parentSelections[parent + 1];

				//KP 
				parent1_chrom = output.u16_chrom[parent1_idx].Kp;
				parent2_chrom = output.u16_chrom[parent2_idx].Kp;

				output.u16_chrom[parent1_idx].Kp = (parent1_chrom & upper_mask) | (parent2_chrom & lower_mask);
				output.u16_chrom[parent2_idx].Kp = (parent1_chrom & lower_mask) | (parent2_chrom & upper_mask);


				//Ki 
				parent1_chrom = output.u16_chrom[parent1_idx].Ki;
				parent2_chrom = output.u16_chrom[parent2_idx].Ki;

				output.u16_chrom[parent1_idx].Ki = (parent1_chrom & upper_mask) | (parent2_chrom & lower_mask);
				output.u16_chrom[parent2_idx].Ki = (parent1_chrom & lower_mask) | (parent2_chrom & upper_mask);


				//Kd 
				parent1_chrom = output.u16_chrom[parent1_idx].Kd;
				parent2_chrom = output.u16_chrom[parent2_idx].Kd;

				output.u16_chrom[parent1_idx].Kd = (parent1_chrom & upper_mask) | (parent2_chrom & lower_mask);
				output.u16_chrom[parent2_idx].Kd = (parent1_chrom & lower_mask) | (parent2_chrom & upper_mask);
			}
		}
		else
		{
			/*
			* EXAMPLE: (swap_lower_chrom_half == true)
			* Chromosome 1: [a b c d | e f g h]
			* Chromosome 2: [i j k l | m n o p]
			* --CROSSOVER--
			* Chromosome 1: [a b c d | m n o p]
			* Chromosome 2: [i j k l | e f g h]
			**/

			/* Reverses the gene swapping directions if true. Default is to only swap the LSB masked bits. */
			if (!input.swap_lower_chrom_half)
			{
				upper_mask = ~upper_mask;
				lower_mask = ~lower_mask;
			}

			for (int parent = 0; parent < input.parentSelections.size() - 1; parent += 2)
			{
				parent1_idx = input.parentSelections[parent];
				parent2_idx = input.parentSelections[parent + 1];

				//KP 
				parent1_chrom = output.u16_chrom[parent1_idx].Kp;
				parent2_chrom = output.u16_chrom[parent2_idx].Kp;

				//										save upper bits					replace lower bits 
				output.u16_chrom[parent1_idx].Kp = (parent1_chrom & upper_mask) | (parent2_chrom & lower_mask);
				output.u16_chrom[parent2_idx].Kp = (parent2_chrom & upper_mask) | (parent1_chrom & lower_mask);


				//Ki 
				parent1_chrom = output.u16_chrom[parent1_idx].Ki;
				parent2_chrom = output.u16_chrom[parent2_idx].Ki;

				output.u16_chrom[parent1_idx].Ki |= (parent1_chrom & upper_mask) | (parent2_chrom & lower_mask);
				output.u16_chrom[parent2_idx].Ki |= (parent2_chrom & upper_mask) | (parent1_chrom & lower_mask);

				//Kd 
				parent1_chrom = output.u16_chrom[parent1_idx].Kd;
				parent2_chrom = output.u16_chrom[parent2_idx].Kd;

				output.u16_chrom[parent1_idx].Kd |= (parent1_chrom & upper_mask) | (parent2_chrom & lower_mask);
				output.u16_chrom[parent2_idx].Kd |= (parent2_chrom & upper_mask) | (parent1_chrom & lower_mask);
			}
		}
	}
	else if (input.chromType == MAPPING_TYPE_BIT_FIELD)
	{
		std::cout << "Using bit field map input for breeding...not currently implemented!!" << std::endl;
	}
}

/*-----------------------------------------------
* Private Functions
*-----------------------------------------------*/
void FixedPointCrossover::breedKp()
{

}

void FixedPointCrossover::breedKi()
{

}

void FixedPointCrossover::breedKd()
{

}


///////////////////////////////////////////////////
/* CLASS:  SimulatedBinaryCrossover */
///////////////////////////////////////////////////
/*-----------------------------------------------
* Constructors/Destructor
*-----------------------------------------------*/
SimulatedBinaryCrossover::SimulatedBinaryCrossover()
{
}

SimulatedBinaryCrossover::~SimulatedBinaryCrossover()
{
}

/*-----------------------------------------------
* Public Functions
*-----------------------------------------------*/
void SimulatedBinaryCrossover::breed(GA_BreedingDataInput input, GA_BreedingDataOutput& output)
{
	std::cout << "Hello from the simulated binary crossover breeding function!" << std::endl;
}

/*-----------------------------------------------
* Private Functions
*-----------------------------------------------*/
void SimulatedBinaryCrossover::breedKp()
{

}

void SimulatedBinaryCrossover::breedKi()
{

}

void SimulatedBinaryCrossover::breedKd()
{

}


// void SimpleCrossover::calculate_cpu_single_threaded()
// {
// 	breedKp();
// 	breedKi();
// 	breedKd();
// }
// 
// void SimpleCrossover::calculate_cpu_multi_threaded()
// {
// 	boost::thread_group tgroup;
// 
// 	tgroup.create_thread(boost::bind(&SimpleCrossover::breedKp, this));
// 	tgroup.create_thread(boost::bind(&SimpleCrossover::breedKi, this));
// 	tgroup.create_thread(boost::bind(&SimpleCrossover::breedKd, this));
// 
// 	tgroup.join_all();
// }
// 
// void SimpleCrossover::calculate_gpu_single_threaded()
// {
// 	std::cout << "GPU Single Threaded Mode Currently Not Supported." << std::endl;
// }
// 
// void SimpleCrossover::calculate_gpu_multi_threaded()
// {
// 	std::cout << "GPU Multi Threaded Mode Currently Not Supported." << std::endl;
// }
// 
// void SimpleCrossover::breedKp()
// {
// 	uint16_t upper_mask = 0xFF00;
// 	uint16_t lower_mask = 0x00FF;
// 	uint16_t parent1_chrom, parent2_chrom, child1_chrom, child2_chrom;
// 
// 	int parent1_idx, parent2_idx;
// 	double parent1_data, parent2_data;
// 
// 	for (int i = 0; i < sc_data->Kp.size() - 1; i += 2)
// 	{
// 		/* Grab the two breeding parent index identifiers */
// 		parent1_idx = sc_parents[i];
// 		parent2_idx = sc_parents[i + 1];
// 
// 		//Convert the current parent data to a chromosome
// 		parent1_chrom = data2Chromosome(mapCF_Kp, sc_data->Kp.data()[parent1_idx]);
// 		parent2_chrom = data2Chromosome(mapCF_Kp, sc_data->Kp.data()[parent2_idx]);
// 
// 		//Create new children
// 		sc_chrom->Kp.data()[parent1_idx] = (parent1_chrom & upper_mask) | (parent2_chrom & lower_mask);
// 		sc_chrom->Kp.data()[parent2_idx] = (parent1_chrom & lower_mask) | (parent2_chrom & upper_mask);
// 	}
// }
// 
// void SimpleCrossover::breedKi()
// {
// 	uint16_t upper_mask = 0xFF00;
// 	uint16_t lower_mask = 0x00FF;
// 	uint16_t parent1_chrom, parent2_chrom, child1_chrom, child2_chrom;
// 
// 	int parent1_idx, parent2_idx;
// 	double parent1_data, parent2_data;
// 
// 	for (int i = 0; i < sc_data->Ki.size() - 1; i += 2)
// 	{
// 		/* Grab the two breeding parent index identifiers */
// 		parent1_idx = sc_parents[i];
// 		parent2_idx = sc_parents[i + 1];
// 
// 		//Convert the current parent data to a chromosome
// 		parent1_chrom = data2Chromosome(mapCF_Ki, sc_data->Ki.data()[parent1_idx]);
// 		parent2_chrom = data2Chromosome(mapCF_Ki, sc_data->Ki.data()[parent2_idx]);
// 
// 		//Create new children
// 		sc_chrom->Ki.data()[parent1_idx] = (parent1_chrom & upper_mask) | (parent2_chrom & lower_mask);
// 		sc_chrom->Ki.data()[parent2_idx] = (parent1_chrom & lower_mask) | (parent2_chrom & upper_mask);
// 	}
// }
// 
// void SimpleCrossover::breedKd()
// {
// 	uint16_t upper_mask = 0xFF00;
// 	uint16_t lower_mask = 0x00FF;
// 	uint16_t parent1_chrom, parent2_chrom, child1_chrom, child2_chrom;
// 
// 	int parent1_idx, parent2_idx;
// 	double parent1_data, parent2_data;
// 
// 	for (int i = 0; i < sc_data->Kd.size() - 1; i += 2)
// 	{
// 		/* Grab the two breeding parent index identifiers */
// 		parent1_idx = sc_parents[i];
// 		parent2_idx = sc_parents[i + 1];
// 
// 		//Convert the current parent data to a chromosome
// 		parent1_chrom = data2Chromosome(mapCF_Kd, sc_data->Kd.data()[parent1_idx]);
// 		parent2_chrom = data2Chromosome(mapCF_Kd, sc_data->Kd.data()[parent2_idx]);
// 
// 		//Create new children
// 		sc_chrom->Kd.data()[parent1_idx] = (parent1_chrom & upper_mask) | (parent2_chrom & lower_mask);
// 		sc_chrom->Kd.data()[parent2_idx] = (parent1_chrom & lower_mask) | (parent2_chrom & upper_mask);
// 	}
// }
// 
// ///////////////////////////////////////////////////
// /* CLASS:  DynamicCrossover */
// ///////////////////////////////////////////////////
// /*-----------------------------------------------
// * Constructors/Destructor
// *-----------------------------------------------*/
// DynamicCrossover::DynamicCrossover()
// {
// }
// 
// DynamicCrossover::~DynamicCrossover()
// {
// }
// /*-----------------------------------------------
// * Public Functions
// *-----------------------------------------------*/
// void DynamicCrossover::breed(iVec in_parents, hPID_DataVector* in_pidVals, hPID_Chromosomes* out_breed, FCSOptimizer_MappingCoeff* in_kp, FCSOptimizer_MappingCoeff* in_ki, FCSOptimizer_MappingCoeff* in_kd)
// {
// 	dc_parents = in_parents;
// 	dc_data = in_pidVals;
// 	dc_chrom = out_breed;
// 
// 	mapCF_Kp = in_kp;
// 	mapCF_Ki = in_ki;
// 	mapCF_Kd = in_kd;
// 
// // 	if (executionType == SINGLE_THREADED)
// // 		calculate_cpu_single_threaded();
// // 
// // 	if (executionType == MULTI_THREADED)
// // 		calculate_cpu_multi_threaded();
// // 
// // 	if (executionType == SINGLE_THREADED_WITH_CUDA)
// // 		calculate_gpu_single_threaded();
// // 
// // 	if (executionType == MULTI_THREADED_WITH_CUDA)
// // 		calculate_gpu_multi_threaded();
// }
// 
// /*-----------------------------------------------
// * Private Functions
// *-----------------------------------------------*/
// void DynamicCrossover::calculate_cpu_single_threaded()
// {
// 	breedKp();
// 	breedKi();
// 	breedKd();
// }
// 
// void DynamicCrossover::calculate_cpu_multi_threaded()
// {
// 	boost::thread_group tgroup;
// 
// 	tgroup.create_thread(boost::bind(&DynamicCrossover::breedKp, this));
// 	tgroup.create_thread(boost::bind(&DynamicCrossover::breedKi, this));
// 	tgroup.create_thread(boost::bind(&DynamicCrossover::breedKd, this));
// 
// 	tgroup.join_all();
// }
// 
// void DynamicCrossover::calculate_gpu_single_threaded()
// {
// 	std::cout << "GPU Single Threaded Mode Currently Not Supported." << std::endl;
// }
// 
// void DynamicCrossover::calculate_gpu_multi_threaded()
// {
// 	std::cout << "GPU Multi Threaded Mode Currently Not Supported." << std::endl;
// }
// 
// void DynamicCrossover::breedKp()
// {
// 	/* Instantiate the RNG device */
// 	std::random_device rd;
// 	std::mt19937 rng(rd());
// 	std::uniform_real_distribution<double> ratioGen(0.0, 1.0);
// 
// 	double dc_ratio = 0.0;
// 	uint16_t upper_mask = 0;
// 	uint16_t lower_mask = 0;
// 	uint16_t parent1_chrom, parent2_chrom, child1_chrom, child2_chrom;
// 
// 	int parent1_idx, parent2_idx;
// 	double parent1_data, parent2_data;
// 
// 	for (int i = 0; i < dc_data->Kp.size() - 1; i += 2)
// 	{
// 		/* Choose genetic exchange ratio */
// 		dc_ratio = ratioGen(rng);
// 		upper_mask = (uint16_t)floor(dc_ratio*((double)0xFFFF));
// 		lower_mask = ~upper_mask;
// 
// 		/* Grab the two breeding parent index identifiers */
// 		parent1_idx = dc_parents[i];
// 		parent2_idx = dc_parents[i + 1];
// 
// 		//Convert the current parent data to a chromosome
// 		parent1_chrom = data2Chromosome(mapCF_Kp, dc_data->Kp.data()[parent1_idx]);
// 		parent2_chrom = data2Chromosome(mapCF_Kp, dc_data->Kp.data()[parent2_idx]);
// 
// 		//Create new children
// 		dc_chrom->Kp.data()[parent1_idx] = (parent1_chrom & upper_mask) | (parent2_chrom & lower_mask);
// 		dc_chrom->Kp.data()[parent2_idx] = (parent1_chrom & lower_mask) | (parent2_chrom & upper_mask);
// 	}
// }
// 
// void DynamicCrossover::breedKi()
// {
// 	/* Instantiate the RNG device */
// 	std::random_device rd;
// 	std::mt19937 rng(rd());
// 	std::uniform_real_distribution<double> ratioGen(0.0, 1.0);
// 
// 	double dc_ratio = 0.0;
// 	uint16_t upper_mask = 0;
// 	uint16_t lower_mask = 0;
// 	uint16_t parent1_chrom, parent2_chrom, child1_chrom, child2_chrom;
// 
// 	int parent1_idx, parent2_idx;
// 	double parent1_data, parent2_data;
// 
// 	for (int i = 0; i < dc_data->Ki.size() - 1; i += 2)
// 	{
// 		/* Choose genetic exchange ratio */
// 		dc_ratio = ratioGen(rng);
// 		upper_mask = (uint16_t)floor(dc_ratio*((double)0xFFFF));
// 		lower_mask = ~upper_mask;
// 
// 		/* Grab the two breeding parent index identifiers */
// 		parent1_idx = dc_parents[i];
// 		parent2_idx = dc_parents[i + 1];
// 
// 		//Convert the current parent data to a chromosome
// 		parent1_chrom = data2Chromosome(mapCF_Ki, dc_data->Ki.data()[parent1_idx]);
// 		parent2_chrom = data2Chromosome(mapCF_Ki, dc_data->Ki.data()[parent2_idx]);
// 
// 		//Create new children
// 		dc_chrom->Ki.data()[parent1_idx] = (parent1_chrom & upper_mask) | (parent2_chrom & lower_mask);
// 		dc_chrom->Ki.data()[parent2_idx] = (parent1_chrom & lower_mask) | (parent2_chrom & upper_mask);
// 	}
// }
// 
// void DynamicCrossover::breedKd()
// {
// 	/* Instantiate the RNG device */
// 	std::random_device rd;
// 	std::mt19937 rng(rd());
// 	std::uniform_real_distribution<double> ratioGen(0.0, 1.0);
// 
// 	double dc_ratio = 0.0;
// 	uint16_t upper_mask = 0;
// 	uint16_t lower_mask = 0;
// 	uint16_t parent1_chrom, parent2_chrom, child1_chrom, child2_chrom;
// 
// 	int parent1_idx, parent2_idx;
// 	double parent1_data, parent2_data;
// 
// 	for (int i = 0; i < dc_data->Kd.size() - 1; i += 2)
// 	{
// 		/* Choose genetic exchange ratio */
// 		dc_ratio = ratioGen(rng);
// 		upper_mask = (uint16_t)floor(dc_ratio*((double)0xFFFF));
// 		lower_mask = ~upper_mask;
// 
// 		/* Grab the two breeding parent index identifiers */
// 		parent1_idx = dc_parents[i];
// 		parent2_idx = dc_parents[i + 1];
// 
// 		//Convert the current parent data to a chromosome
// 		parent1_chrom = data2Chromosome(mapCF_Kd, dc_data->Kd.data()[parent1_idx]);
// 		parent2_chrom = data2Chromosome(mapCF_Kd, dc_data->Kd.data()[parent2_idx]);
// 
// 		//Create new children
// 		dc_chrom->Kd.data()[parent1_idx] = (parent1_chrom & upper_mask) | (parent2_chrom & lower_mask);
// 		dc_chrom->Kd.data()[parent2_idx] = (parent1_chrom & lower_mask) | (parent2_chrom & upper_mask);
// 	}
// }
// 
// 
// ///////////////////////////////////////////////////
// /* CLASS:  FixedRatioCrossover */
// ///////////////////////////////////////////////////
// /*-----------------------------------------------
// * Constructors/Destructor
// *-----------------------------------------------*/
// FixedRatioCrossover::FixedRatioCrossover(double ratio)
// {
// 	fr_ratio = ratio;
// }
// 
// FixedRatioCrossover::~FixedRatioCrossover()
// {
// }
// /*-----------------------------------------------
// * Public Functions
// *-----------------------------------------------*/
// void FixedRatioCrossover::breed(iVec in_parents, hPID_DataVector* in_pidVals, hPID_Chromosomes* out_breed, FCSOptimizer_MappingCoeff* in_kp, FCSOptimizer_MappingCoeff* in_ki, FCSOptimizer_MappingCoeff* in_kd)
// {
// 	fr_parents = in_parents;
// 	fr_data = in_pidVals;
// 	fr_chrom = out_breed;
// 
// 	mapCF_Kp = in_kp;
// 	mapCF_Ki = in_ki;
// 	mapCF_Kd = in_kd;
// 
// // 	if (executionType == SINGLE_THREADED)
// // 		calculate_cpu_single_threaded();
// // 
// // 	if (executionType == MULTI_THREADED)
// // 		calculate_cpu_multi_threaded();
// // 
// // 	if (executionType == SINGLE_THREADED_WITH_CUDA)
// // 		calculate_gpu_single_threaded();
// // 
// // 	if (executionType == MULTI_THREADED_WITH_CUDA)
// // 		calculate_gpu_multi_threaded();
// }
// 
// /*-----------------------------------------------
// * Private Functions
// *-----------------------------------------------*/
// void FixedRatioCrossover::calculate_cpu_single_threaded()
// {
// 	breedKp();
// 	breedKi();
// 	breedKd();
// }
// 
// void FixedRatioCrossover::calculate_cpu_multi_threaded()
// {
// 	boost::thread_group tgroup;
// 
// 	tgroup.create_thread(boost::bind(&FixedRatioCrossover::breedKp, this));
// 	tgroup.create_thread(boost::bind(&FixedRatioCrossover::breedKi, this));
// 	tgroup.create_thread(boost::bind(&FixedRatioCrossover::breedKd, this));
// 
// 	tgroup.join_all();
// }
// 
// void FixedRatioCrossover::calculate_gpu_single_threaded()
// {
// 	std::cout << "GPU Single Threaded Mode Currently Not Supported." << std::endl;
// }
// 
// void FixedRatioCrossover::calculate_gpu_multi_threaded()
// {
// 	std::cout << "GPU Multi Threaded Mode Currently Not Supported." << std::endl;
// }
// 
// void FixedRatioCrossover::breedKp()
// {
// 	uint16_t upper_mask = (uint16_t)floor(fr_ratio*((double)0xFFFF));
// 	uint16_t lower_mask = ~upper_mask;
// 
// 	uint16_t parent1_chrom, parent2_chrom, child1_chrom, child2_chrom;
// 
// 	int parent1_idx, parent2_idx;
// 	double parent1_data, parent2_data;
// 
// 	for (int i = 0; i < fr_data->Kp.size() - 1; i += 2)
// 	{
// 		/* Grab the two breeding parent index identifiers */
// 		parent1_idx = fr_parents[i];
// 		parent2_idx = fr_parents[i + 1];
// 
// 		//Convert the current parent data to a chromosome
// 		parent1_chrom = data2Chromosome(mapCF_Kp, fr_data->Kp.data()[parent1_idx]);
// 		parent2_chrom = data2Chromosome(mapCF_Kp, fr_data->Kp.data()[parent2_idx]);
// 
// 		//Create new children
// 		fr_chrom->Kp.data()[parent1_idx] = (parent1_chrom & upper_mask) | (parent2_chrom & lower_mask);
// 		fr_chrom->Kp.data()[parent2_idx] = (parent1_chrom & lower_mask) | (parent2_chrom & upper_mask);
// 	}
// }
// 
// void FixedRatioCrossover::breedKi()
// {
// 	uint16_t upper_mask = (uint16_t)floor(fr_ratio*((double)0xFFFF));
// 	uint16_t lower_mask = ~upper_mask;
// 
// 	uint16_t parent1_chrom, parent2_chrom, child1_chrom, child2_chrom;
// 
// 	int parent1_idx, parent2_idx;
// 	double parent1_data, parent2_data;
// 
// 	for (int i = 0; i < fr_data->Ki.size() - 1; i += 2)
// 	{
// 		/* Grab the two breeding parent index identifiers */
// 		parent1_idx = fr_parents[i];
// 		parent2_idx = fr_parents[i + 1];
// 
// 		//Convert the current parent data to a chromosome
// 		parent1_chrom = data2Chromosome(mapCF_Ki, fr_data->Ki.data()[parent1_idx]);
// 		parent2_chrom = data2Chromosome(mapCF_Ki, fr_data->Ki.data()[parent2_idx]);
// 
// 		//Create new children
// 		fr_chrom->Ki.data()[parent1_idx] = (parent1_chrom & upper_mask) | (parent2_chrom & lower_mask);
// 		fr_chrom->Ki.data()[parent2_idx] = (parent1_chrom & lower_mask) | (parent2_chrom & upper_mask);
// 	}
// }
// 
// void FixedRatioCrossover::breedKd()
// {
// 	uint16_t upper_mask = (uint16_t)floor(fr_ratio*((double)0xFFFF));
// 	uint16_t lower_mask = ~upper_mask;
// 
// 	uint16_t parent1_chrom, parent2_chrom, child1_chrom, child2_chrom;
// 
// 	int parent1_idx, parent2_idx;
// 	double parent1_data, parent2_data;
// 
// 	for (int i = 0; i < fr_data->Kd.size() - 1; i += 2)
// 	{
// 		/* Grab the two breeding parent index identifiers */
// 		parent1_idx = fr_parents[i];
// 		parent2_idx = fr_parents[i + 1];
// 
// 		//Convert the current parent data to a chromosome
// 		parent1_chrom = data2Chromosome(mapCF_Kd, fr_data->Kd.data()[parent1_idx]);
// 		parent2_chrom = data2Chromosome(mapCF_Kd, fr_data->Kd.data()[parent2_idx]);
// 
// 		//Create new children
// 		fr_chrom->Kd.data()[parent1_idx] = (parent1_chrom & upper_mask) | (parent2_chrom & lower_mask);
// 		fr_chrom->Kd.data()[parent2_idx] = (parent1_chrom & lower_mask) | (parent2_chrom & upper_mask);
// 	}
// }
// 
// ///////////////////////////////////////////////////
// /* CLASS:  SimulatedBinaryCrossover */
// ///////////////////////////////////////////////////
// /*-----------------------------------------------
// * Constructors/Destructor
// *-----------------------------------------------*/
// SimulatedBinaryCrossover::SimulatedBinaryCrossover()
// {
// }
// 
// SimulatedBinaryCrossover::~SimulatedBinaryCrossover()
// {
// }
// /*-----------------------------------------------
// * Public Functions
// *-----------------------------------------------*/
// 
// /*-----------------------------------------------
// * Private Functions
// *-----------------------------------------------*/
// void SimulatedBinaryCrossover::calculate_cpu_single_threaded()
// {
// }
// 
// void SimulatedBinaryCrossover::calculate_cpu_multi_threaded()
// {
// 	std::cout << "CPU Multi Threaded Mode Currently Not Supported." << std::endl;
// }
// 
// void SimulatedBinaryCrossover::calculate_gpu_single_threaded()
// {
// 	std::cout << "GPU Single Threaded Mode Currently Not Supported." << std::endl;
// }
// 
// void SimulatedBinaryCrossover::calculate_gpu_multi_threaded()
// {
// 	std::cout << "GPU Multi Threaded Mode Currently Not Supported." << std::endl;
// }