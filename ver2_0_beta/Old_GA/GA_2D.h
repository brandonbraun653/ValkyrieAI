/* Array Indexers */
#define X 0	//X-Coordinate
#define Y 1 //Y-Coordinate
#define F 2 //Fitness Val
#define NUM_MODELS 8u


//////////////////////////////////////////////////////////////////
/* CLASS: GAClass_2D */
//////////////////////////////////////////////////////////////////
typedef double(*fitfunc_2d_t)(double, double);

struct GA2D_Constraints
{
	double x_upper = 0.0;
	double x_lower = 0.0;
	double y_upper = 0.0;
	double y_lower = 0.0;
	uint32_t x_gridDiv = 10;
	uint32_t y_gridDiv = 10;
	uint32_t numBytes_precision = 2;
};

struct GA2D_RuntimeData
{
	bool locked = false;
	uint32_t iterationNumber;
	double bestFit;
	double best_x;
	double best_y;
};

struct GA2D_Control
{
	bool locked = false;
	uint32_t populationSize;
	uint32_t generationLimit;
	uint32_t mutation_severity; //Can be a value from 1-32
	double mutation_rate;	//TODO: Change naming to something more appropriate
};

typedef struct
{
	char *model_name;
	fitfunc_2d_t function;
	GA2D_Constraints limits;
	GA_BreedStyle_TypeDef breedStyle;
	GA_SelectStyle_TypeDef mateSelectStyle;
	GA_MutateDistribution_TypeDef mutationStyle;
} GA2D_FunctionBenchMark;

class GAClass_2D
{
public:
	GA2D_RuntimeData	algRuntimeData;

	void assignAlgorithm(GA2D_FunctionBenchMark benchMark);
	void localSearch(double x_center, double y_center, uint32_t axis_ticks, double *results, double percentage);

	GA_StatusTypeDef begin(uint32_t populationSize, uint32_t generationLimit, double mutationProbabilityThreshold);
	GA_StatusTypeDef reset();
	GA_StatusTypeDef iterateGeneration();
	GA_StatusTypeDef peekStatus(GA2D_RuntimeData *container);
	GA_StatusTypeDef peekStatus();

	
	GAClass_2D();
	~GAClass_2D();

private:
	/*-----------------------------
	* User interface data
	*-----------------------------*/
	fitfunc_2d_t fitnessFunction;
	GA_StatusTypeDef algStatus;
	GA2D_Constraints	algConstraints;
	GA2D_Control		algControl;
	GA_BreedStyle_TypeDef algBreedStyle;
	GA_SelectStyle_TypeDef algMateSelectStyle;
	GA_MutateDistribution_TypeDef algMutationStyle;

	std::random_device rnd_dev;
	uint32_t generation_counter = 0;
	double *currentFitness;
	double **currentPopulation;
	uint16_t **breedSelections;

	GA_StatusTypeDef initGeneration();
	GA_StatusTypeDef evaluateFitness();
	GA_StatusTypeDef selectMates(GA_SelectStyle_TypeDef selection_approach);
	GA_StatusTypeDef breedGeneration(GA_BreedStyle_TypeDef breeding_approach);
	GA_StatusTypeDef mutateChromosomes(GA_MutateDistribution_TypeDef mutation_approach);

	/*-----------------------------
	* Constants for Conversion Mapping
	*-----------------------------*/
	double bytes_precision;
	double x_offset, y_offset;
	double x_lo, y_lo;
	double x_hi, y_hi;
	double x_dPow, y_dPow;
	double x_bPow, y_bPow;
	double x_sF, y_sF;
	double x_sR, y_sR;

	/*-----------------------------
	* Helper Functions
	*-----------------------------*/
	void assignFitnessFunction(fitfunc_2d_t fitfunc);
	void assignSearchConstraint(GA2D_Constraints problemLimits);
	double uniformRandomNumber(double lower_bound, double upper_bound);
	uint16_t coordinate2Chromosome(double coordinate, uint32_t axis);
	double chromosome2Coordinate(uint16_t chromosome, uint32_t axis);
	double mutationProbability(GA_MutateDistribution_TypeDef distribution);
};
