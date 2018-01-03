#include "ga_helper.h"

uint16_t data2Chromosome(mapCoeff_t *mapping, double data)
{
	return (uint16_t)(trunc(mapping->x_sF*(data - mapping->x_lo) + mapping->x_offset));
}

double chromosome2Data(mapCoeff_t *mapping, uint16_t chromosome)
{
	return (mapping->x_sR*(((float)chromosome) - mapping->x_offset) + mapping->x_lo);
}

double uniformRandomNumber(double lower_bound, double upper_bound)
{
	/* Create a new RNG using the Mersenne Twister engine. A new object must
	be created each time to guarantee a good level of non-repeating randomness.*/
	std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_real_distribution<double> uniform_number(lower_bound, upper_bound);
	return uniform_number(rng);
}

double enforceResolution(double in, GA_Resolution res)
{
	double fracPart = 0.0;
	double intPart = 0.0;

	/* Decompose the data into integral and fractional parts */
	fracPart = std::modf(in, &intPart);

	/* Shift up, truncate, shift down */
	fracPart *= std::pow(10.0, (int)res);
	fracPart = floor(fracPart);
	fracPart /= std::pow(10.0, (int)res);

	return (intPart + fracPart);
}
