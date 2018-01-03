#pragma once
#ifndef HELPER_H_
#define HELPER_H_

/* C/C++ Includes */
#include <stdlib.h>
#include <stdint.h>
#include <iostream>
#include <random>

/* Local Includes */
#include "ga_config.h"
#include "data.h"

extern uint16_t data2Chromosome(mapCoeff_t *mapping, double data);
extern double chromosome2Data(mapCoeff_t *mapping, uint16_t chromosome);
extern double uniformRandomNumber(double lower_bound, double upper_bound);
extern double enforceResolution(double in, GA_Resolution res);
#endif