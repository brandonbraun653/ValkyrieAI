#pragma once
#ifndef HELPER_H_
#define HELPER_H_

/* C/C++ Includes */
#include <stdlib.h>
#include <stdint.h>
#include <iostream>
#include <random>

/* Local Includes */
#include "config.h"
#include "types.h"

extern uint16_t data2Chromosome(FCSOptimizer_MappingCoeff *mapping, double data);
extern double chromosome2Data(FCSOptimizer_MappingCoeff *mapping, uint16_t chromosome);
extern double uniformRandomNumber(double lower_bound, double upper_bound);

#endif