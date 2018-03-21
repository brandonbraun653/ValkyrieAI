#pragma once
#ifndef USER_OPTIMIZER_CONFIG_H
#define USER_OPTIMIZER_CONFIG_H

/* The sole purpose of this header is to provide functions that return initialization 
settings for various tuners. */


#include "valkyrie_engine.h"

extern FCSOptimizer_Init_t rollTunerVer1_Init();
extern FCSOptimizer_Init_t pitchTunerVer1_Init();
extern FCSOptimizer_Init_t yawTunerVer1_Init();

extern FCSOptimizer_Init_t rollTunerVer2_Init();

extern FCSOptimizer_Init_t dynamicInit();
extern FCSOptimizer_Init_t staticInit1();
extern FCSOptimizer_Init_t staticInit2();
extern FCSOptimizer_Init_t staticInit3();
#endif