#pragma once
#ifndef MODEL_H_
#define MODEL_H_

/* C/C++ Includes */
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <memory>

/* Eigen Includes */
#include <eigen/Eigen>
#include <eigen/StdVector>

/* Boost Includes */
#include <boost/thread.hpp>
#include <boost/chrono.hpp>
#include <boost/container/vector.hpp>
#include <boost/timer/timer.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen.hpp>

/* Local Includes */
#include "config.h"
#include "debugger.h"
#include "types.h"
#include "model_simulation.h"
#include "signal_analysis.h"

namespace odeint = boost::numeric::odeint;
void model_error_exit(std::string error_msg, int line = __LINE__, std::string file = __FILE__);
bool dim_assert(size_t row_act, size_t row_exp, size_t col_act, size_t col_exp);


#endif /* !MODEL_H_ */