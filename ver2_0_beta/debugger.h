#pragma once
#ifndef DEBUGGER_H_
#define DEBUGGER_H_

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
#include <boost/container/vector.hpp>

/* Local Includes (Mainly just so functions can understand types) */
#include "data.h"

template <typename T>
void prettyPrint(const Eigen::EigenBase<T> matrix, std::string matrix_name);


void prettyPrint(boost::container::vector<double> data, std::string name);

template <typename T>
void writeCSV(T matrix, std::string filename)
{
	std::ofstream csvFile;
	csvFile.open(filename);
	size_t num_rows = matrix.rows();
	size_t num_cols = matrix.cols();

	for (int row = 0; row < num_rows; row++)
	{
		for (int col = 0; col < num_cols; col++)
			csvFile << std::to_string(matrix(row, col)*1.0) << ",";

		csvFile << "\n";
	}

	csvFile.close();
}



#endif