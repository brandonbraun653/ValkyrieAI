#include "debugger.h"

template <typename T>
void prettyPrint(const Eigen::EigenBase<T> matrix, std::string matrix_name)
{
	size_t num_col = matrix.cols();
	size_t num_row = matrix.rows();
	double data = 0;
	std::cout << std::fixed << std::setprecision(3);
	std::cout << "Matrix: " << matrix_name << std::endl;

	// PRINT COLUMN HEADER
	std::cout << "\t";
	for (int col = 0; col < num_col; col++)
		std::cout << "[" << col << "]\t";
	std::cout << std::endl;

	// PRINT ROWS 
	for (int row = 0; row < num_row; row++)
	{
		std::cout << "[" << row << "]" << "\t";
		for (int col = 0; col < num_col; col++)
		{

			std::cout << matrix(row, col)*1.0 << "\t";

		}
		std::cout << std::endl;
	}
}


void prettyPrint(boost::container::vector<double> data, std::string name)
{
	std::cout << name << std::endl;

	for (int i = 0; i < data.size(); i++)
	{


		std::cout << data[i] << ", ";

	}
	std::cout << std::endl;
}

