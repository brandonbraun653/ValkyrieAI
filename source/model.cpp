#include "model.h"
 
void model_error_exit(std::string error_msg, int line, std::string file)
{
	std::cout << "YOU FOOL! You've caused an error!!" << std::endl;
	std::cout << error_msg << std::endl;
	std::cout << "Error at line " << std::to_string(line) << " in file " << file << std::endl;
	exit(EXIT_FAILURE);
}

bool dim_assert(size_t row_act, size_t row_exp, size_t col_act, size_t col_exp)
{
	if ((row_act != row_exp) || (col_act != col_exp))
		return false;

	return true;
}

template <typename T>
void writeCSV(T matrix, StepPerformance data, std::string filename)
{
	std::ofstream csvFile;
	csvFile.open(filename);
	size_t num_rows = matrix.rows();
	size_t num_cols = matrix.cols();

	//Write each full row
	for (int row = 0; row < num_rows; row++)
	{
		for (int col = 0; col < num_cols; col++)
			if (col == num_cols - 1)
				csvFile << std::to_string(matrix(row, col)*1.0);
			else
				csvFile << std::to_string(matrix(row, col)*1.0) << ",";

		csvFile << "\n";
	}

	//Write the performance data
	csvFile << data.steadyStateValue_performance << "," << data.settlingPcntRange << "\n";
	csvFile << data.delta_overshoot_performance << "\n";
	csvFile << data.percentOvershoot_performance << "\n";
	csvFile << data.riseTime_performance << "," << data.riseTime_Idx[0] << "," << data.riseTime_Idx[1] << "\n";
	csvFile << data.settlingTime_performance << "," << data.settlingTime_Idx << "\n";
	csvFile << data.steadyStateError_performance << "\n";

	csvFile.close();
}

void writeCSV_StepPerfData(StepPerformance data, std::string filename)
{
	std::ofstream csvFile;
	csvFile.open(filename);
	size_t num_rows = data.data.rows();
	size_t num_cols = data.data.cols();

	//Write each full row
	for (int row = 0; row < num_rows; row++)
	{
		for (int col = 0; col < num_cols; col++)

			if (col == num_cols - 1)
				csvFile << std::to_string(data.data(row, col)*1.0);
			else
				csvFile << std::to_string(data.data(row, col)*1.0) << ",";

		csvFile << "\n";
	}

	//Write the performance data
	csvFile << data.steadyStateValue_performance << "," << data.settlingPcntRange << "\n";
	csvFile << data.delta_overshoot_performance << "\n";
	csvFile << data.percentOvershoot_performance << "\n";
	csvFile << data.riseTime_performance << "," << data.riseTime_Idx[0] << "," << data.riseTime_Idx[1] << "\n";
	csvFile << data.settlingTime_performance << "," << data.settlingTime_Idx << "\n";
	csvFile << data.steadyStateError_performance << "\n";

	csvFile.close();
}


// StepPerformance StateSpaceModel::stepResponseSingleThreaded(int member_num, SS_NLTIVModel &model, double Kp, double Ki, double Kd)
// {
// 	StepPerformance temp;
// 	Eigen::MatrixXd raw_data = step_simulator.stepResponse(model, Kp, Ki, Kd);
// 
// 	/* Check the extreme values in the simulation */
// 	double maxValue = raw_data.maxCoeff();
// 	double minValue = raw_data.minCoeff();
// 	double threshold = 100.0;
// 
// 	/* Assuming we didn't go to crazy values, analyze the data returned */
// 	if ((maxValue < threshold) && (minValue > -threshold))
// 	{
// 		temp = step_analyzer.analyze(raw_data);
// 		temp.performance_simulation_data_is_valid = true;
// 	}
// 	else
// 		temp.performance_simulation_data_is_valid = false;
// 
// 	/* Log the settings simulated under */
// 	temp.Kp = Kp;
// 	temp.Ki = Ki;
// 	temp.Kd = Kd;
// 	temp.performance_simulation_data = raw_data;
// 
// 	return temp;
// 
// 	#ifdef SS_TRACE_LOG_CSV
// 	std::string filename = "member_data_" + std::to_string(member_num) + ".csv";
// 	//writeCSV<matrixXd>(raw_data, step_performance, filename);
// 	#endif
// 
// 	#ifdef SS_TRACE_LOG_COUT
// 	std::cout << "Model evaluated!" << std::endl;
// 	#endif
//  
// }
// 
// void StateSpaceModel::stepResponseMultiThreaded(int member_num, SS_NLTIVModel model, StepPerformance_Vec& output_data, boost::mutex& output_data_mutex,
// 	double Kp, double Ki, double Kd)
// {
// 	#ifdef SS_TRACE
// 	#ifdef SS_THREAD_TRACE_LOG_COUT
// 	cout_mutex.lock();
// 	std::cout << "Start: TID " << boost::this_thread::get_id() << std::endl;
// 	cout_mutex.unlock();
// 	#endif
// 	auto start = boost::chrono::high_resolution_clock::now();
// 	#endif
// 
// 	/*-----------------------------------------------
// 	* Setup local thread objects
// 	*-----------------------------------------------*/
// 	SS_NLTIVModel thread_model(model);			//Must be a deep copy?
// 	StateSpaceSimulator thread_step_simulator;
// 	StepResponseAnalyzer thread_step_analyzer;
// 
// 	/*-----------------------------------------------
// 	* Analyze
// 	*-----------------------------------------------*/
// 	StepPerformance temp;
// 	Eigen::MatrixXd thread_raw_data = thread_step_simulator.stepResponse(model, Kp, Ki, Kd);
// 	
// 	/* Check the extreme values in the simulation */
// 	double maxValue = thread_raw_data.maxCoeff();
// 	double minValue = thread_raw_data.minCoeff();
// 	double threshold = 100.0;
// 
// 	/* Assuming we didn't go to crazy values, analyze the data returned */
// 	if ((maxValue < threshold) && (minValue > -threshold))
// 	{
// 		temp = thread_step_analyzer.analyze(thread_raw_data);
// 		temp.performance_simulation_data_is_valid = true;
// 	}
// 
// 	/* Log the settings simulated under */
// 	temp.Kp = Kp;
// 	temp.Ki = Ki;
// 	temp.Kd = Kd;
// 	temp.performance_simulation_data = thread_raw_data;
// 
// 	/* Return the data */
// 	output_data_mutex.lock();
// 	output_data[member_num] = temp;
// 	output_data_mutex.unlock();
// 
// 	#ifdef SS_TRACE
// 	#ifdef SS_THREAD_TRACE_LOG_COUT
// 	auto end = boost::chrono::high_resolution_clock::now();
// 
// 	cout_mutex.lock();
// 	std::cout << "End  : TID "
// 		<< boost::this_thread::get_id() << "  Time Taken: "
// 		<< (end - start).count() * ((double)boost::chrono::high_resolution_clock::period::num / boost::chrono::high_resolution_clock::period::den)
// 		<< std::endl;
// 	cout_mutex.unlock();
// 	#endif
// 
// 	#ifdef SS_THREAD_TRACE_LOG_CSV
// 	std::string filename = "member_data_" + std::to_string(member_num) + ".csv";
// 
// 	csv_mutex.lock();
// 	writeCSV_StepPerfData(*output_data_address, filename);
// 	csv_mutex.unlock();
// 
// 	#endif
// #endif
// }
