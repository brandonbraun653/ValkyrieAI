#pragma once
#ifndef MODEL_H_
#define MODEL_H_

/* If there are build errors about winsock redefinitions, make sure
to in include the WIN32_LEAN_AND_MEAN macro in the project preprocessor
properties page. */
#if defined(_WIN32) && !defined(_WIN64)
#ifndef WIN32_LEAN_AND_MEAN
#error Please define "WIN32_LEAN_AND_MEAN" in the project preprocessor settings
#endif

#include <Winsock2.h>
#include <WS2tcpip.h>	//Needed for InetPton
#include <Windows.h>

#include <sys/types.h>
#include <tchar.h>		//Needed for _T
#include <string>

#pragma comment(lib, "Ws2_32.lib")
#define WINSOCKVERSION MAKEWORD(2,2 ) 
#endif

/* C/C++ Includes */
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <memory>

/* Boost Includes */
#include <boost/shared_ptr.hpp>
#include <boost/container/vector.hpp>
#include <boost/chrono.hpp>
#include <boost/thread.hpp>
#include <boost/filesystem.hpp>

/* JSON */
#include <nlohmann/json.hpp>
using json = nlohmann::json;

/* Local Includes */
#include "config.h"
#include "debugger.h"
#include "types.h"
#include "matlab.h"

struct SimCommand
{
	std::string axis;
	std::string simType;
	int numTimeSteps = 0;
	float dt = 0.0;
	float start_time = 0.0;
	float end_time = 0.0;
	float stepMagnitude = 0.0;
	float angle_kp = 0.0;
	float angle_ki = 0.0;
	float angle_kd = 0.0;
	float rate_kp = 0.0;
	float rate_ki = 0.0;
	float rate_kd = 0.0;
};

struct SimResults
{
	/* Results structure taken from here:
	* https://www.mathworks.com/help/control/ref/stepinfo.html
	**/

	float riseTime = 0.0;
	float settlingTime = 0.0;
	float settlingMin = 0.0;
	float settlingMax = 0.0;
	float overshoot = 0.0;
	float undershoot = 0.0;
	float peak = 0.0;
	float peakTime = 0.0;
};

#if defined(_WIN32) && !defined(_WIN64)
extern std::string parseStruct(SimCommand& input);
extern SimResults parseResults(std::string& results);
extern std::vector<std::string> splitString2String(const std::string& s, char delimiter);
extern boost::container::vector<double> splitString2Double(const std::string& str, char delimiter);


class NN_TCPModel : public NN_ModelBase
{
public:
	typedef struct
	{
		void* hFileMap;
		void* pData;
		char MapName[256];
		size_t Size;
		bool valid = false;
	} SharedMemory;

	NN_TCPModel(std::string host_ip, uint32_t port)
	{
		HOST = host_ip;
		PORT = port;
	}

	~NN_TCPModel()
	{
		closeConnection();

		if (simFile_Pitch)
			FreeMemoryMap(simFile_Pitch);
	}

	int initialize() override;
	void openConnection();
	void closeConnection();
	int send_data(std::string& input);
	std::string recv_data(int length = MAX_BUFFER);

	static const int MAX_BUFFER = 4096;

	SharedMemory* simFile_Pitch = nullptr;

private:
	bool connectionOpen = false;
	char buffer[MAX_BUFFER + 1];
	std::string HOST;
	uint32_t PORT;
	int connectionFd;

	struct sockaddr_in servAddr;
	struct sockaddr_in localAddr;


	bool CreateMemoryMap(SharedMemory* shm);
	bool FreeMemoryMap(SharedMemory* shm);
};
typedef boost::shared_ptr<NN_TCPModel> NN_TCPModel_sPtr;
#endif

class NN_JSONModel : public NN_ModelBase
{
public:

	int initialize() override;

	json executeModel(json cmd);


	NN_JSONModel(std::string working_directory)
	{
		root_dir = working_directory;
	}

	~NN_JSONModel() = default;

private:
	std::string root_dir = "";
	std::string cmd_file = "sim_cmd.json";
	std::string rcv_file = "sim_result.json";

	std::time_t last_write_time = 0;
};
typedef boost::shared_ptr<NN_JSONModel> NN_JSONModel_sPtr;



class MatlabModel : public MatlabSession, public ML_ModelBase
{
public:
	int initialize() override;

	void setModelRoot(const std::string path);
	void setModelFunction(const std::string path);
	void setInitFunction(const std::string path);

	matlab::data::StructArray blankInput();

	MatlabModel() = default;
	~MatlabModel() = default;

	/* Used to create data types for model I/O */
	matlab::data::ArrayFactory factory;
	std::string model_path = "";

	/* Simulation Parameters */
	double startTime = 0.0;
	double endTime = 5.0;
	double stepMagnitude = 10.0;
	double stepTimeEnable = 0.1;
	std::string simAxis = "";


private:
	std::string root_path = "";
	std::string init_path = "";

	
};
typedef boost::shared_ptr<MatlabModel> MatlabModel_sPtr;

#endif /* !MODEL_H_ */