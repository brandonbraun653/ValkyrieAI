#pragma once
#ifndef MODEL_H_
#define MODEL_H_

/* If there are build errors about winsock redefinitions, make sure
to in include the WIN32_LEAN_AND_MEAN macro in the project preprocessor
properties page. */
#ifdef WIN32
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

/* Local Includes */
#include "config.h"
#include "debugger.h"
#include "types.h"

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

#endif /* !MODEL_H_ */