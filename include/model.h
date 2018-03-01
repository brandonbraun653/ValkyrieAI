#pragma once
#ifndef MODEL_H_
#define MODEL_H_
/* Winsock2 Includes */
#ifdef WIN32
#include <WS2tcpip.h>	//Needed for InetPton
#include <Windows.h>

#include <sys/types.h>
#include <Winsock2.h>
#include <tchar.h>		//Needed for _T
#include <string>
#include <iostream>

#pragma comment(lib, "Ws2_32.lib")
#define WINSOCKVERSION MAKEWORD(2,2 ) 
#endif

struct sockaddr_in servAddr;
struct sockaddr_in localAddr;


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


class NN_TCPModel : public NN_ModelBase
{
public:
	/* Copy Constructor */
	NN_TCPModel(const NN_ModelBase_sPtr& base)
	{

	}

	NN_TCPModel(std::string host_ip, uint32_t port, uint32_t bufferSize)
	{

	}


	void setup() override;
	void openConnection();
	void closeConnection();




private:
	std::string host;
	uint32_t port;

};
typedef boost::shared_ptr<NN_TCPModel> NN_TCPModel_sPtr;

#endif /* !MODEL_H_ */