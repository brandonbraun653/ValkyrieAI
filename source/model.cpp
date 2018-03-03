#include "model.h"
 
const int CONNECT_RETRY_MAX = 100;

std::string parseStruct(SimCommand& input)
{
	//The "pack" and "pop" is only for helping the python side out
	//so that any extra characters appended on to the first and last 
	//values can be easily ignored. This is pure laziness at its finest.
	std::string output =
		input.axis + ',' + \
		input.simType + ',' + \
		std::to_string(input.numTimeSteps) + ',' + \
		std::to_string(input.dt) + ',' + \
		std::to_string(input.start_time) + ',' + \
		std::to_string(input.end_time) + ',' + \
		std::to_string(input.stepMagnitude) + ',' + \
		std::to_string(input.angle_kp) + ',' + \
		std::to_string(input.angle_ki) + ',' + \
		std::to_string(input.angle_kd) + ',' + \
		std::to_string(input.rate_kp) + ',' + \
		std::to_string(input.rate_ki) + ',' + \
		std::to_string(input.rate_kd) + ',' + \
		"pop" + '\n';

	return output;
}

SimResults parseResults(std::string& results)
{
	SimResults output;

	/* Finds the end of the return data using '\n' as the delimiter*/
	std::string data = results.substr(0, results.find("\n"));

	std::vector<std::string> values = splitString(data, ',');

	output.riseTime = stof(values[0]);
	output.settlingTime = stof(values[1]);
	output.settlingMin = stof(values[2]);
	output.settlingMax = stof(values[3]);
	output.overshoot = stof(values[4]);
	output.undershoot = stof(values[5]);
	output.peak = stof(values[6]);
	output.peakTime = stof(values[7]);

	return output;
}

std::vector<std::string> splitString(const std::string& s, char delimiter)
{
	std::stringstream ss(s);
	std::string item;
	std::vector<std::string> tokens;

	while (std::getline(ss, item, delimiter))
		tokens.push_back(item);

	return tokens;
}


int NN_TCPModel::initialize()
{
	WSADATA wsaData;
	if (WSAStartup(WINSOCKVERSION, &wsaData) != 0)
		return ERROR;

	memset(&servAddr, 0, sizeof(servAddr));
	servAddr.sin_family = AF_INET;
	servAddr.sin_port = htons(PORT);

	/* Conversion of the HOST ip into the weird LPCWSTR type:
	https://stackoverflow.com/questions/27220/how-to-convert-stdstring-to-lpcwstr-in-c-unicode */
	std::wstring stemp = std::wstring(HOST.begin(), HOST.end());
	LPCWSTR sw = stemp.c_str();

	InetPton(AF_INET, sw, &servAddr.sin_addr.s_addr);

	// Create socket
	connectionFd = socket(AF_INET, SOCK_STREAM, 0);

	/* bind any port number */
	localAddr.sin_family = AF_INET;
	localAddr.sin_addr.s_addr = htonl(INADDR_ANY);
	localAddr.sin_port = htons(0);

	bind(connectionFd, (struct sockaddr *)&localAddr, sizeof(localAddr));

	return 1;
}

void NN_TCPModel::openConnection()
{
	int connectTries = 0;
	while (connectTries < CONNECT_RETRY_MAX)
	{
		std::cout << "Attempt " << connectTries << " of " << CONNECT_RETRY_MAX << std::endl;
		if (connect(connectionFd, (struct sockaddr *)&servAddr, sizeof(servAddr)) == 0)
		{
			connectionOpen = true;
			std::cout << "Connected!" << std::endl;
			break;
		}

		std::cout << "Couldn't connect. Will retry..." << std::endl;
		Sleep(2000);
		++connectTries;

		if (connectTries == CONNECT_RETRY_MAX)
		{
			throw std::runtime_error("Please start the python NN TCP server so we have something to talk to.");
			exit(1);
		}	
	}
}

void NN_TCPModel::closeConnection()
{
	if (connectionOpen)
		closesocket(connectionFd);
}

int NN_TCPModel::send_data(std::string& input)
{
	sprintf_s(buffer, "%s", input.data());
	return send(connectionFd, buffer, strlen(buffer), 0);
}

std::string NN_TCPModel::recv_data(int length)
{
	sprintf_s(buffer, "%s", "");
	recv(connectionFd, buffer, MAX_BUFFER, 0);
	std::string recv_buffer = buffer;

	return recv_buffer;
}