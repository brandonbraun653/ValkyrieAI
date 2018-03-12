#include "model.h"
 


#if defined(_WIN32) && !defined(_WIN64)
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

	std::vector<std::string> values = splitString2String(data, ',');

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

// Eventually, build a template version like this:
// https://gist.github.com/mark-d-holmberg/862733 

std::vector<std::string> splitString2String(const std::string& s, char delimiter)
{
	std::stringstream ss(s);
	std::string item;
	std::vector<std::string> tokens;

	while (std::getline(ss, item, delimiter))
		tokens.push_back(item);

	return tokens;
}

boost::container::vector<double> splitString2Double(const std::string& str, char delimiter)
{
	std::stringstream ss(str);
	std::string item;
	double value;
	boost::container::vector<double> tokens;
	tokens.resize(str.size());

	int index = 0;
	while (std::getline(ss, item, delimiter))
	{
		value = std::stod(item);
		tokens[index] = value;
		++index;
	}


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

	
	InetPton(AF_INET, HOST.data(), &servAddr.sin_addr.s_addr);

	// Create socket
	connectionFd = socket(AF_INET, SOCK_STREAM, 0);

	/* Bind any port number */
	localAddr.sin_family = AF_INET;
	localAddr.sin_addr.s_addr = htonl(INADDR_ANY);
	localAddr.sin_port = htons(0);

	bind(connectionFd, (struct sockaddr *)&localAddr, sizeof(localAddr));

	/*-----------------------------
	* Setup the memory mapped file 
	*-----------------------------*/
	simFile_Pitch = new SharedMemory;

	simFile_Pitch->Size = 256 * 1024 * sizeof(char);
	sprintf_s(simFile_Pitch->MapName, "Local\\PitchData");

	if (CreateMemoryMap(simFile_Pitch))
	{
		simFile_Pitch->valid = true;
		memset(simFile_Pitch->pData, 0, simFile_Pitch->Size);
	}

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

bool NN_TCPModel::CreateMemoryMap(SharedMemory* shm)
{
	if ((shm->hFileMap = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, shm->Size, shm->MapName)) == NULL)
		return false;

	if ((shm->pData = MapViewOfFile(shm->hFileMap, FILE_MAP_ALL_ACCESS, 0, 0, shm->Size)) == NULL)
	{
		CloseHandle(shm->hFileMap);
		return false;
	}
	return true;
}

bool NN_TCPModel::FreeMemoryMap(SharedMemory* shm)
{
	if (shm && shm->hFileMap)
	{
		if (shm->pData)
			UnmapViewOfFile(shm->pData);

		if (shm->hFileMap)
			CloseHandle(shm->hFileMap);

		return true;
	}
	return false;
}
#endif

int MatlabModel::initialize()
{
	if (!this->start())
		return 0;

	if (!root_path.empty())
		this->addPath(root_path);

	if (!init_path.empty())
	{
		/* I don't need a configuration option, but the Matlab API doesn't let you call
		functions with zero arguments. */
		auto configOption = factory.createScalar<int>(0);
		auto result = matlabPtr->feval(matlab::engine::convertUTF8StringToUTF16String(init_path), configOption);

		return (int)result[0];
	}
	else
		return 0;
}

void MatlabModel::setModelRoot(const std::string path)
{
	root_path = path;
}

void MatlabModel::setModelFunction(const std::string path)
{
	model_path = path;
}

void MatlabModel::setInitFunction(const std::string path)
{
	init_path = path;
}