#include "../include/matlab.h"

#ifdef _WIN64
#include "engine.h"

MatlabSession::~MatlabSession()
{

}

int MatlabSession::start()
{
	try
	{
		std::cout << "Starting Matlab Engine...";
		matlabPtr = std::move(matlab::engine::startMATLAB());
		std::cout << "Done" << std::endl;
	}
	catch (matlab::engine::EngineException err)
	{
		return 0;
	}
	return 1;
}

void MatlabSession::end()
{
	if (matlabPtr)
	{
		//terminateEngineClient();
	}
}

void MatlabSession::addPath(const std::string& path)
{
	std::string fCall = "addpath('" + path + "');";
	matlabPtr->eval(matlab::engine::convertUTF8StringToUTF16String(fCall));
}


#endif /* _WIN64 */