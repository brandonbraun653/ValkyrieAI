#pragma once

#ifndef AI_MATLAB_H
#define AI_MATLAB_H
/* In order to use the Matlab interface, several steps must be taken. This assumes x64 build.

1. Add these paths to the environment PATH variable: (Restart computer for changes to take effect)
	<MATLAB_ROOT>\extern\lib\win64\microsoft
	<MATLAB_ROOT>\extern\bin\win64
	<MATLAB_ROOT>\bin\win64

2. Add these paths to the "Additional Include Directories" project property:
	<MATLAB_ROOT>\extern\include

3. Add paths from (1) to the "Additional Library Directories" project property:
	<MATLAB_ROOT>\extern\lib\win64\microsoft
	
*/
#ifdef _WIN64
//Must be in PATH variable, per(1)
#pragma comment (lib, "libmx.lib")
#pragma comment (lib, "libeng.lib")
#pragma comment (lib, "libMatlabDataArray.lib")
#pragma comment (lib, "libMatlabEngine.lib")

#include <stdlib.h>
#include <string>
#include <iostream>


#include "MatlabEngine.hpp"
#include "MatlabDataArray.hpp"

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>


class MatlabSession
{
public:
	int start();
	void end();

	void addPath(const std::string& path);

	boost::shared_ptr<matlab::engine::MATLABEngine> matlabPtr;

	MatlabSession() = default;
	~MatlabSession();

private:

};
#endif /* _WIN64 */

#endif 