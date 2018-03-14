#ifndef LOGGER_H_
#define LOGGER_H_

/* C/C++ Includes */
#include <stdlib.h>
#include <stdint.h>
#include <string>
#include <iostream>

/* Boost Includes */
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/log/common.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/attributes.hpp>
#include <boost/log/sinks/sync_frontend.hpp>
#include <boost/log/sinks/text_multifile_backend.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/support/date_time.hpp>

enum severity_level
{
	normal,
	notification,
	warning,
	error,
	critical
};

class SimpleLogger
{
public:
	void logInit(std::string filename);

	void logNormal(std::string msg);
	void logNotification(std::string msg);
	void logWarning(std::string msg);
	void logError(std::string msg);
	void logCritical(std::string msg);
	
	SimpleLogger();
	~SimpleLogger() = default;

private:
	typedef boost::log::sources::severity_channel_logger<severity_level, std::string> logger_type;
	logger_type* logger;

	static const int idLen = 32;
	std::string channel;
};


#endif /* !LOGGER_H_ */	