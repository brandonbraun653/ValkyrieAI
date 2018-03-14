#include "logger.h"

namespace logging = boost::log;
namespace src = boost::log::sources;
namespace sinks = boost::log::sinks;
namespace keywords = boost::log::keywords;
namespace expr = boost::log::expressions;
namespace attrs = boost::log::attributes;

BOOST_LOG_ATTRIBUTE_KEYWORD(a_timestamp, "TimeStamp", boost::posix_time::ptime)
BOOST_LOG_ATTRIBUTE_KEYWORD(a_channel, "Channel", std::string)
BOOST_LOG_ATTRIBUTE_KEYWORD(tag_attr, "Tag", std::string)
BOOST_LOG_ATTRIBUTE_KEYWORD(sev_attr, "Severity", severity_level)

//Prints a user friendly version of the severity level
std::ostream& operator<< (std::ostream& strm, severity_level level)
{
	static const char* strings[] =
	{
		"normal",
		"notification",
		"warning",
		"error",
		"critical"
	};

	if (static_cast<std::size_t>(level) < sizeof(strings) / sizeof(*strings))
		strm << strings[level];
	else
		strm << static_cast<int>(level);

	return strm;
}


SimpleLogger::SimpleLogger()
{
	/* Generate a unique channel ID for the log file */
	static const char alphanum[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
	
	for (int i = 0; i < idLen; i++)
		channel.push_back(alphanum[rand() % (sizeof(alphanum) - 1)]);

	logger = new logger_type(keywords::channel = channel);
}

void SimpleLogger::logInit(std::string filename)
{
	if (filename.empty())
		throw std::runtime_error("Please specify a logging filename");

	logging::add_common_attributes();

	logging::add_file_log(
		keywords::file_name = filename,
		keywords::filter = (a_channel == channel),
		keywords::format = (expr::stream 
			//<< "[" << expr::format_date_time<boost::posix_time::ptime>("TimeStamp", "%Y-%m-%d %H:%M:%S") << "]"
			<< "<" << sev_attr << ">\t\t\t" << expr::smessage));
	
	BOOST_LOG_SEV(*logger, normal) << "----Log Opened----";
}

void SimpleLogger::logNormal(std::string msg)
{
	BOOST_LOG_SEV(*logger, normal) << msg;
}

void SimpleLogger::logNotification(std::string msg)
{
	BOOST_LOG_SEV(*logger, notification) << msg;
}

void SimpleLogger::logWarning(std::string msg)
{
	BOOST_LOG_SEV(*logger, warning) << msg;
}

void SimpleLogger::logError(std::string msg)
{
	BOOST_LOG_SEV(*logger, error) << msg;
}

void SimpleLogger::logCritical(std::string msg)
{
	BOOST_LOG_SEV(*logger, critical) << msg;
}
