#ifndef LOGGER_H_
#define LOGGER_H_

/* C/C++ Includes */
#include <stdio.h>
#include <iomanip>
#include <stdlib.h>
#include <string>
#include <iostream>

/* Boost Includes */
#include <boost/smart_ptr/shared_ptr.hpp>
#include <boost/log/common.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/attributes.hpp>
#include <boost/log/sinks/sync_frontend.hpp>
#include <boost/log/sinks/text_multifile_backend.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/support/date_time.hpp>

enum severity_level
{
	normal,
	notification,
	warning,
	error,
	critical
};

extern boost::log::sources::severity_logger<severity_level> logger;

extern void init_boost_logger(const char *filename);


#endif /* !LOGGER_H_ */	