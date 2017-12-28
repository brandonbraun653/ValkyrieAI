#include "logger.h"

namespace logging = boost::log;
namespace src = boost::log::sources;
namespace sinks = boost::log::sinks;
namespace keywords = boost::log::keywords;
namespace expr = boost::log::expressions;
namespace attrs = boost::log::attributes;

/* NOTE:
1. Interestingly enough, the name_ parameter in BOOST_LOG_ATTRIBUTE_KEYWORD
   is required to follow a specific naming scheme or else logs won't even 
   be written! "Severity" must have a capital S. "Tag" must have a capital T.
   Why? Great question...I have no clue....TODO
*/
BOOST_LOG_ATTRIBUTE_KEYWORD(tag_attr, "Tag", std::string)
BOOST_LOG_ATTRIBUTE_KEYWORD(sev_attr, "Severity", severity_level)


boost::log::sources::severity_logger<severity_level> logger;


// The operator puts a human-friendly representation of the severity level to the stream
/* TODO: How does this even work? I never call "operator<<"???....*sigh. This c++ stuff
is complex...*/
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

void add_custom_attributes()
{
	boost::shared_ptr< logging::core > core = logging::core::get();
	core->add_global_attribute("FileName", attrs::constant<std::string>("FileName"));
}

void init_boost_logger(const char *filename)
{
	/* Generic formatter for the output */
	logging::formatter fmt = expr::stream
		<< expr::format_date_time< boost::posix_time::ptime >("TimeStamp", "%Y-%m-%d %H:%M:%S")
		<< ":<" << sev_attr << "> "
		<< expr::if_(expr::has_attr(tag_attr))
		[
			expr::stream << "[" << tag_attr << "] "
		]
		<< expr::smessage;


	boost::shared_ptr< sinks::text_multifile_backend > backend =
		boost::make_shared< sinks::text_multifile_backend >();

	/* Setup the file naming pattern */
	backend->set_file_name_composer(sinks::file::as_file_name_composer(expr::stream << filename));

	/* Wrap the pattern into the front end and register to the core. Without
	synchronizing both, the sink doesn't work.*/
	typedef sinks::synchronous_sink< sinks::text_multifile_backend > sink_t;
	boost::shared_ptr< sink_t > sink(new sink_t(backend));

	/* Set the formatter */
	sink->set_filter
	(
		(sev_attr >= normal) ||
		(expr::has_attr(tag_attr))
	);
	sink->set_formatter(fmt);

	/* Add the sink to the core */
	logging::core::get()->add_sink(sink);
	logging::add_common_attributes();

	BOOST_LOG_SEV(logger, normal) << "----Log Opened----";
}

