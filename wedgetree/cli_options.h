#ifndef HEADER_42D591DFA61520B
#define HEADER_42D591DFA61520B


#include <stdexcept>
#include <boost/program_options.hpp>

namespace program_options = boost::program_options;

class CLIOptions : public program_options::variables_map {
private:
	static CLIOptions instance;
	CLIOptions() = default;
public:
	static bool initialized;
	static void init(int argc, char *argv[]) {
		if (!initialized) {
			std::string cmd;
			// Declare the supported options.
			program_options::options_description desc("Allowed options");
			desc.add_options()
				("verbose", "output debug info")
				("K",
				program_options::value<unsigned int>()->default_value(100)->value_name("INT"),
				"set the number of top candidates to store per query")
				("M",
				program_options::value<unsigned int>()->default_value(100)->value_name("INT"),
				"subsequence/motif length")
				("timeseries",
				program_options::value<std::string>()->required()->value_name("FILE"),
				"input timeseries file path (one value per line)")
				;

			try {
				program_options::store(program_options::command_line_parser(argc, argv)
					.options(desc)
					.run(),
					instance);

				program_options::notify(instance);
			}
			catch (program_options::error& e) {
				std::cout << desc << std::endl;
				exit(0);
			}
			initialized = true;
		}
	}
	static CLIOptions const& get_instance() {
		if (initialized)
			return instance;
		throw std::logic_error("You must initialize the options class before retrieving an instance");
	}
	bool verbose() {
		return this->count("verbose");
	}
};


#endif // header guard
