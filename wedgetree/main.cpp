#include "node.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <boost/program_options.hpp>
#include "wedgetree.h"
#include "cli_options.h"

using namespace std;

namespace program_options = boost::program_options;

int main(int argc, char *argv[])
{
	CLIOptions::init(argc, argv);
	CLIOptions opts = CLIOptions::get_instance();

	std::string const& ts_filepath = opts["timeseries"].as<std::string>();
	std::ifstream ts_in(ts_filepath);
	if (ts_in.bad()) {
		throw std::invalid_argument("Timeseries file not found");
	}
	std::vector<double> timeseries;
	double val;
	while (ts_in >> val)
		timeseries.push_back(val);

	unsigned int K = opts["K"].as<unsigned int>();
	unsigned int M = opts["M"].as<unsigned int>();
	bool verbose = opts.verbose();

	size_t B = 5;
	double r = 300;
	WedgeTree(timeseries, M, B, r, verbose);
	return 0;
}
