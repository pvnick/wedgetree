#define DETECT_MEMORY_LEAKS

#ifdef DETECT_MEMORY_LEAKS
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <malloc.h>
#include <crtdbg.h>
#ifdef _DEBUG
#define new new(_CLIENT_BLOCK,__FILE__, __LINE__)
#endif  // _DEBUG
#endif


#include "node.h"
#include <chrono>
#include <iostream>
#include <fstream>
#include <vector>
#include <boost/program_options.hpp>
#include "wedgetree.h"
#include "cli_options.h"

using namespace std;

namespace program_options = boost::program_options;

void do_analysis() {
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
	size_t R = 10;
	WedgeTree tree(timeseries, M, B, r, R, verbose);
	tree.get_merged_candidate_nodes();
}

int main(int argc, char *argv[])
{
#ifdef DETECT_MEMORY_LEAKS
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
#endif
	CLIOptions::init(argc, argv);
	auto start = std::chrono::high_resolution_clock::now();
	do_analysis();
	auto end = std::chrono::high_resolution_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
	std::cout << "Analysis took " << elapsed.count() << " seconds" << std::endl;
	//unsorted 2131
	//sorted 1484 (30% faster)
	return 0;
}
