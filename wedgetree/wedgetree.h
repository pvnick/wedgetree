#ifndef HEADER_C0C1BBB0C1214EB89842E70002FBDFCB
#define HEADER_C0C1BBB0C1214EB89842E70002FBDFCB

#include <iostream>
#include <vector>
#include <stdexcept>
#include <limits>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <memory>
#include "candidate.h"
#include "node.h"

class WedgeTree {
private:
	std::vector<double> const& timeseries;
	size_t M;
	size_t B;
	double r;
	size_t R;
	std::unique_ptr<InternalWedgeNode> root;
	bool verbose;
public:
	WedgeTree(std::vector<double> const& timeseries, size_t M, size_t B, double r, size_t R, bool verbose);
	void insert_timeseries(size_t candidate_position);
	std::list<CandidateNode> get_merged_candidate_nodes() const;
};

#endif
