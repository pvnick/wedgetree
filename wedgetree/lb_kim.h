#ifndef HEADER_D90ECF27CACF41FCAD6422205B32C42A
#define HEADER_D90ECF27CACF41FCAD6422205B32C42A

#include <vector>
#include <cmath>
#include "candidate.h"
#include "boundaryposition.h"

namespace LowerBound {
	class KimWedge {
	private:
		static double dist(BoundaryPosition const& bound, double x);
	public:
		static double lb(std::vector<BoundaryPosition> const& boundary_positions, Candidate const& C, double abandon_after);
	};
}

#endif
