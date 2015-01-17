#ifndef HEADER_D90ECF27CACF41FCAD6422205B32C42A
#define HEADER_D90ECF27CACF41FCAD6422205B32C42A

#include <vector>
#include <cmath>
#include "candidate.h"
#include "boundaryposition.h"
#include "wedge.h"

namespace LowerBound {
	class KeoughWedge { 
	//tight O(N) lower bound
	//lower bounds LB_Keough between C and every candidate held in the wedge
	private:
		static double dist(double x, BoundaryPosition const& bound);
	public:
		static double lb(Wedge const& W, Candidate const& C, double abandon_after);
	};
}

#endif
