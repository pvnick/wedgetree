#include <vector>
#include <cmath>
#include <cassert>
#include <algorithm>
#include "candidate.h"
#include "boundaryposition.h"
#include "wedge.h"
#include "lb_keough.h"
#include "mathutils.h"

namespace LowerBound {
	double KeoughWedge::dist(double x, BoundaryPosition const& bound) {
		if (x > bound.DTW_Ui) {
			return MathUtils::fastdist(bound.DTW_Ui, x);
		}
		else if (x < bound.DTW_Li) {
			return MathUtils::fastdist(bound.DTW_Li, x);
		}
		return 0;
	}
	double KeoughWedge::lb(Wedge const& W, Candidate const& C, double abandon_after) {
		std::vector<BoundaryPosition> const& boundary_positions = W.get_boundary_positions();
		std::vector<size_t> const& bound_pos_indices_LB_keough_sorted = W.get_bound_pos_indices_LB_keough_sorted();
		const size_t len = C.length;
		double result = 0;
		assert(boundary_positions.size() == len);
		for (size_t i : bound_pos_indices_LB_keough_sorted) {
			result += dist(C.series_normalized[i], boundary_positions[i]);
			if (result > abandon_after)
				return result;
		}
		return result;
	};
}
