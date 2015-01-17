#include <vector>
#include <cmath>
#include <cassert>
#include <algorithm>
#include "candidate.h"
#include "boundaryposition.h"
#include "wedge.h"
#include "lb_kim.h"
#include "mathutils.h"

namespace LowerBound {
	double KimWedge::dist(double x, BoundaryPosition const& bound) {
		if (x > bound.Ui) {
			return MathUtils::fastdist(bound.Ui, x);
		}
		else if (x < bound.Li) {
			return MathUtils::fastdist(bound.Li, x);
		}
		return 0;
	}
	double KimWedge::lb(Wedge const& W, Candidate const& C, double abandon_after) {
		std::vector<BoundaryPosition> const& boundary_positions = W.get_boundary_positions();
		const size_t len = C.length;
		assert(boundary_positions.size() == len);

		/// 1 point at front and back
		double x0 = C.series_normalized[0];
		double y0 = C.series_normalized[len - 1];
		double result = dist(x0, boundary_positions[0]) + dist(y0, boundary_positions[len - 1]);
		//in ucr_dtw , the following check uses greater than *or equal to*, but we want this function to indicate that the candidate
		//cant possibly be similar to any of the candidates in the wedge, so we use greater than
		if (result > abandon_after) return result;

		/// 2 points at front
		double x1 = C.series_normalized[1];
		result += std::min({
			dist(x1, boundary_positions[0]),
			dist(x0, boundary_positions[1]),
			dist(x1, boundary_positions[1])
		});
		if (result > abandon_after) return result;

		/// 2 points at back
		double y1 = C.series_normalized[len - 2];
		result += std::min({
			dist(y1, boundary_positions[len - 1]),
			dist(y0, boundary_positions[len - 2]),
			dist(y1, boundary_positions[len - 2])
		});
		if (result > abandon_after) return result;

		/// 3 points at front
		double x2 = C.series_normalized[2];
		result += std::min({
			dist(x0, boundary_positions[2]),
			dist(x1, boundary_positions[2]),
			dist(x2, boundary_positions[2]),
			dist(x2, boundary_positions[1]),
			dist(x2, boundary_positions[0])
		});
		if (result > abandon_after) return result;

		/// 3 points at back
		double y2 = C.series_normalized[len - 3];
		result += std::min({
			dist(y0, boundary_positions[len - 3]),
			dist(y1, boundary_positions[len - 3]),
			dist(y2, boundary_positions[len - 3]),
			dist(y2, boundary_positions[len - 2]),
			dist(y2, boundary_positions[len - 1])
		});
		return result;
	}
}