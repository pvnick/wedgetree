#ifndef HEADER_EBBA17D58C1D41C998C62E20DAB2D1B7
#define HEADER_EBBA17D58C1D41C998C62E20DAB2D1B7

#include <limits>
#include <vector>

struct BoundaryPosition {
	size_t index;
	double ED = 0;
	double Ui = std::numeric_limits<double>::min();
	double Li = std::numeric_limits<double>::max();
	struct Compare {
		std::vector<BoundaryPosition> const& boundary_positions;
		bool operator()(size_t const& c_index1, size_t const& c_index2) const {
			return boundary_positions[c_index1].ED < boundary_positions[c_index2].ED;
		}
		Compare(std::vector<BoundaryPosition> const& boundary_positions) : boundary_positions(boundary_positions) { }
	};
};

#endif