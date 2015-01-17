#ifndef HEADER_EBBA17D58C1D41C998C62E20DAB2D1B7
#define HEADER_EBBA17D58C1D41C998C62E20DAB2D1B7

#include <limits>
#include <vector>

struct BoundaryPosition {
	size_t index;
	double ED = 0;
	double Ui = std::numeric_limits<double>::min();
	double Li = std::numeric_limits<double>::max();
	double DTW_Ui = std::numeric_limits<double>::min();
	double DTW_Li = std::numeric_limits<double>::max();
	struct Compare_UL_ED {
		//for sorting the upper and lower wedge boundary positions in order of increasing euclidean distance
		std::vector<BoundaryPosition> const& boundary_positions;
		bool operator()(size_t const& c_index1, size_t const& c_index2) const {
			return boundary_positions[c_index1].ED < boundary_positions[c_index2].ED;
		}
		Compare_UL_ED(std::vector<BoundaryPosition> const& boundary_positions) : boundary_positions(boundary_positions) { }
	};
	struct Compare_DTW_UL {
		//for sorting the upper and lower lb_keough envelope in order of increasing euclidean distance
		std::vector<BoundaryPosition> const& boundary_positions;
		bool operator()(size_t const& c_index1, size_t const& c_index2) const {
			return (boundary_positions[c_index1].DTW_Ui - boundary_positions[c_index1].DTW_Li) < (boundary_positions[c_index2].DTW_Ui - boundary_positions[c_index2].DTW_Li);
		}
		Compare_DTW_UL(std::vector<BoundaryPosition> const& boundary_positions) : boundary_positions(boundary_positions) { }
	};
};

#endif