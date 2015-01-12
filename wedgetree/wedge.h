#ifndef HEADER_A3E13BE6A09A41989CB88235F7155A45
#define HEADER_A3E13BE6A09A41989CB88235F7155A45

#include <memory>
#include <vector>
#include <boost/heap/d_ary_heap.hpp>
#include "candidate.h"

class Wedge {
	//every node stores a wedge, which tightly bounds all of its descendent candidates
private:
	struct BoundaryPosition {
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
	std::vector<double> const& timeseries;
	size_t M;
	size_t B;
	double r;
	double ED = 0; //squared euclidean distance between upper and lower bounds
	bool enlarged = false;
	static size_t id_counter; //holds the number of instantiated wedges
	size_t id; //set to the value of id_counter when this wedge was created
	std::vector<size_t> bound_pos_indices_sorted;
	std::vector<BoundaryPosition> boundary_positions;
	BoundaryPosition::Compare bound_comparator;
public:
	Wedge(std::vector<double> const& timeseries, size_t M, size_t B, double r);
	double enlargement_necessary(Candidate const& C, double abandon_after) const;
	double enlargement_necessary(Wedge const& other_wedge, double abandon_after) const;
	double get_new_ED(Candidate const& C, double abandon_after) const;
	double get_new_ED(Wedge const& other_wedge, double abandon_after) const;
	void enlarge(Candidate const& C);
	void enlarge(Wedge const& other_wedge);
};

#endif