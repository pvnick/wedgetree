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

#include "wedge.h"

size_t Wedge::id_counter = 0;

Wedge::Wedge(std::vector<double> const& timeseries, size_t M, size_t B, double r) :
	timeseries(timeseries),
	M(M),
	B(B),
	r(r),
	id(++id_counter),
	boundary_positions(M),
	bound_comparator(boundary_positions)
{ 
	bound_pos_indices_sorted.reserve(M + 1); //reserve an extra slot for insertion before erase when resorting
	for (size_t i = 0; i != M; ++i) {
		bound_pos_indices_sorted.push_back(i);
	}
}

double Wedge::enlargement_necessary(Candidate const& C, double abandon_after = std::numeric_limits<double>::max()) const {
	double enlargement = 0;
	if (enlarged) {
		auto boundary_iter_end = bound_pos_indices_sorted.end();
		for (auto boundary_iter = bound_pos_indices_sorted.begin(); boundary_iter != boundary_iter_end; ++boundary_iter) {
			size_t c_index = *boundary_iter;
			BoundaryPosition const& bound = boundary_positions[c_index];
			double val = C.series_normalized[c_index];
			if (val < bound.Li)
				enlargement += bound.Li - val;
			else if (val > bound.Ui)
				enlargement += val - bound.Ui;
			if (enlargement > abandon_after)
				break;
		}
	}
	return enlargement;
}

double Wedge::enlargement_necessary(Wedge const& other_wedge, double abandon_after = std::numeric_limits<double>::max()) const {
	double enlargement = 0;
	if (enlarged) {
		auto boundary_iter_end = bound_pos_indices_sorted.end();
		for (auto boundary_iter = bound_pos_indices_sorted.begin(); boundary_iter != boundary_iter_end; ++boundary_iter) {
			size_t c_index = *boundary_iter;
			BoundaryPosition const& bound = boundary_positions[c_index];
			BoundaryPosition const& other_bound = other_wedge.boundary_positions[c_index];
			if (other_bound.Li < bound.Li)
				enlargement += bound.Li - other_bound.Li;
			if (other_bound.Ui > bound.Ui)
				enlargement += other_bound.Ui - bound.Ui;
			if (enlargement > abandon_after)
				break;
		}
	}
	return enlargement;
}

double Wedge::get_new_ED(Candidate const& C, double abandon_after = std::numeric_limits<double>::max()) const {
	//return the ED of the upper and lower bound if a candidate is added to this wedge
	double new_ED = ED;
	if (enlarged) {
		auto boundary_iter_end = bound_pos_indices_sorted.end();
		for (auto boundary_iter = bound_pos_indices_sorted.begin(); boundary_iter != boundary_iter_end; ++boundary_iter) {
			size_t c_index = *boundary_iter;
			BoundaryPosition const& bound = boundary_positions[c_index];
			double val = C.series_normalized[c_index];
			bool outside_boundary = false;
			double new_U = bound.Ui, new_L = bound.Li;
			if (C.series_normalized[c_index] > bound.Ui) {
				new_U = C.series_normalized[c_index];
				outside_boundary = true;
			}
			if (C.series_normalized[c_index] < bound.Li) {
				new_L = C.series_normalized[c_index];
				outside_boundary = true;
			}
			if (outside_boundary) {
				new_ED -= bound.ED;
				new_ED += std::pow(new_U - new_L, 2);
			}
			if (new_ED > abandon_after)
				break;
		}
	}
	return new_ED;
}

double Wedge::get_new_ED(Wedge const& other_wedge, double abandon_after = std::numeric_limits<double>::max()) const {
	//return the ED of the upper and lower bound if a candidate is added to this wedge
	double new_ED = ED;
	if (enlarged) {
		auto boundary_iter_end = bound_pos_indices_sorted.end();
		for (auto boundary_iter = bound_pos_indices_sorted.begin(); boundary_iter != boundary_iter_end; ++boundary_iter) {
			size_t c_index = *boundary_iter;
			BoundaryPosition const& bound = boundary_positions[c_index];
			BoundaryPosition const& other_bound = other_wedge.boundary_positions[c_index];
			bool outside_boundary = false;
			double new_U = bound.Ui, new_L = bound.Li;
			if (other_bound.Ui > bound.Ui) {
				new_U = other_bound.Ui;
				outside_boundary = true;
			}
			if (other_bound.Li < bound.Li) {
				new_L = other_bound.Li;
				outside_boundary = true;
			}
			if (outside_boundary) {
				new_ED -= bound.ED;
				new_ED += std::pow(new_U - new_L, 2);
			}
			if (new_ED > abandon_after)
				break;
		}
	}
	return new_ED;
}

void Wedge::enlarge(Candidate const& C) {
	for (size_t c_index = 0; c_index != M; ++c_index) {
		bool outside_boundary = false;
		BoundaryPosition& bound = boundary_positions[c_index];
		if (C.series_normalized[c_index] > bound.Ui) {
			bound.Ui = C.series_normalized[c_index];
			outside_boundary = true;
		}
		if (C.series_normalized[c_index] < bound.Li) {
			bound.Li = C.series_normalized[c_index];
			outside_boundary = true;
		}
		if (outside_boundary && bound.Ui - bound.Li > 0) {
			ED -= bound.ED; //subtract the old distance between the two bounds at this point
			bound.ED = std::pow(bound.Ui - bound.Li, 2); //store the new distance
			ED += bound.ED; //add the new distance between the two bounds at this point
			//boundary position may get moved back in the sorted vector
			/*auto reposition_iter = bound_pos_indices_sorted.begin() + c_index;
			auto insertion_position = std::upper_bound(reposition_iter, bound_pos_indices_sorted.end(), bound.ED, bound_comparator);
			if (reposition_iter != insertion_position) {
				bound_pos_indices_sorted.insert(insertion_position, c_index);
				bound_pos_indices_sorted.erase(bound_pos_indices_sorted.begin() + c_index);
			}*/
			//need to sort at this point
			//sorted_boundary.increase(heap_handles[i]);
		}
	}
	std::sort(bound_pos_indices_sorted.begin(), bound_pos_indices_sorted.end(), bound_comparator);
	/*for (size_t i : bound_pos_indices_sorted)
		std::cout << boundary_positions[i].ED << " ";
	std::cout << std::endl;*/
	enlarged = true;
}

void Wedge::enlarge(Wedge const& other_wedge) {
	for (size_t c_index = 0; c_index != M; ++c_index) {
		bool outside_boundary = false;
		BoundaryPosition& bound = boundary_positions[c_index];
		BoundaryPosition const& other_bound = other_wedge.boundary_positions[c_index];
		if (other_bound.Ui > bound.Ui) {
			bound.Ui = other_bound.Ui;
			outside_boundary = true;
		}
		if (other_bound.Li < bound.Li) {
			bound.Li = other_bound.Li;
			outside_boundary = true;
		}
		if (outside_boundary) {
			ED -= bound.ED; //subtract the old distance between the two bounds at this point
			bound.ED = std::pow(bound.Ui - bound.Li, 2); //store the new distance
			ED += bound.ED; //add the new distance between the two bounds at this point
			//need to sort here
			//sorted_boundary.increase(heap_handles[i]);
		}
	}
	std::sort(bound_pos_indices_sorted.begin(), bound_pos_indices_sorted.end(), bound_comparator);
	/*for (size_t i : bound_pos_indices_sorted)
		std::cout << boundary_positions[i].ED << " ";
	std::cout << std::endl;*/
	enlarged = true;
}