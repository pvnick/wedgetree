#include "wedge.h"

size_t Wedge::id_counter = 0;

Wedge::Wedge(std::vector<double> const& timeseries, size_t M, size_t B, double r) :
	timeseries(timeseries),
	M(M),
	B(B),
	r(r),
	id(++id_counter)
{ 
	boundary_handles_by_c_index.reserve(M);
	for (size_t i = 0; i != M; ++i) {
		std::shared_ptr<BoundaryPosition> bound = std::make_shared<BoundaryPosition>(i);
		boundary_handles_by_c_index.push_back(bound);
		sorted_boundary.push_back(bound);
	}
}

double Wedge::enlargement_necessary(Candidate const& C, double abandon_after = std::numeric_limits<double>::max()) const {
	double enlargement = 0;
	if (enlarged) {
		auto boundary_iter_end = sorted_boundary.end();
		for (auto boundary_iter = sorted_boundary.begin(); boundary_iter != boundary_iter_end; ++boundary_iter) {
			BoundaryPosition const& bound = **boundary_iter;
			double val = C.series_normalized[bound.c_index];
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
		auto boundary_iter_end = sorted_boundary.end();
		for (auto boundary_iter = sorted_boundary.begin(); boundary_iter != boundary_iter_end; ++boundary_iter) {
			BoundaryPosition const& bound = **boundary_iter;
			BoundaryPosition const& other_bound = *(other_wedge.boundary_handles_by_c_index[bound.c_index]);
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
		sorter::const_iterator boundary_iter_end = sorted_boundary.end();
		for (auto boundary_iter = sorted_boundary.begin(); boundary_iter != boundary_iter_end; ++boundary_iter) {
			BoundaryPosition const& bound = **boundary_iter;
			double val = C.series_normalized[bound.c_index];
			bool outside_boundary = false;
			double new_U = bound.Ui, new_L = bound.Li;
			if (C.series_normalized[bound.c_index] > bound.Ui) {
				new_U = C.series_normalized[bound.c_index];
				outside_boundary = true;
			}
			if (C.series_normalized[bound.c_index] < bound.Li) {
				new_L = C.series_normalized[bound.c_index];
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
		sorter::const_iterator boundary_iter_end = sorted_boundary.end();
		for (auto boundary_iter = sorted_boundary.begin(); boundary_iter != boundary_iter_end; ++boundary_iter) {
			BoundaryPosition const& bound = **boundary_iter;
			BoundaryPosition const& other_bound = *(other_wedge.boundary_handles_by_c_index[bound.c_index]);
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
	for (size_t i = 0; i != M; ++i) {
		bool outside_boundary = false;
		BoundaryPosition& bound = *(boundary_handles_by_c_index[i]);
		if (C.series_normalized[i] > bound.Ui) {
			bound.Ui = C.series_normalized[i];
			outside_boundary = true;
		}
		if (C.series_normalized[i] < bound.Li) {
			bound.Li = C.series_normalized[i];
			outside_boundary = true;
		}
		if (outside_boundary) {
			ED -= bound.ED; //subtract the old distance between the two bounds at this point
			bound.ED = std::pow(bound.Ui - bound.Li, 2); //store the new distance
			ED += bound.ED; //add the new distance between the two bounds at this point
			//need to sort at this point
			//sorted_boundary.increase(heap_handles[i]);
		}
	}
	enlarged = true;
}

void Wedge::enlarge(Wedge const& other_wedge) {
	for (size_t i = 0; i != M; ++i) {
		bool outside_boundary = false;
		BoundaryPosition& bound = *boundary_handles_by_c_index[i];
		BoundaryPosition& other_bound = *(other_wedge.boundary_handles_by_c_index[i]);
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
	enlarged = true;
}