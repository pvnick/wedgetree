#ifndef _WEDGE_H_
#define _WEDGE_H_

#include <memory>
#include <vector>
#include "candidate.h"

struct Wedge {
	//every node stores a wedge, which tightly bounds all of its descendent candidates
	std::vector<double> const& timeseries;
	size_t M;
	size_t B;
	double r;
	std::vector<double> U;
	std::vector<double> L;
	std::vector<double> UL_ED; //holds squared euclidean distance between individual U/L points
	double ED = 0; //squared euclidean distance between upper and lower bounds
	bool enlarged = false;
	static size_t id_counter; //holds the number of instantiated wedges
	size_t id; //set to the value of id_counter when this wedge was created
	Wedge(std::vector<double> const& timeseries, size_t M, size_t B, double r) :
		timeseries(timeseries),
		M(M),
		B(B),
		r(r),
		U(M, std::numeric_limits<double>::min()),
		L(M, std::numeric_limits<double>::max()),
		UL_ED(M, 0),
		id(++id_counter)
	{ }
	double enlargement_necessary(std::shared_ptr<Candidate> C, double abandon_after = std::numeric_limits<double>::max()) const {
		double enlargement = 0;
		if (enlarged) {
			for (size_t i = 0; i != M; ++i) {
				double val = C->series_normalized[i];
				if (val < L[i])
					enlargement += L[i] - val;
				else if (val > U[i])
					enlargement += val - U[i];
				if (enlargement > abandon_after)
					break;
			}
		}
		return enlargement;
	}
	double enlargement_necessary(Wedge const& other_wedge, double abandon_after = std::numeric_limits<double>::max()) const {
		double enlargement = 0;
		if (enlarged) {
			for (size_t i = 0; i != M; ++i) {
				if (other_wedge.L[i] < L[i])
					enlargement += L[i] - other_wedge.L[i];
				if (other_wedge.U[i] > U[i])
					enlargement += other_wedge.U[i] - U[i];
				if (enlargement > abandon_after)
					break;
			}
		}
		return enlargement;
	}
	double get_new_ED(std::shared_ptr<Candidate> C, double abandon_after = std::numeric_limits<double>::max()) const {
		//return the ED of the upper and lower bound if a candidate is added to this wedge
		double new_ED = ED;
		if (enlarged) {
			for (size_t i = 0; i != M; ++i) {
				bool outside_boundary = false;
				double new_U = U[i], new_L = L[i];
				if (C->series_normalized[i] > U[i]) {
					new_U = C->series_normalized[i];
					outside_boundary = true;
				}
				if (C->series_normalized[i] < L[i]) {
					new_L = C->series_normalized[i];
					outside_boundary = true;
				}
				if (outside_boundary) {
					new_ED -= UL_ED[i];
					new_ED += std::pow(new_U - new_L, 2);
				}
				if (new_ED > abandon_after)
					break;
			}
		}
		return new_ED;
	}
	double get_new_ED(Wedge const& other_wedge, double abandon_after = std::numeric_limits<double>::max()) const {
		//return the ED of the upper and lower bound of the union between this wedge and another wedge
		double new_ED = ED;
		if (enlarged) {
			for (size_t i = 0; i != M; ++i) {
				bool outside_boundary = false;
				double new_U = U[i], new_L = L[i];
				if (other_wedge.U[i] > U[i]) {
					new_U = other_wedge.U[i];
					outside_boundary = true;
				}
				if (other_wedge.L[i] < L[i]) {
					new_L = other_wedge.L[i];
					outside_boundary = true;
				}
				if (outside_boundary) {
					new_ED -= UL_ED[i];
					new_ED += std::pow(new_U - new_L, 2);
				}
				if (new_ED > abandon_after)
					break;
			}
		}
		return new_ED;
	}
	void enlarge(std::shared_ptr<Candidate> C) {
		for (size_t i = 0; i != M; ++i) {
			bool outside_boundary = false;
			if (C->series_normalized[i] > U[i]) {
				U[i] = C->series_normalized[i];
				outside_boundary = true;
			}
			if (C->series_normalized[i] < L[i]) {
				L[i] = C->series_normalized[i];
				outside_boundary = true;
			}
			if (outside_boundary) {
				ED -= UL_ED[i]; //subtract the old distance between the two bounds at this point
				UL_ED[i] = std::pow(U[i] - L[i], 2); //store the new distance
				ED += UL_ED[i]; //add the new distance between the two bounds at this point
			}
		}
		enlarged = true;
	}
	void enlarge(Wedge const& other_wedge) {
		for (size_t i = 0; i != M; ++i) {
			bool outside_boundary = false;
			if (other_wedge.U[i] > U[i]) {
				U[i] = other_wedge.U[i];
				outside_boundary = true;
			}
			if (other_wedge.L[i] < L[i]) {
				L[i] = other_wedge.L[i];
				outside_boundary = true;
			}
			if (outside_boundary) {
				ED -= UL_ED[i]; //subtract the old distance between the two bounds at this point
				UL_ED[i] = std::pow(U[i] - L[i], 2); //store the new distance
				ED += UL_ED[i]; //add the new distance between the two bounds at this point
			}
		}
		enlarged = true;
	}
};

#endif