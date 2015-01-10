#ifndef HEADER_A3E13BE6A09A41989CB88235F7155A45
#define HEADER_A3E13BE6A09A41989CB88235F7155A45

#include <memory>
#include <vector>
#include "candidate.h"

class Wedge {
	//every node stores a wedge, which tightly bounds all of its descendent candidates
private:
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