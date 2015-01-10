#ifndef HEADER_211D3D8E8FC34194A31569F870F581A0
#define HEADER_211D3D8E8FC34194A31569F870F581A0

#include <iostream>
#include <vector>
#include <stdexcept>
#include <limits>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <memory>
#include "wedge.h"

class Node {
	//abstract base Node class
protected:
	std::vector<double> const& timeseries;
	static size_t id_counter;
	size_t id;
	size_t M;
	size_t B;
	double r;
	Wedge W;
	enum NodeType { NODE_CANDIDATE, NODE_INTERNAL_WEDGE, NODE_LEAF_WEDGE };
	NodeType node_type;
public:
	virtual bool insert_timeseries(Candidate&& C) = 0;
	virtual size_t get_height() const = 0;
	Node(std::vector<double> const& timeseries, size_t M, size_t B, double r, NodeType node_type);
	Wedge const& get_wedge() const;
	NodeType get_node_type() const;
};

class CandidateNode : public Node {
	//exists as a true leaf in the wedge tree, held by LeafWedgeNodes. 
	//can hold unlimited number of candidates, but wedge ED must be under r.
	//not splittable
private:
	std::vector<Candidate> C_set; //holds the starting positions for each of the candidates
public:
	CandidateNode(std::vector<double> const& timeseries, size_t M, size_t B, double r);
	bool can_contain_candidate(Candidate const& C) const;
	bool insert_timeseries(Candidate&& C);
	size_t get_height() const;
};

class WedgeNode : public Node {
	//wedge tree internal nodes, holding up to B child nodes (wedge or candidate nodes)
protected:
	std::vector<std::shared_ptr<Node>> entries;
public:
	WedgeNode(std::vector<double> const& timeseries, size_t M, size_t B, double r, NodeType node_type);
	void clear_entries();
	std::vector<std::shared_ptr<Node>> const& get_entries() const;
	void add_entry(std::shared_ptr<Node> entry);
	size_t get_min_enlargement_insertion_target(Candidate const& C) const;
	size_t get_height() const;
	//std::shared_ptr<std::vector<CandidateNode>> getMergedCandidateNodes() = 0;
};

class LeafWedgeNode : public WedgeNode {
	//a "leaf" insofar as this type of node is at the edge of the tree and can only hold candidate nodes, not other wedge nodes
    //not splittable
public:
	LeafWedgeNode(std::vector<double> const& timeseries, size_t M, size_t B, double r);
	bool insert_timeseries(Candidate&& C);
	//std::shared_ptr<std::vector<CandidateNode>> getMergedCandidateNodes();
};

class InternalWedgeNode : public WedgeNode {
	//internal wedge node which can only hold other wedge nodes
	//splittable
public:
	InternalWedgeNode(std::vector<double> const& timeseries, size_t M, size_t B, double r);
	void split_child(size_t target_entry_index);
	bool insert_timeseries(Candidate&& C);
	//std::shared_ptr<std::vector<CandidateNode>> getMergedCandidateNodes();
};

#endif