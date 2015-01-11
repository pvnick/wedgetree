#ifndef HEADER_211D3D8E8FC34194A31569F870F581A0
#define HEADER_211D3D8E8FC34194A31569F870F581A0

#include <iostream>
#include <vector>
#include <stdexcept>
#include <limits>
#include <iomanip>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <memory>
#include <list>
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
	bool is_tree_root = false;
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
	std::list<Candidate> C_set; //holds the starting positions for each of the candidates
public:
	CandidateNode(std::vector<double> const& timeseries, size_t M, size_t B, double r);
	bool can_contain_candidate(Candidate const& C) const;
	bool can_contain_wedge(Wedge const& wedge) const;
	bool insert_timeseries(Candidate&& C);
	std::list<Candidate> const& get_candidates() const;
	void merge_and_destroy_source(CandidateNode* src);
	void move_candidates(std::list<Candidate>* dest_list);
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
	virtual std::list<CandidateNode> get_merged_candidate_nodes(bool show_progress, size_t total_num_leaf_wedge_nodes, size_t* leaf_wedge_node_counter) const = 0;
	virtual size_t count_leaf_wedge_nodes() = 0;
	virtual size_t count_candidate_nodes() = 0;
};

class LeafWedgeNode : public WedgeNode {
	//a "leaf" insofar as this type of node is at the edge of the tree and can only hold candidate nodes, not other wedge nodes
    //not splittable
public:
	LeafWedgeNode(std::vector<double> const& timeseries, size_t M, size_t B, double r);
	bool insert_timeseries(Candidate&& C);
	std::list<CandidateNode> get_merged_candidate_nodes(bool show_progress, size_t total_num_leaf_wedge_nodes, size_t* leaf_wedge_node_counter) const;
	size_t count_leaf_wedge_nodes();
	size_t count_candidate_nodes();
};

class InternalWedgeNode : public WedgeNode {
	//internal wedge node which can only hold other wedge nodes
	//splittable
public:
	InternalWedgeNode(std::vector<double> const& timeseries, size_t M, size_t B, double r);
	void split_child(size_t target_entry_index);
	bool insert_timeseries(Candidate&& C);
	std::list<CandidateNode> get_merged_candidate_nodes(bool show_progress, size_t total_num_leaf_wedge_nodes, size_t* leaf_wedge_node_counter) const;
	size_t count_leaf_wedge_nodes();
	size_t count_candidate_nodes();
};

#endif