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
#include <string>
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
	size_t R;
	Wedge W;
	enum NodeType { NODE_CANDIDATE, NODE_INTERNAL_WEDGE, NODE_LEAF_WEDGE };
	NodeType node_type;
public:
	bool is_tree_root = false;
	virtual bool insert_timeseries(Candidate&& C) = 0;
	virtual size_t get_height() const = 0;
	virtual void recalculate_wedge_lb_keough_envelopes() = 0;
	Node(std::vector<double> const& timeseries, size_t M, size_t B, double r, size_t R, NodeType node_type);
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
	CandidateNode(std::vector<double> const& timeseries, size_t M, size_t B, double r, size_t R);
	virtual bool insert_timeseries(Candidate&& C);
	virtual void recalculate_wedge_lb_keough_envelopes();
	bool can_contain_candidate(Candidate const& C) const;
	bool can_contain_wedge(Wedge const& wedge) const;
	std::list<Candidate> const& get_candidates() const;
	void merge_and_destroy_source(CandidateNode* src);
	void move_candidates(std::list<Candidate>* dest_list);
	size_t get_height() const;
};

class WedgeNode : public Node {
	//wedge tree internal nodes, holding up to B child nodes (wedge or candidate nodes)
	//splittable by its parent (always an internal wedge node)
protected:
	std::vector<std::unique_ptr<Node>> entries;
public:
	WedgeNode(std::vector<double> const& timeseries, size_t M, size_t B, double r, size_t R, NodeType node_type);
	virtual std::list<CandidateNode> get_merged_candidate_nodes(bool show_progress, unsigned int depth, size_t total_num_leaf_wedge_nodes, size_t* leaf_wedge_node_counter) const = 0;
	virtual size_t count_leaf_wedge_nodes() const = 0;
	virtual size_t count_candidate_nodes() const = 0;
	virtual void recalculate_wedge_lb_keough_envelopes();
	void clear_entries();
	std::vector<std::unique_ptr<Node>>& get_entries();
	void add_entry(std::unique_ptr<Node> entry);
	size_t get_min_enlargement_insertion_target(Candidate const& C) const;
	size_t get_height() const;
};

class LeafWedgeNode : public WedgeNode {
	//a "leaf" insofar as this type of node is at the edge of the tree and can only hold candidate nodes, not other wedge nodes
public:
	LeafWedgeNode(std::vector<double> const& timeseries, size_t M, size_t B, double r, size_t R);
	virtual bool insert_timeseries(Candidate&& C);
	std::list<CandidateNode> get_merged_candidate_nodes(bool show_progress, unsigned int depth, size_t total_num_leaf_wedge_nodes, size_t* leaf_wedge_node_counter) const;
	size_t count_leaf_wedge_nodes() const;
	size_t count_candidate_nodes() const;
};

class InternalWedgeNode : public WedgeNode {
	//internal wedge node which can only hold other wedge nodes
public:
	InternalWedgeNode(std::vector<double> const& timeseries, size_t M, size_t B, double r, size_t R);
	virtual bool insert_timeseries(Candidate&& C);
	void split_child(size_t target_entry_index);
	std::list<CandidateNode> get_merged_candidate_nodes(bool show_progress, unsigned int depth, size_t total_num_leaf_wedge_nodes, size_t* leaf_wedge_node_counter) const;
	size_t count_leaf_wedge_nodes() const;
	size_t count_candidate_nodes() const;
};

#endif