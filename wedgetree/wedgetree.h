#ifndef _WEDGETREE_H_
#define _WEDGETREE_H_

#include <iostream>
#include <vector>
#include <stdexcept>
#include <limits>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <memory>
#include "candidate.h"
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
	virtual bool insert_timeseries(std::shared_ptr<Candidate> C) = 0;
	virtual size_t get_height() = 0;
	Node(std::vector<double> const& timeseries, size_t M, size_t B, double r, NodeType node_type) :
		timeseries(timeseries),
		id(++id_counter),
		M(M),
		B(B),
		r(r),
		W(timeseries, M, B, r),
		node_type(node_type)
	{ }
	Wedge const& get_wedge() {
		return W;
	}
	NodeType get_node_type() {
		return node_type;
	}
};

class CandidateNode : public Node {
	//exists as a true leaf in the wedge tree, held by LeafWedgeNodes. 
	//can hold unlimited number of candidates, but wedge ED must be under r.
	//not splittable
private:
	std::vector<std::shared_ptr<Candidate>> C_set; //holds the starting positions for each of the candidates
public:
	CandidateNode(std::vector<double> const& timeseries, size_t M, size_t B, double r) :
		Node(timeseries, M, B, r, NODE_CANDIDATE)
	{ }
	bool can_contain_candidate(std::shared_ptr<Candidate> C) {
		double ED_after_insertion = W.get_new_ED(C, r);
		return ED_after_insertion <= r;
	}
	bool insert_timeseries(std::shared_ptr<Candidate> C) {
		C_set.push_back(C);
		W.enlarge(C);
		return false;
	}
	size_t get_height() {
		return 1;
	}
};

class WedgeNode : public Node {
	//wedge tree internal nodes, holding up to B child nodes (wedge or candidate nodes)
	//splittable
protected:
	std::vector<std::shared_ptr<Node>> entries;
public:
	WedgeNode(std::vector<double> const& timeseries, size_t M, size_t B, double r, NodeType node_type) :
		Node(timeseries, M, B, r, node_type)
	{ }
	void clear_entries() {
		entries.clear();
	}
	std::vector<std::shared_ptr<Node>> const& get_entries() {
		return entries;
	}
	void add_entry(std::shared_ptr<Node> entry) {
		//note: this function does not check for entry count overflows
		entries.push_back(entry);
		W.enlarge(entry->get_wedge());
	}
	size_t get_min_enlargement_insertion_target(std::shared_ptr<Candidate> C) {
		//this call assumes at least one entry present and will demonstrate undefined behavior otherwise
		size_t target_index = 0;
		double min_enlargement = entries[0]->get_wedge().enlargement_necessary(C, std::numeric_limits<double>::max());
		for (size_t i = 1; i != entries.size(); ++i) {
			//insert to the entry which requires the minimum enlargement to contain the new candidate
			double enlargement = entries[i]->get_wedge().enlargement_necessary(C, min_enlargement);
			if (enlargement < min_enlargement) {
				min_enlargement = enlargement;
				target_index = i;
			}
		}
		return target_index;
	}
	size_t get_height() {
		size_t max_height = 0;
		for (std::shared_ptr<Node> n : entries) {
			size_t height = n->get_height();
			if (height > max_height)
				max_height = height;
		}
		return 1 + max_height;
	}
};

class LeafWedgeNode : public WedgeNode {
	//a "leaf" insofar as this type of node can only hold candidate nodes, not other wedge nodes
public:
	LeafWedgeNode(std::vector<double> const& timeseries, size_t M, size_t B, double r) :
		WedgeNode(timeseries, M, B, r, NODE_LEAF_WEDGE)
	{ }
	bool insert_timeseries(std::shared_ptr<Candidate> C) {
		//returns true if current node's entry count exceeds B and must be split, false otherwise
		bool split_required = false;
		if (entries.empty()) {
			std::shared_ptr<CandidateNode>n = std::make_shared<CandidateNode>(timeseries, M, B, r);
			n->insert_timeseries(C);
			entries.push_back(n);
		}
		else {
			size_t target_index = get_min_enlargement_insertion_target(C);
			std::shared_ptr<CandidateNode> target = std::static_pointer_cast<CandidateNode>(entries[target_index]);
			if (target->can_contain_candidate(C)) {
				target->insert_timeseries(C);
			} else {
				//insertion to entry at a leaf node can be done only when the enlargement satisfies the condition W_UL <= r; otherwise
				//we have to create a new entry at that leaf node
				std::shared_ptr<CandidateNode> n = std::make_shared<CandidateNode>(timeseries, M, B, r);
				n->insert_timeseries(C);
				entries.push_back(n);
			}
		}
		if (entries.size() > B) {
			//not enough space to contain the new entry in the current node, need to split
			split_required = true;
		}
		W.enlarge(C);
		return split_required;
	}
};

class InternalWedgeNode : public WedgeNode {
	//internal wedge node which can only hold other wedge nodes
public:
	InternalWedgeNode(std::vector<double> const& timeseries, size_t M, size_t B, double r) :
		WedgeNode(timeseries, M, B, r, NODE_INTERNAL_WEDGE)
	{ }
	void split_child(size_t target_entry_index) {
		std::shared_ptr<WedgeNode> target = std::static_pointer_cast<WedgeNode>(entries[target_entry_index]);
		size_t E1_index, E2_index;
		double max_union_ED = 0;
		//choose the two entries whose combined wedges have the largest euclidean distance
		std::vector<std::shared_ptr<Node>> const& target_entries = target->get_entries();
		//create a flattened 2-dimensional array which holds euclidean distances between entry wedges
		std::vector<double> union_EDs(target_entries.size() * target_entries.size(), 0); 
		size_t num_entries = target_entries.size();
		for (size_t i = 0; i != num_entries; ++i) {
			std::shared_ptr<WedgeNode> E1_tmp = std::static_pointer_cast<WedgeNode>(target_entries[i]);
			for (size_t j = i + 1; j != num_entries; ++j) {
				std::shared_ptr<WedgeNode> E2_tmp = std::static_pointer_cast<WedgeNode>(target_entries[j]);
				double ED_tmp = E1_tmp->get_wedge().get_new_ED(E2_tmp->get_wedge(), max_union_ED);
				size_t index1 = i * num_entries + j;
				size_t index2 = j * num_entries + i;
				union_EDs[index1] = union_EDs[index2] = ED_tmp;
				if (ED_tmp >= max_union_ED) { //GTE here guarantees that E1_index and E2_index will be set at least once before these loops complete
					E1_index = i;
					E2_index = j;
					max_union_ED = ED_tmp;
				}
			}
		}
		//create two new entries to contain E1 and E2
		std::shared_ptr<WedgeNode> new_entry1;
		std::shared_ptr<WedgeNode> new_entry2;
		switch (target->get_node_type()) {
		case NODE_INTERNAL_WEDGE:
			new_entry1 = std::make_shared<InternalWedgeNode>(timeseries, M, B, r);
			new_entry2 = std::make_shared<InternalWedgeNode>(timeseries, M, B, r);
			break;
		case NODE_LEAF_WEDGE:
			new_entry1 = std::make_shared<LeafWedgeNode>(timeseries, M, B, r);
			new_entry2 = std::make_shared<LeafWedgeNode>(timeseries, M, B, r);
			break;
		default:
			throw new std::logic_error("Unexpected node type");
		}
		new_entry1->add_entry(target_entries[E1_index]);
		new_entry2->add_entry(target_entries[E2_index]);
		//add the rest of target's entries to either of the two new entries, based on their enlargement of the entry's original wedge
		for (size_t i = 0; i != target_entries.size(); ++i) {
			if (i != E1_index && i != E2_index) {
				std::shared_ptr<Node> entry_to_add = target_entries[i];
				double new_entry1_ED = union_EDs[i * num_entries + E1_index];
				double new_entry2_ED = union_EDs[i * num_entries + E2_index];
				if (new_entry1_ED < new_entry2_ED) {
					new_entry1->add_entry(entry_to_add);
				}
				else {
					new_entry2->add_entry(entry_to_add);
				}
			}
		}
		//add the new entries, replacing the old one (increases entry count by one)
		//this operation does not enlarge the current wedge
		entries[target_entry_index] = new_entry1;
		entries.push_back(new_entry2);
	}
	bool insert_timeseries(std::shared_ptr<Candidate> C) {
		//returns true if current node's entry count exceeds B and must be split, false otherwise
		bool split_required = false;
		if (entries.empty()) {
			std::shared_ptr<LeafWedgeNode> n = std::make_shared<LeafWedgeNode>(timeseries, M, B, r);
			//no need to check for if a split is necessary since the candidate is guaranteed to insertable into a fresh node
			n->insert_timeseries(C);
			add_entry(n);
		}
		else {
			size_t target_index = get_min_enlargement_insertion_target(C);
			std::shared_ptr<WedgeNode> target = std::static_pointer_cast<WedgeNode>(entries[target_index]);
			if (target->insert_timeseries(C)) {
				//target experienced overflow and must be split
				split_child(target_index);
			}
			if (entries.size() > B) {
				//not enough space to contain the new entry in the current node, need to split
				split_required = true;
			}
		}
		W.enlarge(C);
		return split_required;
	}
};

class WedgeTree {
private:
	std::vector<double> const& timeseries;
	size_t M;
	size_t B;
	double r;
	std::shared_ptr<InternalWedgeNode> root;
	bool verbose;
public:
	WedgeTree(std::vector<double> const& timeseries, size_t M, size_t B, double r, bool verbose = false) :
		timeseries(timeseries),
		M(M),
		B(B),
		r(r),
		root(new InternalWedgeNode(timeseries, M, B, r)),
		verbose(verbose)
	{
		if (timeseries.size() < M)
			throw new std::invalid_argument("Motif length M must be smaller than the length of the input timeseries");
		if (verbose) {
			std::cout << "Loading timeseries into wedge tree" << std::endl;
		}
		size_t ts_size = timeseries.size();
		for (size_t C_index = 0; C_index + M - 1 != ts_size; ++C_index) {
			std::shared_ptr<Candidate> C = std::make_shared<Candidate>(timeseries, C_index, M);
			insert_timeseries(C);
			/*std::vector<std::shared_ptr<Node>> entries = root->get_entries();
			if ((C_index % 100) == 0) {
				std::cout << C_index << " - ";
				for (size_t i = 0; i != entries.size(); ++i)
					std::cout << "i:" << i << " (" << entries[i]->get_height() << ") | ";
				std::cout << std::endl;
			}*/
			if (verbose && (C_index % 1000) == 0) {
				printf("\r%.2f%% completed", (100.0 * C_index / ts_size));
				std::flush(std::cout);
			}
		}
		if (verbose)
			std::cout << "\r100.0% completed" << std::endl;
	}
	void insert_timeseries(std::shared_ptr<Candidate> C) {
		if (root->insert_timeseries(C)) {
			//root overflow. add a new root to contain the old root, then split the old root
			std::shared_ptr<InternalWedgeNode> old_root = root;
			root = std::make_shared<InternalWedgeNode>(timeseries, M, B, r);
			root->add_entry(old_root);
			root->split_child(0);
		}
	}
};

#endif
