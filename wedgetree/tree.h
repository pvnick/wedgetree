#ifndef HEADER_4EC913BBBD921EB5
#define HEADER_4EC913BBBD921EB5

#include <iostream>
#include <vector>
#include <stdexcept>
#include <limits>
#include <iomanip>
#include <cmath>
#include <algorithm>

struct Wedge {
	std::vector<double> const& timeseries;
	size_t M;
	size_t B;
	double r;
	std::vector<double> U;
	std::vector<double> L;
	double ED; //squared euclidean distance between upper and lower bounds
	bool enlarged = false;
	Wedge(std::vector<double> const& timeseries, size_t M, size_t B, double r) :
		timeseries(timeseries),
		M(M),
		B(B),
		r(r),
		U(M, std::numeric_limits<double>::min()),
		L(M, std::numeric_limits<double>::max())
	{ }
	double enlargement_necessary(size_t C_index) const {
		double enlargement = 0;
		if (enlarged) {
			for (size_t i = 0; i != M; ++i) {
				double val = timeseries[C_index + i];
				if (val < L[i])
					enlargement += L[i] - val;
				else if (val > U[i])
					enlargement += val - U[i];
			}
		}
		return enlargement;
	}
	double enlargement_necessary(Wedge const& other_wedge) const {
		double enlargement = 0;
		if (enlarged) {
			for (size_t i = 0; i != M; ++i) {
				if (other_wedge.L[i] < L[i])
					enlargement += L[i] - other_wedge.L[i];
				if (other_wedge.U[i] > U[i])
					enlargement += other_wedge.U[i] - U[i];
			}
		}
		return enlargement;
	}
	double get_new_ED(size_t C_index) const {
		double new_ED = 0;
		if (enlarged) {
			for (size_t i = 0; i != M; ++i) {
				double val = timeseries[C_index + i];
				double U_val = std::max(U[i], val);
				double L_val = std::min(L[i], val);
				new_ED += std::pow(U_val, 2) - std::pow(L_val, 2);
			}
		}
		return new_ED;
	}
	double get_new_ED(Wedge const& other_wedge) const {
		double new_ED = 0;
		if (enlarged) {
			for (size_t i = 0; i != M; ++i) {
				double U_val = std::max(U[i], other_wedge.U[i]);
				double L_val = std::min(L[i], other_wedge.L[i]);
				new_ED += std::pow(U_val, 2) - std::pow(L_val, 2);
			}
		}
		return new_ED;
	}
	void enlarge(size_t C_index) {
		ED = 0;
		for (size_t i = 0; i != M; ++i) {
			double val = timeseries[C_index + i];
			U[i] = std::max(U[i], val);
			L[i] = std::min(L[i], val);
			ED += std::pow(U[i], 2) - std::pow(L[i], 2);
		}
		enlarged = true;
	}
	void enlarge(Wedge const& other_wedge) {
		ED = 0;
		for (size_t i = 0; i != M; ++i) {
			U[i] = std::max(U[i], other_wedge.U[i]);
			L[i] = std::min(L[i], other_wedge.L[i]);
			ED += std::pow(U[i], 2) - std::pow(L[i], 2);
		}
		enlarged = true;
	}
};

class Node {
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
	virtual bool insert_timeseries(size_t C_index) = 0;
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
	//exists as a leaf in the wedge tree. can hold unlimited number of candidates, but wedge ED must be under r
	//not splittable
private:
	std::vector<size_t> C_set; //holds the starting positions for each of the candidates
public:
	CandidateNode(std::vector<double> const& timeseries, size_t M, size_t B, double r) :
		Node(timeseries, M, B, r, NODE_CANDIDATE)
	{ }
	bool can_contain_candidate(size_t C_index) {
		return W.get_new_ED(C_index) <= r;
	}
	bool insert_timeseries(size_t C_index) {
		C_set.push_back(C_index);
		W.enlarge(C_index);
		return false;
	}
	size_t get_height() {
		return 1;
	}
};

class WedgeNode : public Node {
	//exists as a leaf in the wedge tree, holding up to B candidate nodes
	//splittable
protected:
	std::vector<Node*> entries;
public:
	WedgeNode(std::vector<double> const& timeseries, size_t M, size_t B, double r, NodeType node_type) :
		Node(timeseries, M, B, r, node_type)
	{ }
	void clear_entries() {
		entries.clear();
	}
	std::vector<Node*> const& get_entries() {
		return entries;
	}
	void add_entry(Node* entry) {
		//note: this function does not check for entry count overflows
		entries.push_back(entry);
		W.enlarge(entry->get_wedge());
	}
	size_t get_height() {
		size_t max_height = 0;
		for (Node* n : entries) {
			size_t height = n->get_height();
			if (height > max_height)
				max_height = height;
		}
		return 1 + max_height;
	}
};

class LeafWedgeNode : public WedgeNode {
public:
	LeafWedgeNode(std::vector<double> const& timeseries, size_t M, size_t B, double r) :
		WedgeNode(timeseries, M, B, r, NODE_LEAF_WEDGE)
	{ }
	bool insert_timeseries(size_t C_index) {
		//returns true if current node's entry count exceeds B and must be split, false otherwise
		bool split_required = false;
		if (entries.empty()) {
			CandidateNode* n = new CandidateNode(timeseries, M, B, r);
			n->insert_timeseries(C_index);
			entries.push_back(n);
		}
		else {
			size_t target_index;
			double min_enlargement = std::numeric_limits<double>::max();
			for (size_t i = 0; i != entries.size(); ++i) {
				//insert to the entry which requires the minimum enlargement to contain the new candidate
				double enlargement = entries[i]->get_wedge().enlargement_necessary(C_index);
				if (enlargement < min_enlargement) {
					min_enlargement = enlargement;
					target_index = i;
				}
			}
			CandidateNode* target = static_cast<CandidateNode*>(entries[target_index]);
			if (target->can_contain_candidate(C_index)) {
				target->insert_timeseries(C_index);
			} else {
				//insertion to entry at a leaf node can be done only when the enlargement satisfies the condition W_UL <= r; otherwise
				//we have to create a new entry at that leaf node
				CandidateNode* n = new CandidateNode(timeseries, M, B, r);
				n->insert_timeseries(C_index);
				entries.push_back(n);
			}
		}
		if (entries.size() > B) {
			//not enough space to contain the new entry in the current node, need to split
			split_required = true;
		}
		W.enlarge(C_index);
		return split_required;
	}
};

class InternalWedgeNode : public WedgeNode {
public:
	InternalWedgeNode(std::vector<double> const& timeseries, size_t M, size_t B, double r) :
		WedgeNode(timeseries, M, B, r, NODE_INTERNAL_WEDGE)
	{ }
	void split_child(size_t target_entry_index) {
		WedgeNode* target = static_cast<WedgeNode*>(entries[target_entry_index]);
		size_t E1_index, E2_index;
		double max_union_ED = 0;
		//choose the two entries whose combined wedges have the largest euclidean distance
		std::vector<Node*> const& target_entries = target->get_entries();
		for (size_t i = 0; i != target_entries.size(); ++i) {
			WedgeNode* E1_tmp = static_cast<WedgeNode*>(target_entries[i]);
			for (size_t j = i + 1; j != target_entries.size(); ++j) {
				WedgeNode* E2_tmp = static_cast<WedgeNode*>(target_entries[j]);
				double ED_tmp = E1_tmp->get_wedge().get_new_ED(E2_tmp->get_wedge());
				if (ED_tmp > max_union_ED) {
					E1_index = i;
					E2_index = j;
					max_union_ED = ED_tmp;
				}
			}
		}
		//create two new entries to contain E1 and E2
		WedgeNode* new_entry1;
		WedgeNode* new_entry2;
		switch (target->get_node_type()) {
		case NODE_INTERNAL_WEDGE:
			new_entry1 = new InternalWedgeNode(timeseries, M, B, r);
			new_entry2 = new InternalWedgeNode(timeseries, M, B, r);
			break;
		case NODE_LEAF_WEDGE:
			new_entry1 = new LeafWedgeNode(timeseries, M, B, r);
			new_entry2 = new LeafWedgeNode(timeseries, M, B, r);
			break;
		default:
			throw new std::logic_error("Unexpected node type");
		}
		new_entry1->add_entry(target_entries[E1_index]);
		new_entry2->add_entry(target_entries[E2_index]);
		//add the rest of target's entries to either of the two new entries to minimize enlargement
		for (size_t i = 0; i != target_entries.size(); ++i) {
			if (i != E1_index && i != E2_index) {
				Node* entry_to_add = target_entries[i];
				if (new_entry1->get_wedge().enlargement_necessary(entry_to_add->get_wedge()) > new_entry2->get_wedge().enlargement_necessary(entry_to_add->get_wedge())) {
					new_entry2->add_entry(entry_to_add);
				}
				else {
					new_entry1->add_entry(entry_to_add);
				}
			}
		}
		//add the new entries, replacing the old one (increases entry count by one)
		//this operation does not enlarge the current wedge
		target->clear_entries();
		entries[target_entry_index] = new_entry1;
		entries.push_back(new_entry2);
		//the old target is no longer needed
		delete target;
	}
	bool insert_timeseries(size_t C_index) {
		//returns true if current node's entry count exceeds B and must be split, false otherwise
		bool split_required = false;
		if (entries.empty()) {
			LeafWedgeNode* n = new LeafWedgeNode(timeseries, M, B, r);
			//no need to check for if a split is necessary since the candidate is guaranteed to insertable into a fresh node
			n->insert_timeseries(C_index);
			add_entry(n);
		}
		else {
			size_t target_index;
			double min_enlargement = std::numeric_limits<double>::max();
			for (size_t i = 0; i != entries.size(); ++i) {
				//insert to the entry which requires the minimum enlargement to contain the new candidate
				double enlargement = entries[i]->get_wedge().enlargement_necessary(C_index);
				if (enlargement < min_enlargement) {
					min_enlargement = enlargement;
					target_index = i;
				}
			}
			WedgeNode* target = static_cast<WedgeNode*>(entries[target_index]);
			if (target->insert_timeseries(C_index)) {
				//target experienced overflow and must be split
				split_child(target_index);
			}
			if (entries.size() > B) {
				//not enough space to contain the new entry in the current node, need to split
				split_required = true;
			}
		}
		W.enlarge(C_index);
		return split_required;
	}
};

class WedgeTree {
private:
	std::vector<double> const& timeseries;
	size_t M;
	size_t B;
	double r;
	InternalWedgeNode* root;
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
			insert_timeseries(C_index);
			//std::vector<Node*> entries = root->get_entries();
			//for (size_t i = 0; i != entries.size(); ++i)
			//    std::cout << "i:" << i << " (" << entries[i]->get_height() << ") | ";
			//std::cout << std::endl;
			if (verbose && (C_index % 1000) == 0) {
				std::cout << "\r" << std::setprecision(4) << (100.0 * C_index / ts_size) << "% completed       ";
				std::flush(std::cout);
			}
		}
		if (verbose)
			std::cout << "\r100% completed" << std::endl;
	}
	void insert_timeseries(size_t C_index) {
		if (root->insert_timeseries(C_index)) {
			//root overflow. add a new root to contain the old root, then split the old root
			InternalWedgeNode* old_root = root;
			root = new InternalWedgeNode(timeseries, M, B, r);
			root->add_entry(old_root);
			root->split_child(0);
		}
	}
};

#endif // header guard
