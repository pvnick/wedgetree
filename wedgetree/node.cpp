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

#include "node.h"

size_t Node::id_counter = 0;

Node::Node(std::vector<double> const& timeseries, size_t M, size_t B, double r, size_t R, NodeType node_type) :
	timeseries(timeseries),
	id(++id_counter),
	M(M),
	B(B),
	r(r),
	R(R),
	W(timeseries, M, B, r, R),
	node_type(node_type)
{ }

Wedge const& Node::get_wedge() const {
	return W;
}

Node::NodeType Node::get_node_type() const {
	return node_type;
}

CandidateNode::CandidateNode(std::vector<double> const& timeseries, size_t M, size_t B, double r, size_t R) :
	Node(timeseries, M, B, r, R, NODE_CANDIDATE)
{ }

bool CandidateNode::can_contain_candidate(Candidate const& C) const {
	double ED_after_insertion = W.get_new_ED(C, r);
	return ED_after_insertion <= r;
}

bool CandidateNode::can_contain_wedge(Wedge const& wedge) const {
	double ED_after_insertion = W.get_new_ED(wedge, r);
	return ED_after_insertion <= r;
}

bool CandidateNode::insert_timeseries(Candidate&& C) {
	W.enlarge(C);
	C_set.push_back(std::forward<Candidate>(C));
	return false;
}

void CandidateNode::move_candidates(std::list<Candidate>* dest_list) {
	std::move(C_set.begin(), C_set.end(), std::back_inserter(*dest_list));
	C_set.clear();
}

void CandidateNode::merge_and_destroy_source(CandidateNode* src) {
	W.enlarge(src->get_wedge());
	src->move_candidates(&C_set);
}

std::list<Candidate> const& CandidateNode::get_candidates() const {
	return C_set;
}

size_t CandidateNode::get_height() const {
	return 1;
}

void CandidateNode::recalculate_wedge_lb_keough_envelopes() {
	W.recalculate_DTW_envelope();
}

WedgeNode::WedgeNode(std::vector<double> const& timeseries, size_t M, size_t B, double r, size_t R, NodeType node_type) :
	Node(timeseries, M, B, r, R, node_type)
{ 
	entries.reserve(B);
}

void WedgeNode::clear_entries() {
	entries.clear();
}

std::vector<std::unique_ptr<Node>>& WedgeNode::get_entries() {
	return entries;
}

void WedgeNode::add_entry(std::unique_ptr<Node> entry) {
	//note: this function does not check for entry count overflows
	W.enlarge(entry->get_wedge());
	entries.push_back(std::move(entry));
}

size_t WedgeNode::get_min_enlargement_insertion_target(Candidate const& C) const {
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

size_t WedgeNode::get_height() const {
	size_t max_height = 0;
	for (std::unique_ptr<Node> const& n : entries) {
		size_t height = n->get_height();
		if (height > max_height)
			max_height = height;
	}
	return 1 + max_height;
}

void WedgeNode::recalculate_wedge_lb_keough_envelopes() {
	W.recalculate_DTW_envelope();
	for (std::unique_ptr<Node> const& node : entries) {
		node->recalculate_wedge_lb_keough_envelopes();
	}
}

LeafWedgeNode::LeafWedgeNode(std::vector<double> const& timeseries, size_t M, size_t B, double r, size_t R) :
	WedgeNode(timeseries, M, B, r, R, NODE_LEAF_WEDGE)
{ }

bool LeafWedgeNode::insert_timeseries(Candidate&& C) {
	//returns true if current node's entry count exceeds B and must be split, false otherwise
	bool split_required = false;
	W.enlarge(C);
	if (entries.empty()) {
		std::unique_ptr<CandidateNode> n = std::unique_ptr<CandidateNode>(new CandidateNode(timeseries, M, B, r, R));
		n->insert_timeseries(std::forward<Candidate>(C));
		entries.push_back(std::move(n));
	}
	else {
		size_t target_index = get_min_enlargement_insertion_target(C);
		//aquire ownership of candidate node
		std::unique_ptr<Node>& target = entries[target_index];
		if (dynamic_cast<CandidateNode*>(target.get())->can_contain_candidate(C)) {
			target->insert_timeseries(std::forward<Candidate>(C));
		}
		else {
			//insertion to entry at a leaf node can be done only when the enlargement satisfies the condition W_UL <= r; otherwise
			//we have to create a new entry at that leaf node
			std::unique_ptr<CandidateNode> n = std::unique_ptr<CandidateNode>(new CandidateNode(timeseries, M, B, r, R));
			n->insert_timeseries(std::forward<Candidate>(C));
			entries.push_back(std::move(n));
		}
	}
	if (entries.size() > B) {
		//not enough space to contain the new entry in the current node, need to split
		split_required = true;
	}
	return split_required;
}

std::list<CandidateNode> LeafWedgeNode::get_merged_candidate_nodes(bool show_progress, unsigned int depth, size_t total_num_leaf_wedge_nodes, size_t* leaf_wedge_node_counter) const {
	std::list<CandidateNode> merged_nodes;
	for (std::unique_ptr<Node> const& node_ptr : entries) {
		//obtain a copy of the candidate node since we need to be able to modify it
		CandidateNode n = *dynamic_cast<CandidateNode*>(node_ptr.get()); //xxx is this the right way to do this?
		merged_nodes.push_back(std::move(n));
	}
	if (show_progress && merged_nodes.size() > 1000)
		std::cout << std::string(depth, '.') << "Starting time-consuming merge of " << merged_nodes.size() << " nodes" << std::endl;
	auto node_iter_end = merged_nodes.end();
	for (auto node_outer_iter = merged_nodes.begin(); node_outer_iter != node_iter_end; ++node_outer_iter) {
		auto node_inner_iter = node_outer_iter;
		++node_inner_iter;
		while (node_inner_iter != node_iter_end) {
			if (node_outer_iter->can_contain_wedge(node_inner_iter->get_wedge())) {
				node_outer_iter->merge_and_destroy_source(&*node_inner_iter);
				node_inner_iter = merged_nodes.erase(node_inner_iter);
			}
			else {
				++node_inner_iter;
			}
		}
	}
	++*leaf_wedge_node_counter;
	if (show_progress && (*leaf_wedge_node_counter % 100) == 0) {
		printf("\r%.2f%% completed", (100.0 * *leaf_wedge_node_counter / total_num_leaf_wedge_nodes));
		std::flush(std::cout);
	}
	return merged_nodes;
}

size_t LeafWedgeNode::count_leaf_wedge_nodes() const {
	return 1;
}

size_t LeafWedgeNode::count_candidate_nodes() const {
	return entries.size();
}

InternalWedgeNode::InternalWedgeNode(std::vector<double> const& timeseries, size_t M, size_t B, double r, size_t R) :
	WedgeNode(timeseries, M, B, r, R, NODE_INTERNAL_WEDGE)
{ }

void InternalWedgeNode::split_child(size_t target_entry_index) {
	WedgeNode* target = dynamic_cast<WedgeNode*>(entries[target_entry_index].release());
	size_t E1_index, E2_index;
	double max_union_ED = 0;
	//choose the two entries whose combined wedges have the largest euclidean distance
	std::vector<std::unique_ptr<Node>>& target_entries = target->get_entries();
	//create a flattened 2-dimensional array which holds euclidean distances between entry wedges
	std::vector<double> union_EDs(target_entries.size() * target_entries.size(), 0);
	size_t num_entries = target_entries.size();
	for (size_t i = 0; i != num_entries; ++i) {
		std::unique_ptr<Node> const& E1_tmp = target_entries[i];
		for (size_t j = i + 1; j != num_entries; ++j) {
			std::unique_ptr<Node> const& E2_tmp = target_entries[j];
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
	WedgeNode* new_entry1;
	WedgeNode* new_entry2;
	switch (target->get_node_type()) {
	case NODE_INTERNAL_WEDGE:
		new_entry1 = new InternalWedgeNode(timeseries, M, B, r, R);
		new_entry2 = new InternalWedgeNode(timeseries, M, B, r, R);
		break;
	case NODE_LEAF_WEDGE:
		new_entry1 = new LeafWedgeNode(timeseries, M, B, r, R);
		new_entry2 = new LeafWedgeNode(timeseries, M, B, r, R);
		break;
	default:
		throw new std::logic_error("Unexpected node type");
	}
	new_entry1->add_entry(std::move(target_entries[E1_index]));
	new_entry2->add_entry(std::move(target_entries[E2_index]));
	//add the rest of target's entries to either of the two new entries, based on their enlargement of the entry's original wedge
	for (size_t i = 0; i != target_entries.size(); ++i) {
		if (i != E1_index && i != E2_index) {
			std::unique_ptr<Node>& entry_to_move = target_entries[i];
			double new_entry1_ED = union_EDs[i * num_entries + E1_index];
			double new_entry2_ED = union_EDs[i * num_entries + E2_index];
			if (new_entry1_ED < new_entry2_ED) {
				new_entry1->add_entry(std::move(entry_to_move));
			}
			else {
				new_entry2->add_entry(std::move(entry_to_move));
			}
		}
	}
	//add the new entries, replacing the old one (increases entry count by one)
	//this operation does not enlarge the current wedge
	delete target; //this pointer is no longer managed and must be manually deleted
	entries[target_entry_index] = std::unique_ptr<Node>(new_entry1);
	entries.push_back(std::unique_ptr<Node>(new_entry2));
}

bool InternalWedgeNode::insert_timeseries(Candidate&& C) {
	//returns true if current node's entry count exceeds B and must be split, false otherwise
	bool split_required = false;
	W.enlarge(C);
	if (entries.empty()) {
		std::unique_ptr<LeafWedgeNode> n(new LeafWedgeNode(timeseries, M, B, r, R));
		//no need to check for if a split is necessary since the candidate is guaranteed to insertable into a fresh node
		n->insert_timeseries(std::forward<Candidate>(C));
		add_entry(std::move(n));
	}
	else {
		size_t target_index = get_min_enlargement_insertion_target(C);
		std::unique_ptr<Node>& target = entries[target_index];
		if (target->insert_timeseries(std::forward<Candidate>(C))) {
			//target experienced overflow and must be split
			split_child(target_index);
		}
		if (entries.size() > B) {
			//not enough space to contain the new entry in the current node, need to split
			split_required = true;
		}
	}
	return split_required;
}

std::list<CandidateNode> InternalWedgeNode::get_merged_candidate_nodes(bool show_progress, unsigned int depth, size_t total_num_leaf_wedge_nodes, size_t* leaf_wedge_node_counter) const {
	std::list<CandidateNode> merged_nodes;
	std::list<std::list<CandidateNode>> merged_candidate_node_lists;
	size_t nodes_to_merge = 0;
	for (std::unique_ptr<Node> const& node_ptr : entries) {
		WedgeNode* n = dynamic_cast<WedgeNode*>(node_ptr.get());
		std::list<CandidateNode> merged_descendent_nodes = n->get_merged_candidate_nodes(show_progress, depth + 1, total_num_leaf_wedge_nodes, leaf_wedge_node_counter);
		nodes_to_merge += merged_descendent_nodes.size();
		merged_candidate_node_lists.push_back(std::move(merged_descendent_nodes));
	}
	if (show_progress && nodes_to_merge > 1000)
		std::cout << std::string(depth, '.') << "Starting time-consuming merge of " << nodes_to_merge << " nodes" << std::endl;
	auto node_list_iter_end = merged_candidate_node_lists.end();
	for (auto node_list_iter = merged_candidate_node_lists.begin(); node_list_iter != node_list_iter_end; ++node_list_iter) {
		auto next_node_list_iter = node_list_iter;
		++next_node_list_iter; //all nodes have already been been compared to each other within their own list
		if (next_node_list_iter != node_list_iter_end) {
			auto node_outer_iter_end = node_list_iter->end();
			for (auto node_outer_iter = node_list_iter->begin(); node_outer_iter != node_outer_iter_end; ++node_outer_iter) {
				auto node_inner_iter = next_node_list_iter->begin();
				auto node_inner_iter_end = next_node_list_iter->end();
				while (node_inner_iter != node_inner_iter_end) {
					if (node_outer_iter->can_contain_wedge(node_inner_iter->get_wedge())) {
						node_outer_iter->merge_and_destroy_source(&*node_inner_iter);
						node_inner_iter = next_node_list_iter->erase(node_inner_iter);
					}
					else {
						++node_inner_iter;
					}
				}
			}
		}
	}
	for (std::list<CandidateNode>& merged_node_list : merged_candidate_node_lists) {
		merged_nodes.splice(merged_nodes.end(), merged_node_list);
	}
	return merged_nodes;
}

size_t InternalWedgeNode::count_leaf_wedge_nodes() const {
	size_t count = 0;
	for (std::unique_ptr<Node> const& node : entries)
		count += dynamic_cast<WedgeNode*>(node.get())->count_leaf_wedge_nodes();
	return count;
}

size_t InternalWedgeNode::count_candidate_nodes() const {
	size_t count = 0;
	for (std::unique_ptr<Node> const& node : entries)
		count += dynamic_cast<WedgeNode*>(node.get())->count_candidate_nodes();
	return count;
}
