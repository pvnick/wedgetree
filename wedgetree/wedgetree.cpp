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

#include "wedgetree.h"

WedgeTree::WedgeTree(std::vector<double> const& timeseries, size_t M, size_t B, double r, bool verbose = false) :
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
	//4046418
	root->is_tree_root = true;
	size_t ts_size = timeseries.size();
	for (size_t C_index = 0; C_index + M - 1 != ts_size; ++C_index) {
		insert_timeseries(C_index);
		/*std::vector<std::unique_ptr<Node>> const& entries = root->get_entries();
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

void WedgeTree::insert_timeseries(size_t candidate_position) {
	if (root->insert_timeseries(Candidate(timeseries, candidate_position, M))) {
		//root overflow. add a new root to contain the old root, then split the old root
		std::unique_ptr<InternalWedgeNode> old_root = std::move(root);
		old_root->is_tree_root = false;
		root = std::unique_ptr<InternalWedgeNode>(new InternalWedgeNode(timeseries, M, B, r));
		root->add_entry(std::move(old_root));
		root->is_tree_root = true;
		root->split_child(0);
	}
}

std::list<CandidateNode> WedgeTree::get_merged_candidate_nodes() const {
	if (verbose) {
		std::cout << "Merging candidate nodes (subclusters)" << std::endl;
	}
	size_t total_num_candidate_nodes = root->count_candidate_nodes();
	size_t total_num_leaf_wedge_nodes = root->count_leaf_wedge_nodes();
	std::cout << total_num_candidate_nodes << std::endl;
	size_t leaf_wedge_node_counter = 0;
	std::list<CandidateNode> merged_nodes = root->get_merged_candidate_nodes(verbose, 0, total_num_leaf_wedge_nodes, &leaf_wedge_node_counter);
	if (verbose) {
		std::cout << "\r100.0% completed" << std::endl;
		std::cout << "Started with " << total_num_candidate_nodes << " candidate nodes" << std::endl;
		size_t num_merged_nodes = merged_nodes.size();
		std::cout << "Finished with " << num_merged_nodes << " candidate nodes (";			
		printf("%.2f%% reduction)\n", 100.0 * (total_num_candidate_nodes - num_merged_nodes) / total_num_candidate_nodes);
	}
	return merged_nodes;
}