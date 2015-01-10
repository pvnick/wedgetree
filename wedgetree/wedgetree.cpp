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
	size_t ts_size = timeseries.size();
	for (size_t C_index = 0; C_index + M - 1 != ts_size; ++C_index) {
		insert_timeseries(Candidate(timeseries, C_index, M));
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

void WedgeTree::insert_timeseries(Candidate&& C) {
	if (root->insert_timeseries(std::forward<Candidate>(C))) {
		//root overflow. add a new root to contain the old root, then split the old root
		std::shared_ptr<InternalWedgeNode> old_root = root;
		root = std::make_shared<InternalWedgeNode>(timeseries, M, B, r);
		root->add_entry(old_root);
		root->split_child(0);
	}
}