#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <ctime>
#include <cassert>
#include <sstream>

#include "graph_set.h"
#include "interval_map.h"
#include "combined_graph.h"

using namespace std;

int main(int argc, const char **argv)
{
	/*
	graph_set gs;
	gs.load(argv[1]);
	gs.print();

	combined_graph cb;
	cb.combine(gs.gset[0]);
	cb.combine(gs.gset[1]);
	cb.build_combined_splice_graph();
	cb.gr.print_weights();
	*/

	test_interval_set_map();

	return 0;
}
