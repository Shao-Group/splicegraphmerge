#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <ctime>
#include <cassert>
#include <sstream>

#include "graph_set.h"
#include "combiner.h"

using namespace std;

int main(int argc, const char **argv)
{
	graph_set gs;
	gs.load(argv[1]);
	gs.print();

	combiner cb;
	cb.combine(gs.gset[0]);
	cb.combine(gs.gset[1]);
	cb.build_combined_splice_graph();
	cb.gr.print_weights();

	return 0;
}
