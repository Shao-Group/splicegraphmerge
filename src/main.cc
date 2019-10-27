#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <ctime>
#include <cassert>
#include <sstream>

#include "incubator.h"
#include "combined_graph.h"

using namespace std;

int main(int argc, const char **argv)
{
	vector<splice_graph> v = load(argv[1]);

	combined_graph cb;
	cb.combine(v[0]);
	cb.combine(v[1]);
	cb.build_combined_splice_graph();
	cb.gr.print_weights();

	return 0;
}
