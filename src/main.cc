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

	combiner cb(gs.gset[0], gs.gset[1]);
	cb.combine();
	cb.gr3.print();

	return 0;
}
