#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <ctime>
#include <cassert>
#include <sstream>

#include "graph_set.h"

using namespace std;

int main(int argc, const char **argv)
{
	graph_set gs;
	gs.load(argv[1]);
	gs.print();

	return 0;
}
