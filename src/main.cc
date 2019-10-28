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
#include "interval_map.h"

using namespace std;

int main(int argc, const char **argv)
{
	if(argc == 3)
	{
		incubator icb("");
		icb.binary_merge(argv[1]);
		icb.write(argv[2]);
	}
	else if(argc == 4)
	{
		incubator icb(argv[3]);
		icb.binary_merge(argv[1]);
		icb.write(argv[2]);
	}
	return 0;
}
