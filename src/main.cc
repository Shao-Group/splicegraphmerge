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
	incubator icb;
	icb.merge(argv[1]);
	return 0;
}
