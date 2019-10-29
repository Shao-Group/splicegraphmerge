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
	if(argc == 1)
	{
		cout<<"usage: " <<endl;
		cout<<"       " <<argv[0] << " combine <graph.list> <combined.gr> [<intermediate-dir>]"<<endl;
		cout<<"       " <<argv[0] << " analyze <combined.gr>"<<endl;
		return 0;
	}

	if(string(argv[1]) == "combine")
	{
		if(argc == 4)
		{
			incubator icb("");
			icb.binary_merge(argv[2]);
			icb.write(argv[3]);
		}
		else if(argc == 5)
		{
			incubator icb(argv[4]);
			icb.binary_merge(argv[2]);
			icb.write(argv[3]);
		}
	}

	if(string(argv[1]) == "analyze")
	{
		incubator icb("");
		icb.analyze(argv[2]);
	}
	return 0;
}
