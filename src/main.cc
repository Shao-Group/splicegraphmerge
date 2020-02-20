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
		cout<<"       " <<argv[0] << " combine <graph.list> <combined.gr> <max-combine-num> <max-threads> <merge-ratio>"<<endl;
		cout<<"       " <<argv[0] << " analyze <combined.gr>"<<endl;
		return 0;
	}

	if(string(argv[1]) == "combine")
	{
		if(argc != 7) return 0;
		incubator icb(atoi(argv[4]), atoi(argv[5]));
		icb.load(argv[2]);
		icb.merge(atoi(argv[6]));
		icb.write(argv[3]);
		//icb.binary_merge(argv[2]);
	}

	if(string(argv[1]) == "analyze")
	{
		incubator icb(0, 0);
		icb.analyze(argv[2]);
	}
	return 0;
}
