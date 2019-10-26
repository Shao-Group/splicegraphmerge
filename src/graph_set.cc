#include "graph_set.h"
#include <fstream>
#include <sstream>
#include <iostream>

int graph_set::load(const string &file)
{
	ifstream fin(file.c_str());
	if(fin.fail())
	{
		printf("could not load file %s\n", file.c_str());
		exit(0);
	}

	char line[10240];
	char gid[10240];
	char chrm[10240];
	char tmp[1024];

	while(fin.getline(line, 10240, '\n'))
	{
		if(line[0] != '#') continue;
		stringstream sstr(line);
		sstr >> tmp >> gid >> chrm;

		splice_graph gr;
		gr.build(fin, gid, chrm);
		gset.push_back(gr);
	}
	
	fin.close();
	return 0;
}

int graph_set::print()
{
	for(int k = 0; k < gset.size(); k++)
	{
		printf("graph %d contains %lu vertices and %lu edges\n", k, gset[k].num_vertices(), gset[k].num_edges());
	}
	return 0;
}
