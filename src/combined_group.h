#ifndef __COMBINED_GROUP_H__
#define __COMBINED_GROUP_H__

#include "combined_graph.h"

class combined_group
{
public:
	combined_group(string c, char s);

public:
	vector<combined_graph> gset;
	string chrm;
	char strand;

public:
	int add_graph(const combined_graph &gr);
};

#endif
