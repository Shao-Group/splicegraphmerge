#ifndef __COMBINED_GROUP_H__
#define __COMBINED_GROUP_H__

#include "combined_graph.h"

class combined_group
{
public:
	combined_group();

public:
	vector<combined_graph> gset;
	string chrm;
	char strand;
};

#endif
