#ifndef __GRAPH_SET_H__
#define __GRAPH_SET_H__

#include "splice_graph.h"
#include "interval_map.h"

class graph_set
{
public:
	vector<splice_graph> gset;		// graph set
	split_interval_map imap;

public:
	int load(const string &file);
	int print();
};

#endif
