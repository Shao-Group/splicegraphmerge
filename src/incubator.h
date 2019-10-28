#ifndef __MERGE_H__
#define __MERGE_H__

#include "splice_graph.h"
#include "interval_map.h"
#include "combined_graph.h"

class incubator
{
public:
	vector<combined_graph> gset;		// graph set
	interval_set_map ism;
	vector<bool> merged;

public:
	int merge(const string &file);
	int merge(const splice_graph &gr);
	int merge_final();
	int merge_component(const set<int> &s);
	int write(const string &file);
	int print();
};

vector<splice_graph> load(const string &file);

#endif
