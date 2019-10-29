#ifndef __MERGE_H__
#define __MERGE_H__

#include "splice_graph.h"
#include "interval_map.h"
#include "combined_graph.h"

class incubator
{
public:
	incubator(const string &dir);

public:
	vector<combined_graph> gset;			// graph set
	map< int32_t, set<int> > mis;			// splice map
	interval_set_map ism;
	vector<bool> merged;
	string mdir;

public:
	// single-chain merge
	int merge(const string &file);
	int merge(const splice_graph &gr);
	int merge_final_interval_map();
	int merge_final_splice_map();
	int merge_component(const set<int> &s);

	// binary search
	int binary_merge(const string &file);
	int binary_merge(const vector<string> &files, int low, int high, vector<combined_graph> &vc);
	int build_interval_map();
	int build_splice_map();

	// write and print
	int write(const string &file);
	int print();
};

vector<splice_graph> load(const string &file);

#endif
