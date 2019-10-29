#ifndef __MERGE_H__
#define __MERGE_H__

#include "splice_graph.h"
#include "interval_map.h"
#include "combined_graph.h"

typedef map< int32_t, set<int> > MISI;
typedef pair< int32_t, set<int> > PISI;

class incubator
{
public:
	incubator(const string &dir);

public:
	vector<combined_graph> gset;			// graph set
	interval_set_map ism;					// interval map
	MISI mis;								// splice map
	vector<bool> merged;					// whether gset[k] is merged
	string mdir;							// output dir

public:
	// binary search
	int binary_merge(const string &file);
	int binary_merge(const vector<string> &files, int low, int high, vector<combined_graph> &vc);

	int merge();
	int merge_component(const set<int> &s);

	int build_interval_map();
	int build_splice_map();

	// write and print
	int write(const string &file);
	int print();
};

int load(const string &file, vector<splice_graph> &vs);
int load(const string &file, vector<combined_graph> &vc);

#endif
