#ifndef __MERGE_H__
#define __MERGE_H__

#include "splice_graph.h"
#include "interval_map.h"
#include "combined_graph.h"

typedef map< int32_t, set<int> > MISI;
typedef pair< int32_t, set<int> > PISI;
typedef pair<int, int> PI;
typedef pair<PI, double> PID;

class incubator
{
public:
	incubator(const string &dir);

public:
	vector<combined_graph> gset;			// graph set
	vector<bool> merged;					// whether gset[k] is merged
	string mdir;							// output dir

public:
	// binary search
	int binary_merge(const string &file);
	int binary_merge(const vector<string> &files, int low, int high, vector<combined_graph> &vc);

	int merge();
	int merge_component(const set<int> &s);
	int build_splice_map(MISI &mis);

	// write and print
	int write(const string &file, bool headers = false);
	int print();

	// analysis
	int analyze(const string &file);
};

int load(const string &file, vector<combined_graph> &vc);
bool compare_graph_overlap(const PID &x, const PID &y);

#endif
