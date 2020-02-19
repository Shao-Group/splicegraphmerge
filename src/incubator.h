#ifndef __MERGE_H__
#define __MERGE_H__

#include "splice_graph.h"
#include "interval_map.h"
#include "combined_graph.h"
#include "combined_group.h"

typedef map< int32_t, set<int> > MISI;
typedef pair< int32_t, set<int> > PISI;
typedef pair<int, int> PI;
typedef pair<PI, double> PID;

class incubator
{
public:
	incubator(int m, int t, const string &dir);

public:
	vector<combined_group> groups;			// graph groups

	vector<combined_graph> fixed;			// fixed set of graphs
	int max_combined_num;					// parameter
	string mdir;							// output dir
	int max_threads;

public:
	// multiple-thread load
	int load(const string &file);

	// binary search
	int binary_merge(const string &file);
	int binary_merge(const vector<string> &files, int low, int high, vector<combined_graph> &vc, bool last);
	int merge(const vector<combined_graph> &grset, vector<combined_graph> &vc, bool last);
	int build_splice_map(const vector<combined_graph> &grset, MISI &mis);

	// write and print
	int write(const string &file, bool headers = false);
	int print();

	// analysis
	int analyze(const string &file);
};

int load_multiple(const vector<string> &files, vector<combined_graph> &vc);
int load_single(const string &file, vector<combined_graph> &vc);
bool compare_graph_overlap(const PID &x, const PID &y);

#endif
