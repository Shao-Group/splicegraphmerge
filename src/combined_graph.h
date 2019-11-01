#ifndef __COMBINER_H__
#define __COMBINER_H__

#include "splice_graph.h"
#include "interval_map.h"

typedef pair<int32_t, int32_t> PI32;
typedef pair<double, int> DI;

class combined_graph
{
public:
	combined_graph();

public:
	split_interval_map imap;
	map<PI32, DI> emap;
	map<vector<int32_t>, DI> pmap;
	vector<int32_t> spos;
	int num_combined;
	string chrm;

	map<int32_t, int> lindex;
	map<int32_t, int> rindex;
	splice_graph gr;			// combined splice graph
	map<vector<int>, int> hs;	// hyper-set

public:
	int combine(const combined_graph &gt);
	int combine_vertices(const combined_graph &gt);
	int combine_edges(const combined_graph &gt);
	int combine_paths(const combined_graph &gt);
	int combine_splice_positions(const combined_graph &gt);

	PI32 get_bounds();
	int get_overlapped_splice_positions(const vector<int32_t> &v);

	int build_combined_splice_graph();
	int build_vertices();
	int build_vertex_indices();
	int build_edges();
	int build_paths();

	int build(istream &is, const string &chrm);
	int write(ostream &os, int index, bool headers = false);
	int write(ostream &os);
	int print(int index);

	int analyze(int index);
};

#endif
