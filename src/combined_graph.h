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
	vector<int32_t> spos;
	int num_combined;
	string chrm;

	splice_graph gr;		// combined splice graph
	map<int32_t, int> lindex;
	map<int32_t, int> rindex;

public:
	int combine(const splice_graph &gt);
	int combine_vertices(const splice_graph &gt);
	int combine_edges(const splice_graph &gt);
	int combine_splice_positions(const splice_graph &gt);

	int build_combined_splice_graph();
	int build_vertices();
	int build_vertex_indices();
	int build_edges();

	PI32 get_bounds();
	int get_overlapped_splice_positions(const vector<int32_t> &v);
	int add_edge(splice_graph &gr, int s, int t, double w, int type);
	int draw(splice_graph &gr, const string &file);

	int print(int index);
	int write(int index, ostream &os);
};

#endif
