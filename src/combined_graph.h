#ifndef __COMBINER_H__
#define __COMBINER_H__

#include "splice_graph.h"
#include "interval_map.h"

typedef pair<int32_t, int32_t> PI32;

class combiner
{
public:
	split_interval_map imap;
	map<PI32, double> emap;

	splice_graph gr;		// combined splice graph
	map<int32_t, int> lindex;
	map<int32_t, int> rindex;

public:
	int combine(const splice_graph &gt);
	int combine_vertices(const splice_graph &gt);
	int combine_edges(const splice_graph &gt);

	int build_combined_splice_graph();
	int build_vertices();
	int build_vertex_indices();
	int build_edges();

	int add_edge(splice_graph &gr, int s, int t, double w, int type);
	int draw(splice_graph &gr, const string &file);
};

#endif
