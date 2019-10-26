#ifndef __COMBINER_H__
#define __COMBINER_H__

#include "splice_graph.h"
#include "interval_map.h"

class combiner
{
public:
	combiner(const splice_graph &g1, const splice_graph &g2);

public:
	splice_graph gr1;		// reference sgraph
	splice_graph gr2;		// evaluate sgraph
	splice_graph gr3;		// combination

public:
	split_interval_map imap;
	map<int32_t, int> lindex;
	map<int32_t, int> rindex;

public:
	int combine();

private:
	int build_split_interval_map(splice_graph &gr);
	int add_vertices();
	int build_vertex_indices();
	int add_inner_edges(splice_graph &gt, int type);
	int add_existing_edges(splice_graph &gt, int type);
	int add_edge(splice_graph &gr, int s, int t, double w, int type);
	int search_splice_graph(splice_graph &gr, int32_t p);

	int compare_splice_positions();
	int compare_boundary_edges();
	bool verify_unique_5end_edge(splice_graph &gr, edge_descriptor e);
	bool verify_unique_3end_edge(splice_graph &gr, edge_descriptor e);
	int draw(splice_graph &gr, const string &file);
};

#endif
