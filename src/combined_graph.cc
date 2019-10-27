#include "combined_graph.h"
#include "util.h"
#include "draw.h"

int combined_graph::combine(const splice_graph &gt)
{
	combine_vertices(gt);
	combine_edges(gt);
	return 0;
}

int combined_graph::build_combined_splice_graph()
{
	build_vertices();
	build_vertex_indices();
	build_edges();
	return 0;
}

int combined_graph::combine_vertices(const splice_graph &gt)
{
	for(int i = 1; i < gt.num_vertices() - 1; i++)
	{
		vertex_info vi = gt.get_vertex_info(i);
		imap += make_pair(ROI(vi.lpos, vi.rpos), 1);
	}
	return 0;
}

int combined_graph::combine_edges(const splice_graph &gt)
{
	PEEI pei = gt.edges();
	for(edge_iterator it = pei.first; it != pei.second; it++)
	{
		edge_descriptor e = (*it);
		double w = gt.get_edge_weight(e);

		int s = e->source();
		int t = e->target();
		vertex_info vs = gt.get_vertex_info(s);
		vertex_info vt = gt.get_vertex_info(t);

		PI32 p(vs.rpos, vt.lpos);
		map<PI32, double>::iterator x = emap.find(p);

		if(x == emap.end()) emap.insert(pair<PI32, double>(p, w));
		else x->second += 1;
	}
	return 0;
}

int combined_graph::build_vertices()
{
	if(imap.size() == 0) return 0;

	gr.add_vertex();
	vertex_info vi0;
	SIMI it = imap.begin();
	vi0.lpos = lower(it->first);
	vi0.rpos = lower(it->first);
	gr.set_vertex_info(0, vi0);

	for(it = imap.begin(); it != imap.end(); it++)
	{
		gr.add_vertex();
		vertex_info vi;
		vi.length = upper(it->first) - lower(it->first);
		vi.lpos = lower(it->first);
		vi.rpos = upper(it->first);
		gr.set_vertex_weight(gr.num_vertices() - 1, 1);
		gr.set_vertex_info(gr.num_vertices() - 1, vi);
	}

	it = imap.end();
	it--;
	vertex_info vin;
	vin.lpos = upper(it->first);
	vin.rpos = upper(it->first);
	gr.add_vertex();
	gr.set_vertex_weight(gr.num_vertices() - 1, 1);
	gr.set_vertex_info(gr.num_vertices() - 1, vin);
	return 0;
}

int combined_graph::build_vertex_indices()
{
	lindex.clear();
	rindex.clear();
	int n = gr.num_vertices();
	for(int k = 1; k < n - 1; k++)
	{
		vertex_info vi = gr.get_vertex_info(k);
		int32_t l = vi.lpos;
		int32_t r = vi.rpos;
		assert(l < r);
		assert(lindex.find(l) == lindex.end());
		assert(rindex.find(r) == rindex.end());
		lindex.insert(pair<int32_t, int>(l, k));
		rindex.insert(pair<int32_t, int>(r, k));

		//printf("add %d -> %d to lindex\n", l, k);
		//printf("add %d -> %d to rindex\n", r, k);
	}

	vertex_info vi1 = gr.get_vertex_info(0);
	vertex_info vi2 = gr.get_vertex_info(n - 1);

	assert(rindex.find(vi1.rpos) == rindex.end());
	assert(lindex.find(vi2.lpos) == lindex.end());
	rindex.insert(pair<int32_t, int>(vi1.rpos, 0));
	lindex.insert(pair<int32_t, int>(vi2.lpos, n - 1));

	//printf("add %d -> %d to rindex\n", vi1.rpos, 0);
	//printf("add %d -> %d to lindex\n", vi2.lpos, n - 1);
	return 0;
}

int combined_graph::build_edges()
{
	for(map<PI32, double>::iterator it = emap.begin(); it != emap.end(); it++)
	{
		int32_t s = it->first.first;
		int32_t t = it->first.second;
		double w = it->second;

		map<int32_t, int>::iterator xs = rindex.find(s);
		map<int32_t, int>::iterator xt = lindex.find(t);
		assert(xs != rindex.end());
		assert(xt != lindex.end());

		int ks = xs->second;
		int kt = xt->second;

		PEB p = gr.edge(ks, kt);
		assert(p.second == false);

		edge_descriptor e = gr.add_edge(ks, kt);
		gr.set_edge_weight(e, w);
		edge_info ei;
		gr.set_edge_info(e, ei);
	}

	return 0;
}

int combined_graph::add_edge(splice_graph &gr, int s, int t, double w, int type)
{
	assert(s >= 0 && s < gr.num_vertices());
	assert(t >= 0 && t < gr.num_vertices());

	PEB p = gr.edge(s, t);
	if(p.second == false)
	{
		edge_descriptor e = gr.add_edge(s, t);
		gr.set_edge_weight(e, w);
		edge_info ei;
		ei.type = type;
		gr.set_edge_info(e, ei);
	}
	else
	{
		edge_descriptor e = p.first;
		w += gr.get_edge_weight(e);
		gr.set_edge_weight(e, w);
		edge_info ei;
		ei.type = type;
		gr.set_edge_info(e, ei);
	}
	return 0;
}

int combined_graph::draw(splice_graph &gr, const string &file)
{
	ofstream fout(file.c_str());
	if(fout.fail())
	{
		printf("open file %s error.\n", file.c_str());
		return 0;
	}

	draw_header(fout);

	double len = 3.0;

	fout<<"\\def\\len{"<<len<<"cm}\n";

	// draw file name
	fout<<"\\node[draw, thick, red] at (1.6 * \\len, 0.58 * \\len) {"<<file.c_str()<<"};\n";

	// draw vertices
	char sx[1024];
	char sy[1024];
	double pos = 0;
	for(int i = 0; i < gr.num_vertices(); i++)
	{
		int d = gr.degree(i);
		//if(d == 0) continue;

		vertex_info vi = gr.get_vertex_info(i);

		pos++;

		sprintf(sx, "s%d", i);
		string s = "";
		fout.precision(0);
		fout<<fixed;
		fout<<"\\node[mycircle, \\colx, draw, label = below:{";
		fout<<vi.lpos % 100000<<"-"<<vi.rpos % 100000;
		fout<<"}] ("<<sx<<") at ("<<pos<<" *\\len, 0.0) {"<<i<<"};\n";
	}

	// draw reference edges
	edge_iterator it1, it2;
	PEEI pei;
	for(pei = gr.edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		edge_info ei = gr.get_edge_info(e);

		int d = (int)(gr.get_edge_weight(e) * 2.0);

		if(ei.type == 2 || ei.type == 4) continue;

		int s = e->source();
		int t = e->target();

		sprintf(sx, "s%d", s);
		sprintf(sy, "s%d", t);

		double bend = 0;
		if(ei.type == 3) bend = -40;

		string line = "line width = 0.12cm, gray, ";
		if(e->source() == 0 || e->target() == gr.num_vertices() - 1) line = "line width = 0.12cm, gray, densely dotted, ";

		fout<<"\\draw[->,"<< line.c_str() <<"bend right = "<< bend <<"] ("<<sx<<") to node[gray, label=below:{" << d << "}]{} "<<"("<<sy<<");\n";
	}

	// draw evaluated edges
	for(pei = gr.edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		edge_info ei = gr.get_edge_info(e);

		int d = (int)(gr.get_edge_weight(e));

		if(ei.type == 1 || ei.type == 3) continue;

		int s = e->source();
		int t = e->target();

		sprintf(sx, "s%d", s);
		sprintf(sy, "s%d", t);

		double bend = 0;
		if(ei.type == 4) bend = -40;

		string line = "line width = 0.02cm, red,";
		if(e->source() == 0 || e->target() == gr.num_vertices() - 1) line = "line width = 0.02cm, red, densely dotted, ";

		fout<<"\\draw[->,"<< line.c_str() <<"bend right = "<< bend <<"] ("<<sx<<") to node[red, label=above:{" << d << "}]{} "<<"("<<sy<<");\n";
	}

	draw_footer(fout);

	fout.close();
	return 0;
}