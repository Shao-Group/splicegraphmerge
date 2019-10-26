#include "combiner.h"
#include "util.h"
#include "draw.h"

combiner::combiner(const splice_graph &g1, const splice_graph &g2)
	:gr1(g1), gr2(g2)
{}

int combiner::combine()
{
	imap.clear();
	build_split_interval_map(gr1);
	build_split_interval_map(gr2);

	add_vertices();
	build_vertex_indices();

	add_inner_edges(gr1, 1);
	add_inner_edges(gr2, 2);
	add_existing_edges(gr1, 3);
	add_existing_edges(gr2, 4);

	//if(file != "") draw(gr3, file);
	return 0;
}

int combiner::build_split_interval_map(splice_graph &gr)
{
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		vertex_info vi = gr.get_vertex_info(i);
		imap += make_pair(ROI(vi.lpos, vi.rpos), 1);
	}
	return 0;
}

int combiner::add_vertices()
{
	if(imap.size() == 0) return 0;

	gr3.add_vertex();
	vertex_info vi0;
	SIMI it = imap.begin();
	vi0.lpos = lower(it->first);
	vi0.rpos = lower(it->first);
	gr3.set_vertex_info(0, vi0);

	for(it = imap.begin(); it != imap.end(); it++)
	{
		gr3.add_vertex();
		vertex_info vi;
		vi.length = upper(it->first) - lower(it->first);
		vi.lpos = lower(it->first);
		vi.rpos = upper(it->first);
		gr3.set_vertex_weight(gr3.num_vertices() - 1, 1);
		gr3.set_vertex_info(gr3.num_vertices() - 1, vi);
	}

	it = imap.end();
	it--;
	vertex_info vin;
	vin.lpos = upper(it->first);
	vin.rpos = upper(it->first);
	gr3.add_vertex();
	gr3.set_vertex_weight(gr3.num_vertices() - 1, 1);
	gr3.set_vertex_info(gr3.num_vertices() - 1, vin);
	return 0;
}

int combiner::build_vertex_indices()
{
	lindex.clear();
	rindex.clear();
	int n = gr3.num_vertices();
	for(int k = 1; k < n - 1; k++)
	{
		vertex_info vi = gr3.get_vertex_info(k);
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

	vertex_info vi1 = gr3.get_vertex_info(0);
	vertex_info vi2 = gr3.get_vertex_info(n - 1);

	assert(rindex.find(vi1.rpos) == rindex.end());
	assert(lindex.find(vi2.lpos) == lindex.end());
	rindex.insert(pair<int32_t, int>(vi1.rpos, 0));
	lindex.insert(pair<int32_t, int>(vi2.lpos, n - 1));

	//printf("add %d -> %d to rindex\n", vi1.rpos, 0);
	//printf("add %d -> %d to lindex\n", vi2.lpos, n - 1);

	return 0;
}

int combiner::add_inner_edges(splice_graph &gt, int type)
{
	for(int i = 1; i < gt.num_vertices() - 1; i++)
	{
		vertex_info vi = gt.get_vertex_info(i);
		double w = gt.get_vertex_weight(i);

		assert(lindex.find(vi.lpos) != lindex.end());
		assert(rindex.find(vi.rpos) != rindex.end());

		int k1 = lindex[vi.lpos];
		int k2 = rindex[vi.rpos];

		//int k1 = search_splice_graph(gr3, vi.lpos);
		//int k2 = search_splice_graph(gr3, vi.rpos - 1);
		assert(k1 >= 1 && k1 < gr3.num_vertices() - 1);
		assert(k2 >= 1 && k2 < gr3.num_vertices() - 1);

		for(int k = k1; k < k2; k++)
		{
			//printf("inner edge %d -> %d due to vertex [%d, %d)\n", k, k + 1, vi.lpos, vi.rpos);
			add_edge(gr3, k, k + 1, w, type);
		}
	}
	return 0;
}

int combiner::add_existing_edges(splice_graph &gt, int type)
{
	edge_iterator it1, it2;
	PEEI pei;
	for(pei = gt.edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		double w = gt.get_edge_weight(e);

		int s = e->source();
		int t = e->target();
		vertex_info vs = gt.get_vertex_info(s);
		vertex_info vt = gt.get_vertex_info(t);

		//int ss = search_splice_graph(gr3, vs.rpos - 1);
		//int tt = search_splice_graph(gr3, vt.lpos);

		assert(lindex.find(vs.rpos) != lindex.end());
		assert(rindex.find(vt.lpos) != rindex.end());

		int ss = rindex[vs.rpos];
		int tt = lindex[vt.lpos];

		//if(s == 0) ss = 0;
		//if(t == gt.num_vertices() - 1) tt = gr3.num_vertices() - 1;

		//printf("existing edge %d -> %d due to edge from %d [%d, %d) to %d [%d, %d)\n", ss, tt, s, vs.lpos, vs.rpos, t, vt.lpos, vt.rpos);
		add_edge(gr3, ss, tt, w, type);
	}
	return 0;
}

int combiner::add_edge(splice_graph &gr, int s, int t, double w, int type)
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

int combiner::search_splice_graph(splice_graph &gr, int32_t p)
{
	if(gr.num_vertices() <= 2) return -1;
	int l = 1;
	int r = gr.num_vertices() - 2;
	while(l <= r)
	{
		int m = (l + r) / 2;
		int32_t p1 = (int32_t)(gr.get_vertex_info(m).lpos);
		int32_t p2 = (int32_t)(gr.get_vertex_info(m).rpos);
		assert(p1 < p2);
		if(p >= p1 && p < p2) return m;
		if(p < p1) r = m - 1;
		if(p >= p2) l = m + 1;
	}
	return -1;
}

int combiner::compare_boundary_edges()
{
	int tp = 0, fp = 0, fn = 0;

	edge_iterator it1, it2;
	PEEI pei;
	for(pei = gr3.out_edges(0), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		double w = gr3.get_edge_weight(*it1);
		edge_info ei = gr3.get_edge_info(*it1);
		int t = (*it1)->target();
		int32_t p = gr3.get_vertex_info(t).lpos;
		assert(ei.type == 3 || ei.type == 4);
		bool b = verify_unique_5end_edge(gr3, *it1);
		if(b == true)
		{
			if(ei.type == 3)
			{
				fn++;
				printf("FN 5end boundary, weight = %.2lf, pos = %d\n", w, p);
			}
			else
			{
				fp++;
				printf("FP 5end boundary, weight = %.2lf, pos = %d\n", w, p);
			}
		}
		else
		{
			tp++;
			printf("TP 5end boundary, weight = %.2lf, pos = %d\n", w, p);
		}
	}

	for(pei = gr3.in_edges(gr3.num_vertices() - 1), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		double w = gr3.get_edge_weight(*it1);
		edge_info ei = gr3.get_edge_info(*it1);
		int s = (*it1)->source();
		int32_t p = gr3.get_vertex_info(s).rpos;
		bool b = verify_unique_3end_edge(gr3, *it1);
		if(b == true)
		{
			if(ei.type == 3)
			{
				fn++;
				printf("FN 3end boundary, weight = %.2lf, pos = %d\n", w, p);
			}
			else
			{
				fp++;
				printf("FP 3end boundary, weight = %.2lf, pos = %d\n", w, p);
			}
		}
		else
		{
			tp++;
			printf("TP 3end boundary, weight = %.2lf, pos = %d\n", w, p);
		}
	}

	printf("summary boundary edges: TP = %d FP = %d FN = %d\n", tp, fp, fn);
	return 0;
}

int combiner::compare_splice_positions()
{
	int tp = 0, fp = 0, fn = 0;
	for(int i = 1; i < gr1.num_vertices() - 1; i++)
	{
		vertex_info vi = gr1.get_vertex_info(i);
		PEB p1 = gr1.edge(0, i);
		PEB p2 = gr1.edge(i, gr1.num_vertices() - 1);

		if(p1.second == false && vi.lpos > gr1.get_vertex_info(i - 1).rpos)
		{
			int p = search_splice_graph(gr2, vi.lpos);
			if(p < 0 || gr2.get_vertex_info(p).lpos != vi.lpos)
			{
				fn++;
				printf("FN left position, length = %d, [%d, %d)\n", vi.rpos - vi.lpos, vi.lpos, vi.rpos);
			}
			else
			{
				tp++;
				printf("TP left position, length = %d, [%d, %d)\n", vi.rpos - vi.lpos, vi.lpos, vi.rpos);
			}
		}

		if(p2.second == false && vi.rpos < gr1.get_vertex_info(i + 1).lpos)
		{
			int p = search_splice_graph(gr2, vi.rpos - 1);
			if(p < 0 || gr2.get_vertex_info(p).rpos != vi.rpos)
			{
				fn++;
				printf("FN right position, length = %d, [%d, %d)\n", vi.rpos - vi.lpos, vi.lpos, vi.rpos);
			}
			else
			{
				tp++;
				printf("TP right position, length = %d, [%d, %d)\n", vi.rpos - vi.lpos, vi.lpos, vi.rpos);
			}
		}
	}

	for(int i = 1; i < gr2.num_vertices() - 1; i++)
	{
		vertex_info vi = gr2.get_vertex_info(i);
		PEB p1 = gr2.edge(0, i);
		PEB p2 = gr2.edge(i, gr2.num_vertices() - 1);

		if(p1.second == false && vi.lpos > gr2.get_vertex_info(i - 1).rpos)
		{
			int p = search_splice_graph(gr1, vi.lpos);
			if(p < 0 || gr1.get_vertex_info(p).lpos != vi.lpos)
			{
				fp++;
				printf("FP left position, length = %d, [%d, %d)\n", vi.rpos - vi.lpos, vi.lpos, vi.rpos);
			}
		}

		if(p2.second == false && vi.rpos < gr2.get_vertex_info(i + 1).lpos)
		{
			int p = search_splice_graph(gr1, vi.rpos - 1);
			if(p < 0 || gr1.get_vertex_info(p).rpos != vi.rpos)
			{
				fp++;
				printf("FP right position, length = %d, [%d, %d)\n", vi.rpos - vi.lpos, vi.lpos, vi.rpos);
			}
		}
	}

	printf("summary splice positions: TP = %d FP = %d FN = %d\n", tp, fp, fn);
	return 0;
}

bool combiner::verify_unique_5end_edge(splice_graph &gr, edge_descriptor e)
{
	int t = e->target();
	int type = gr.get_edge_info(e).type;

	// left extend
	int32_t p = gr.get_vertex_info(t).lpos;
	for(int i = t; i >= 1; i--)
	{
		vertex_info vi = gr.get_vertex_info(i);
		int32_t pp = gr.get_vertex_info(i).rpos;
		if(i < t && pp != p) break;

		edge_iterator it1, it2;
		PEEI pei;
		for(pei = gr.in_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			if((*it1)->source() != 0) continue;
			int type2 = gr.get_edge_info(*it1).type;
			if(type2 != type) return false;
		}

		p = gr.get_vertex_info(i).lpos;
	}

	// right extend
	p = gr.get_vertex_info(t).rpos;
	for(int i = t + 1; i < gr.num_vertices() - 1; i++)
	{
		vertex_info vi = gr.get_vertex_info(i);
		int32_t pp = gr.get_vertex_info(i).lpos;
		if(pp != p) break;

		edge_iterator it1, it2;
		PEEI pei;
		for(pei = gr.in_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			if((*it1)->source() != 0) continue;
			int type2 = gr.get_edge_info(*it1).type;
			if(type2 != type) return false;
		}

		p = gr.get_vertex_info(i).rpos;
	}

	return true;
}

bool combiner::verify_unique_3end_edge(splice_graph &gr, edge_descriptor e)
{
	int s = e->source();
	int type = gr.get_edge_info(e).type;

	// left extend
	int32_t p = gr.get_vertex_info(s).lpos;
	for(int i = s; i >= 1; i--)
	{
		vertex_info vi = gr.get_vertex_info(i);
		int32_t pp = gr.get_vertex_info(i).rpos;
		if(i < s && pp != p) break;

		edge_iterator it1, it2;
		PEEI pei;
		for(pei = gr.out_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			if((*it1)->target() != gr.num_vertices() - 1) continue;
			int type2 = gr.get_edge_info(*it1).type;
			if(type2 != type) return false;
		}

		p = gr.get_vertex_info(i).lpos;
	}

	// right extend
	p = gr.get_vertex_info(s).rpos;
	for(int i = s + 1; i < gr.num_vertices() - 1; i++)
	{
		vertex_info vi = gr.get_vertex_info(i);
		int32_t pp = gr.get_vertex_info(i).lpos;
		if(pp != p) break;

		edge_iterator it1, it2;
		PEEI pei;
		for(pei = gr.out_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			if((*it1)->target() != gr.num_vertices() - 1) continue;
			int type2 = gr.get_edge_info(*it1).type;
			if(type2 != type) return false;
		}

		p = gr.get_vertex_info(i).rpos;
	}

	return true;
}

int combiner::draw(splice_graph &gr, const string &file)
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
