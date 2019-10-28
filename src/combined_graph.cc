#include "combined_graph.h"
#include "draw.h"

combined_graph::combined_graph()
{
	num_combined = 0;
}

int combined_graph::combine(const splice_graph &gt)
{
	if(chrm == "") chrm = gt.chrm;
	assert(gt.chrm == chrm);
	combine_vertices(gt);
	combine_edges(gt);
	combine_splice_positions(gt);
	num_combined++;
	return 0;
}

int combined_graph::combine(const combined_graph &gt)
{
	if(chrm == "") chrm = gt.chrm;
	assert(gt.chrm == chrm);
	combine_vertices(gt);
	combine_edges(gt);
	combine_splice_positions(gt);
	num_combined += gt.num_combined;
	return 0;
}

int combined_graph::build_combined_splice_graph()
{
	gr.clear();
	gr.chrm = chrm;
	build_vertices();
	build_vertex_indices();
	build_edges();
	return 0;
}

int combined_graph::combine_vertices(const combined_graph &gt)
{
	for(SIMI it = gt.imap.begin(); it != gt.imap.end(); it++)
	{
		imap += make_pair(it->first, 1);
	}
	return 0;
}

int combined_graph::combine_vertices(const splice_graph &gt)
{
	for(int i = 1; i < gt.num_vertices() - 1; i++)
	{
		vertex_info vi = gt.get_vertex_info(i);
		imap += make_pair(ROI(vi.lpos, vi.rpos), 1);
		//printf("add interval: graph contains %lu vertices, %lu edges, imap.size() = %lu, advance = %ld\n", gt.num_vertices(), gt.num_edges(), imap.size(), std::distance(imap.begin(), imap.end()));
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

		int32_t ss = gt.get_vertex_info(s).rpos;
		int32_t tt = gt.get_vertex_info(t).lpos;

		if(s == 0) ss = -1;
		if(t == gt.num_vertices() - 1) tt = -2;

		PI32 p(ss, tt);
		map<PI32, DI>::iterator x = emap.find(p);

		if(x == emap.end()) emap.insert(pair<PI32, DI>(p, DI(w, 1)));
		else 
		{
			x->second.first += w;
			x->second.second += 1;
		}
	}
	return 0;
}

int combined_graph::combine_edges(const combined_graph &gt)
{
	for(map<PI32, DI>::const_iterator it = gt.emap.begin(); it != gt.emap.end(); it++)
	{
		PI32 p = it->first;
		DI d = it->second;

		map<PI32, DI>::iterator x = emap.find(p);

		if(x == emap.end())
		{
			emap.insert(pair<PI32, DI>(p, d));
		}
		else 
		{
			x->second.first += d.first;
			x->second.second += d.second;
		}
	}
	return 0;
}

int combined_graph::combine_splice_positions(const splice_graph &gt)
{
	vector<int32_t> vt = gt.get_splice_positions();
	vector<int32_t> vv(vt.size() + spos.size(), 0);
	vector<int32_t>::iterator it = set_union(vt.begin(), vt.end(), spos.begin(), spos.end(), vv.begin());
	vv.resize(it - vv.begin());
	spos = vv;
	return 0;
}

int combined_graph::combine_splice_positions(const combined_graph &gt)
{
	vector<int32_t> vv(gt.spos.size() + spos.size(), 0);
	vector<int32_t>::iterator it = set_union(gt.spos.begin(), gt.spos.end(), spos.begin(), spos.end(), vv.begin());
	vv.resize(it - vv.begin());
	spos = vv;
	return 0;
}

int combined_graph::get_overlapped_splice_positions(const vector<int32_t> &v)
{
	vector<int32_t> vv(v.size(), 0);
	vector<int32_t>::iterator it = set_intersection(v.begin(), v.end(), spos.begin(), spos.end(), vv.begin());
	return it - vv.begin();
}

PI32 combined_graph::get_bounds()
{
	if(imap.size() == 0) return PI32(-1, -1);

	SIMI it = imap.begin();
	int32_t p1 = lower(it->first);

	it = imap.end();
	it--;
	int32_t p2 = upper(it->first);

	return PI32(p1, p2);
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
	int n = gr.num_vertices() - 1;
	for(map<PI32, DI>::iterator it = emap.begin(); it != emap.end(); it++)
	{
		int32_t s = it->first.first;
		int32_t t = it->first.second;
		double w = it->second.first;
		int c = it->second.second;

		int ks = -1;
		int kt = -1;
		if(s == -1) 
		{
			ks = 0;
		}
		else
		{
			map<int32_t, int>::iterator xs = rindex.find(s);
			assert(xs != rindex.end());
			ks = xs->second;
		}

		if(t == -2)
		{
			kt = n;
		}
		else
		{
			map<int32_t, int>::iterator xt = lindex.find(t);
			assert(xt != lindex.end());
			kt = xt->second;
		}

		/*
		if(xt == lindex.end())
		{
			for(int k = 0; k < gr.num_vertices(); k++)
			{
				vertex_info vi = gr.get_vertex_info(k);
				printf("node %d: [%d, %d)\n", k, vi.lpos, vi.rpos);
			}
			printf("now test edge %d -> %d\n", s, t);
		}
		*/

		PEB p = gr.edge(ks, kt);
		assert(p.second == false);

		edge_descriptor e = gr.add_edge(ks, kt);
		gr.set_edge_weight(e, w);
		edge_info ei;
		ei.type = c;
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

int combined_graph::print(int index)
{
	PI32 p = get_bounds();
	printf("combined-graph %d: #combined = %d, chrm = %s, #intervals = %lu, #edges = %lu, boundary = [%d, %d)\n", 
			index, num_combined, chrm.c_str(), std::distance(imap.begin(), imap.end()), emap.size(), p.first, p.second);
	return 0;
}

int combined_graph::write(int index, ostream &os)
{
	char name[10240];
	sprintf(name, "graph.%d", index);
	int n = gr.num_vertices();
	int m = gr.num_edges();
	os << "# " << name << " " << chrm.c_str() << " " << n << " " << m << " " << num_combined << endl;
	gr.write(os);
	os << endl;
	return 0;
}
