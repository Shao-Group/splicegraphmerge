#include "combined_graph.h"
#include "draw.h"

combined_graph::combined_graph()
{
	num_combined = 0;
}

int combined_graph::combine(const combined_graph &gt)
{
	if(chrm == "") chrm = gt.chrm;
	assert(gt.chrm == chrm);
	combine_vertices(gt);
	combine_edges(gt);
	combine_paths(gt);
	combine_splice_positions(gt);
	num_combined += gt.num_combined;
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

int combined_graph::combine_paths(const combined_graph &gt)
{
	for(map<vector<int32_t>, DI>::const_iterator it = gt.pmap.begin(); it != gt.pmap.end(); it++)
	{
		const vector<int32_t> &p = it->first;
		DI d = it->second;

		map<vector<int32_t>, DI>::iterator x = pmap.find(p);

		if(x == pmap.end())
		{
			pmap.insert(pair<vector<int32_t>, DI>(p, d));
		}
		else 
		{
			x->second.first += d.first;
			x->second.second += d.second;
		}
	}
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

int combined_graph::build_combined_splice_graph()
{
	gr.clear();
	gr.chrm = chrm;
	build_vertices();
	build_vertex_indices();
	build_edges();
	build_paths();
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

int combined_graph::build_paths()
{
	int n = gr.num_vertices() - 1;
	for(map<vector<int32_t>, DI>::iterator it = pmap.begin(); it != pmap.end(); it++)
	{
		const vector<int32_t> &v = it->first;
		double w = it->second.first;
		int c = it->second.second;

		vector<int> vv;
		bool fail = false;
		for(int k = 0; k < v.size() / 2; k++)
		{
			int32_t s = v[2 * k + 0];
			int32_t t = v[2 * k + 1];

			int ks = -1;
			int kt = -1;
			if(s == -1) 
			{
				fail = true;
				ks = 0;
			}
			else
			{
				map<int32_t, int>::iterator xs = rindex.find(s);
				if(xs == rindex.end()) fail = true;
				else ks = xs->second;
			}

			if(fail == true) break;

			if(t == -2)
			{
				fail = true;
				kt = n;
			}
			else
			{
				map<int32_t, int>::iterator xt = lindex.find(t);
				if(xt == lindex.end()) fail = true;
				else kt = xt->second;
			}

			if(fail == true) break;

			if(vv.size() >= 1)
			{
				int z = vv.back();
				for(int j = z + 1; j < ks; j++) vv.push_back(j);
			}
			vv.push_back(ks);
			vv.push_back(kt);
		}

		assert(hs.find(vv) == hs.end());
		hs.insert(pair<vector<int>, int>(vv, (int)w));
	}

	return 0;
}


int combined_graph::build(istream &is, const string &c)
{
	chrm = c;

	char line[10240];
	char name[10240];

	int n = -1;
	vector<int32_t> vv1;
	vector<int32_t> vv2;
	while(is.getline(line, 10240, '\n'))
	{
		stringstream sstr(line);
		if(string(line).length() == 0) break;
		
		sstr >> name;

		if(string(name) == "node")
		{
			int index;
			double weight;
			int32_t lpos;
			int32_t rpos;
			sstr >> index >> weight >> lpos >> rpos;
			imap += make_pair(ROI(lpos, rpos), 1);

			if(index > n) n = index;
			if(vv1.size() <= n) vv1.resize(n + 1);
			if(vv2.size() <= n) vv2.resize(n + 1);

			vv1[index] = lpos;
			vv2[index] = rpos;
		}
		else if(string(name) == "edge")
		{
			int x, y;
			double w;
			int c;
			sstr >> x >> y >> w >> c;

			assert(x != y);
			assert(x >= 0 && x <= n);
			assert(y >= 0 && y <= n);

			int32_t s = vv2[x];
			int32_t t = vv1[y];

			if(x == 0) s = -1;
			if(y == n) t = -2;

			if(x != 0 && y != n && s < t) spos.push_back(s);
			if(x != 0 && y != n && s < t) spos.push_back(t);

			PI32 p(s, t);
			map<PI32, DI>::iterator it = emap.find(p);

			if(it == emap.end()) 
			{
				emap.insert(pair<PI32, DI>(p, DI(w, c)));
			}
			else 
			{
				it->second.first += w;
				it->second.second += c;
			}
		}
		else if(string(name) == "path" || string(name) == "topo")
		{
			vector<int32_t> v;
			int z;
			sstr >> z;
			assert(z >= 1);
			int x, y;
			sstr >> x;
			for(int k = 1; k < z; k++)
			{
				sstr >> y;

				assert(x != y);
				assert(x >= 0 && x <= n);
				assert(y >= 0 && y <= n);

				int32_t s = vv2[x];
				int32_t t = vv1[y];

				if(x == 0) s = -1;
				if(y == n) t = -2;

				//if(x != 0 && y != n && s < t) spos.push_back(s);
				//if(x != 0 && y != n && s < t) spos.push_back(t);
				v.push_back(s);
				v.push_back(t);
			}

			double w;
			int c;

			sstr >> w >> c;

			map<vector<int32_t>, DI>::iterator it = pmap.find(v);

			if(it == pmap.end()) 
			{
				if(string(name) == "topo") pmap.insert(pair<vector<int32_t>, DI>(v, DI(w, 0)));
				if(string(name) == "path") pmap.insert(pair<vector<int32_t>, DI>(v, DI(w, c)));
			}
			else 
			{
				it->second.first += w;
				it->second.second += c;
			}
		}
		else
		{
			break;
		}
	}

	sort(spos.begin(), spos.end());
	num_combined++;

	return 0;
}

int combined_graph::write(ostream &os)
{
	os<<fixed;
	os.precision(2);

	PI32 p = get_bounds();

	lindex.clear();
	rindex.clear();

	int id = 0;

	rindex.insert(pair<int32_t, int>(p.first, 0));
	os<<"node " << id++ << " "<< 1 << " " << p.first << " " << p.first << endl;

	for(SIMI it = imap.begin(); it != imap.end(); it++)
	{
		int32_t l = lower(it->first);
		int32_t r = upper(it->first);
		lindex.insert(pair<int32_t, int>(l, id));
		rindex.insert(pair<int32_t, int>(r, id));
		os << "node " << id++ << " "<< it->second << " ";
		os << l << " " << r << endl;
	}

	lindex.insert(pair<int32_t, int>(p.second, id));
	os<<"node " << id << " "<< 1 << " " << p.second << " " << p.second << endl;

    assert(std::distance(imap.begin(), imap.end()) + 1 == id);

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
			kt = id;
		}
		else
		{
			map<int32_t, int>::iterator xt = lindex.find(t);
			assert(xt != lindex.end());
			kt = xt->second;
		}

		os << "edge " << ks << " " << kt << " " << w << " " << c << endl;
	}

	for(map<vector<int32_t>, DI>::iterator it = pmap.begin(); it != pmap.end(); it++)
	{
		const vector<int32_t> &v = it->first;
		double w = it->second.first;
		int c = it->second.second;

		vector<int> vv;
		bool fail = false;
		for(int k = 0; k < v.size() / 2; k++)
		{
			int32_t s = v[2 * k + 0];
			int32_t t = v[2 * k + 1];

			int ks = -1;
			int kt = -1;
			if(s == -1) 
			{
				fail = true;
				ks = 0;
			}
			else
			{
				map<int32_t, int>::iterator xs = rindex.find(s);
				if(xs == rindex.end()) fail = true;
				else ks = xs->second;
			}

			if(fail == true) break;

			if(t == -2)
			{
				fail = true;
				kt = id;
			}
			else
			{
				map<int32_t, int>::iterator xt = lindex.find(t);
				if(xt == lindex.end()) fail = true;
				else kt = xt->second;
			}

			if(fail == true) break;

			if(vv.size() >= 1)
			{
				int z = vv.back();
				for(int j = z + 1; j <= ks; j++) vv.push_back(j);
			}
			vv.push_back(kt);
		}
		if(fail == true) continue;
		os << "path " << vv.size();
		for(int i = 0; i < vv.size(); i++) os << " " << vv[i];
		os << " " << w << " " << c << endl;
	}
	return 0;
}

int combined_graph::write(ostream &os, int index, bool headers)
{
	char name[10240];
	sprintf(name, "graph.%d", index);
	int n = std::distance(imap.begin(), imap.end()) + 2;
	int m = emap.size();
	os << "# " << name << " " << chrm.c_str() << " " << n << " " << m << " " << num_combined << endl;
	if(headers == false)
	{
		write(os);
		os << endl;
	}
	return 0;
}

int combined_graph::print(int index)
{
	PI32 p = get_bounds();
	printf("combined-graph %d: #combined = %d, chrm = %s, #intervals = %lu, #edges = %lu, boundary = [%d, %d)\n", 
			index, num_combined, chrm.c_str(), std::distance(imap.begin(), imap.end()), emap.size(), p.first, p.second);
	return 0;
}

int combined_graph::analyze(int index)
{
	int num_junctions = 0;
	int total_support = 0;

	PI32 p = get_bounds();
	for(map<PI32, DI>::iterator it = emap.begin(); it != emap.end(); it++)
	{
		int32_t s = it->first.first;
		int32_t t = it->first.second;

		if(s == p.first) continue;
		if(t == p.second) continue;

		if(s >= t) continue;

		double w = it->second.first;
		int c = it->second.second;

		num_junctions++;
		total_support += c;
	}

	printf("analyze-graph %d: %d junctions, %d total-support, %.2lf average-support, chrm = %s, #vertices = %lu, #edges = %lu, boundary = [%d, %d)\n", 
			index, num_junctions, total_support, total_support * 1.0 / num_junctions, chrm.c_str(), std::distance(imap.begin(), imap.end()), emap.size(), p.first, p.second);

	return 0;
}
