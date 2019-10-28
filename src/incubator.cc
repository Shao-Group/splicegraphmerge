#include "incubator.h"
#include "undirected_graph.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

int incubator::merge(const string &file)
{
	ifstream fin(file.c_str());
	if(fin.fail())
	{
		printf("cannot open file %s\n", file.c_str());
		exit(0);
	}

	char line[102400];
	while(fin.getline(line, 10240, '\n'))
	{
		string s(line);
		if(s.size() == 0) continue;

		vector<splice_graph> v = load(s);
		for(int k = 0; k < v.size(); k++)
		{
			merge(v[k]);
		}

		printf("after combining %s\n", line);
		print();
		printf("\n");
	}

	merge_final();

	return 0;
}

int incubator::binary_merge(const string &file)
{
	ifstream fin(file.c_str());
	if(fin.fail())
	{
		printf("cannot open file %s\n", file.c_str());
		exit(0);
	}

	vector<string> files;

	char line[102400];
	while(fin.getline(line, 10240, '\n'))
	{
		string s(line);
		if(s.size() == 0) continue;
		files.push_back(s);
	}

	vector<combined_graph> vc;
	binary_merge(files, 0, files.size(), vc);

	gset = vc;
	merged.assign(gset.size(), false);

	return 0;
}

int incubator::binary_merge(const vector<string> &files, int low, int high, vector<combined_graph> &vc)
{
	//printf("binary merge from %d to %d (total = %lu)\n", low, high, files.size());
	vc.clear();

	if(low >= high) return 0;

	if(low + 1 == high)
	{
		string file = files[low];
		vector<splice_graph> v = load(file);
		for(int k = 0; k < v.size(); k++)
		{
			combined_graph gr;
			gr.combine(v[k]);
			vc.push_back(gr);
		}
		printf("create %lu combined-graphs for file %d (%s)\n", vc.size(), low, files[low].c_str());
		return 0;
	}

	int mid = (low + high) / 2;

	vector<combined_graph> vc1;
	vector<combined_graph> vc2;

	binary_merge(files, low, mid, vc1);
	binary_merge(files, mid, high, vc2);

	gset = vc1;
	gset.insert(gset.end(), vc2.begin(), vc2.end());

	build_interval_set_map();

	printf("merge final with %lu (%lu/%lu) combined-graphs for files [%d, %d)\n", gset.size(), vc1.size(), vc2.size(), low, high);

	merge_final();

	assert(merged.size() == gset.size());

	for(int k = 0; k < gset.size(); k++)
	{
		if(merged[k] == true) continue;
		vc.push_back(gset[k]);
	}
	return 0;
}

int incubator::merge_final()
{
	undirected_graph gr;
	for(int k = 0; k < gset.size(); k++) gr.add_vertex();

	int min_overlapped_splice_position = 2;
	for(ISMI it = ism.begin(); it != ism.end(); it++)
	{
		set<int> s = it->second;
		vector<int> v(s.begin(), s.end());
		for(int xi = 0; xi < v.size(); xi++)
		{
			int i = v[xi];
			for(int xj = xi + 1; xj < v.size(); xj++)
			{
				int j = v[xj];
				int c = gset[j].get_overlapped_splice_positions(gset[i].spos);
				if(c < min_overlapped_splice_position) continue;
				if(gset[i].chrm != gset[j].chrm) continue;
				gr.add_edge(i, j);
			}
		}
	}

	vector< set<int> > vs = gr.compute_connected_components();

	merged.resize(gset.size());
	merged.assign(gset.size(), false);


	for(int k = 0; k < vs.size(); k++)
	{
		merge_component(vs[k]);
	}
	return 0;
}

int incubator::merge_component(const set<int> &s)
{
	if(s.size() <= 1) return 0;
	int x = *(s.begin());

	for(set<int>::iterator it = s.begin(); it != s.end(); it++)
	{
		int k = *it;
		if(k == x) continue;

		//printf("final combine %d and %d with %lu and %lu vertices\n", x, k, gset[x].imap.size(), gset[k].imap.size());

		gset[x].combine(gset[k]);
		merged[k] = true;
	}
	return 0;
}

int incubator::write(const string &file)
{
	ofstream fout(file.c_str());
	if(fout.fail()) exit(1);

	for(int k = 0; k < gset.size(); k++)
	{
		if(merged[k] == true) continue;
		gset[k].build_combined_splice_graph();
		gset[k].write(k, fout);
	}
	fout.close();
	return 0;
}

int incubator::merge(const splice_graph &gr)
{
	int n = gr.num_vertices() - 1;
	int32_t l = gr.get_vertex_info(0).rpos;
	int32_t r = gr.get_vertex_info(n).lpos;

	//print_interval_set_map(ism);
	//printf("query [%d, %d) in ism\n", l, r);

	set<int> s = get_overlapped_set_partial(ism, l, r);

	//printf("get %lu overlapped\n", s.size());

	vector<int32_t> spos = gr.get_splice_positions();

	// TODO parameter
	int min_overlapped_splice_position = 2;

	set<int> ss;
	for(set<int>::iterator it = s.begin(); it != s.end(); it++)
	{
		int k = *it;
		combined_graph &csg = gset[k];
		if(csg.chrm != gr.chrm) continue;
		int overlap = csg.get_overlapped_splice_positions(spos);
		if(overlap < min_overlapped_splice_position) continue;
		ss.insert(k);
		break;
	}

	if(ss.size() == 0)
	{
		combined_graph csg;
		gset.push_back(csg);
		ss.insert(gset.size() - 1);
	}

	for(set<int>::iterator it = ss.begin(); it != ss.end(); it++)
	{
		int k = *it;
		combined_graph &csg = gset[k];
		PI32 p = csg.get_bounds();

		csg.combine(gr);

		set<int> x;
		x.insert(k);

		if(p.first == -1 || p.second == -1)
		{
			assert(l < r);
			ism += make_pair(interval32(l, r), x);
		}
		else
		{
			if(l < p.first) ism += make_pair(interval32(l, p.first), x);
			if(r > p.second) ism += make_pair(interval32(p.second, r), x);
		}
	}
	return 0;
}

int incubator::build_interval_set_map()
{
	ism.clear();
	for(int k = 0; k < gset.size(); k++)
	{
		PI32 p = gset[k].get_bounds();
		set<int> s;
		s.insert(k);
		ism += make_pair(interval32(p.first, p.second), s); 
	}
	return 0;
}

int incubator::print()
{
	for(int k = 0; k < gset.size(); k++) gset[k].print(k);
	return 0;
}

vector<splice_graph> load(const string &file)
{
	ifstream fin(file.c_str());
	if(fin.fail())
	{
		printf("could not load file %s\n", file.c_str());
		exit(0);
	}

	char line[10240];
	char gid[10240];
	char chrm[10240];
	char tmp[1024];

	vector<splice_graph> v;
	while(fin.getline(line, 10240, '\n'))
	{
		if(line[0] != '#') continue;
		stringstream sstr(line);
		sstr >> tmp >> gid >> chrm;

		splice_graph gr;
		gr.build(fin, gid, chrm);
		v.push_back(gr);
	}
	
	fin.close();
	return v;
}
