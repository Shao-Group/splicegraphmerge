#include "incubator.h"
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
	return 0;
}

int incubator::write(const string &file)
{
	ofstream fout(file.c_str());
	if(fout.fail()) exit(1);

	for(int k = 0; k < gset.size(); k++)
	{
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
