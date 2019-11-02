#include "incubator.h"
#include "undirected_graph.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

incubator::incubator(const string &dir)
{
	mdir = dir;
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
		load(file, vc);
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

	merge();
	printf("merge final with %lu (%lu/%lu) combined-graphs for files [%d, %d]\n", gset.size(), vc1.size(), vc2.size(), low, high - 1);

	if(mdir != "" && high - low >= 10)
	{
		char file[10240];
		sprintf(file, "%s/graph-%d-%d.gr", mdir.c_str(), low, high - 1);
		write(file, true);
	}

	assert(merged.size() == gset.size());

	for(int k = 0; k < gset.size(); k++)
	{
		if(merged[k] == true) continue;
		vc.push_back(gset[k]);
	}
	return 0;
}

int incubator::merge()
{
	undirected_graph gr;
	for(int k = 0; k < gset.size(); k++) gr.add_vertex();

	/*
	build_interval_map();
	for(ISMI it = ism.begin(); it != ism.end(); it++)
	{
		set<int> s = it->second;
		*/

	build_splice_map();

	for(MISI::iterator mi = mis.begin(); mi != mis.end(); mi++)
	{
		set<int> &ss = mi->second;
		vector<int> v(ss.begin(), ss.end());
		for(int xi = 0; xi < v.size(); xi++)
		{
			int i = v[xi];
			for(int xj = xi + 1; xj < v.size(); xj++)
			{
				int j = v[xj];
				if(gset[i].chrm != gset[j].chrm) continue;
				if(gset[i].strand != gset[j].strand) continue;

				int c = gset[j].get_overlapped_splice_positions(gset[i].spos);
				// TODO parameter
				double r1 = c * 1.0 / gset[i].spos.size();
				double r2 = c * 1.0 / gset[j].spos.size();
				if(r1 < 0.4 || r2 < 0.4) continue;
				gr.add_edge(i, j);
			}
		}
	}

	//}

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

int incubator::write(const string &file, bool headers)
{
	ofstream fout(file.c_str());
	if(fout.fail()) exit(1);

	for(int k = 0; k < gset.size(); k++)
	{
		if(merged[k] == true) continue;
		//gset[k].build_combined_splice_graph();
		gset[k].write(fout, k, headers);
	}
	fout.close();
	return 0;
}

int incubator::build_splice_map()
{
	mis.clear();
	for(int k = 0; k < gset.size(); k++)
	{
		for(int i = 0; i < gset[k].spos.size(); i++)
		{
			int32_t p = gset[k].spos[i];
			MISI::iterator it = mis.find(p);
			if(it == mis.end())
			{
				set<int> s;
				s.insert(k);
				mis.insert(PISI(p, s));
			}
			else
			{
				it->second.insert(k);
			}
		}
	}
	return 0;
}

int incubator::build_interval_map()
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

int incubator::analyze(const string &file)
{
	vector<combined_graph> vc;
	load(file, vc);
	for(int k = 0; k < vc.size(); k++)
	{
		vc[k].analyze(k);
	}
	return 0;
}

int load(const string &file, vector<splice_graph> &v)
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
	char strand[1024];
	int nodes;

	while(fin.getline(line, 10240, '\n'))
	{
		if(line[0] != '#') continue;
		stringstream sstr(line);
		sstr >> tmp >> gid >> chrm >> nodes >> strand;

		splice_graph gr;
		gr.build(fin, gid, chrm);
		v.push_back(gr);
	}
	
	fin.close();
	return 0;
}

int load(const string &file, vector<combined_graph> &vc)
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
	char strand[1024];
	int nodes;

	while(fin.getline(line, 10240, '\n'))
	{
		if(line[0] != '#') continue;
		stringstream sstr(line);
		sstr >> tmp >> gid >> chrm >> nodes >> strand;

		combined_graph gr;
		gr.build(fin, chrm, strand[0]);
		vc.push_back(gr);
	}
	
	fin.close();
	return 0;
}
