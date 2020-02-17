#include "incubator.h"
#include "undirected_graph.h"
#include "boost/pending/disjoint_sets.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

incubator::incubator(const string &dir)
{
	mdir = dir;
	max_combined_num = 5;
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
	vector<combined_graph> fd;
	binary_merge(files, 0, files.size(), vc);

	grset = vc;
	merged.assign(grset.size(), false);

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

	grset = vc1;
	grset.insert(grset.end(), vc2.begin(), vc2.end());

	merge();

	int n = 0;
	for(int k = 0; k < merged.size(); k++)
	{
		if(merged[k] == false) n++;
	}
	printf("merge final with %d <- %lu + %lu (fixed = %lu) combined-graphs for files [%d, %d]\n", n, vc1.size(), vc2.size(), fixed.size(), low, high - 1);

	if(mdir != "" && high - low >= 10)
	{
		char file[10240];
		sprintf(file, "%s/graph-%d-%d.gr", mdir.c_str(), low, high - 1);
		write(file, true);
	}

	assert(merged.size() == grset.size());

	for(int k = 0; k < grset.size(); k++)
	{
		if(merged[k] == true) continue;
		vc.push_back(grset[k]);
	}
	return 0;
}

int incubator::merge()
{
	// maintain num_combined
	vector<int> csize(grset.size(), 0);
	for(int i = 0; i < grset.size(); i++)
	{
		csize[i] = grset[i].num_combined;
	}

	// disjoint map, maintain clusters
	vector<int> rank(grset.size(), -1);
	vector<int> parent(grset.size(), -1);

	disjoint_sets<int*, int*> ds(&rank[0], &parent[0]);
	for(int k = 0; k < grset.size(); k++)
	{
		ds.make_set(k);
	}

	// maintain the similarity between any two graphs
	vector< map<int, double> > gmap(grset.size());

	vector<PID> vpid;

	MISI mis;
	build_splice_map(mis);
	for(MISI::iterator mi = mis.begin(); mi != mis.end(); mi++)
	{
		set<int> &ss = mi->second;
		vector<int> v(ss.begin(), ss.end());
		for(int xi = 0; xi < v.size(); xi++)
		{
			int i = v[xi];
			assert(grset[i].num_combined < max_combined_num);
			for(int xj = 0; xj < v.size(); xj++)
			{
				int j = v[xj];
				if(i >= j) continue;

				assert(grset[j].num_combined < max_combined_num);

				if(grset[i].chrm != grset[j].chrm) continue;
				if(grset[i].strand != grset[j].strand) continue;

				if(gmap[i].find(j) != gmap[i].end()) continue;

				int c = grset[j].get_overlapped_splice_positions(grset[i].splices);

				double r1 = c * 1.0 / grset[i].splices.size();
				double r2 = c * 1.0 / grset[j].splices.size();
				double r = r1 < r2 ? r1 : r2;

				// TODO parameter
				if(r1 < 0.1 || r2 < 0.1) continue;
				//printf("r1 = %.3lf, r2 = %.3lf, r = %.3lf, size1 = %lu, size2 = %lu\n", r1, r2, r, grset[i].splices.size(), grset[j].splices.size());

				gmap[i].insert(pair<int, double>(j, r));
				vpid.push_back(PID(PI(i, j), r));
			}
		}
	}

	merged.resize(grset.size());
	merged.assign(grset.size(), false);
	sort(vpid.begin(), vpid.end(), compare_graph_overlap);

	for(int i = 0; i < vpid.size(); i++)
	{
		int x = vpid[i].first.first;
		int y = vpid[i].first.second;
		double r = vpid[i].second;
		assert(x < y);

		int px = ds.find_set(x);
		int py = ds.find_set(y);

		if(px == py) continue;
		if(csize[px] >= max_combined_num) continue;
		if(csize[py] >= max_combined_num) continue;

		int sum = csize[px] + csize[py]; 

		printf("combine graph %d (#splices = %lu) and %d (#splices = %lu) with score = %.3lf: %d + %d -> %d\n", 
				x, grset[x].splices.size(), y, grset[y].splices.size(), r, csize[px], csize[py], sum);

		ds.link(px, py);
		int q = ds.find_set(px);
		assert(q == ds.find_set(py));

		if(q == x) 
		{
			grset[x].combine(grset[y]);
			merged[y] = true;
		}
		else if(q == y) 
		{
			grset[y].combine(grset[x]);
			merged[x] = true;
		}
		else assert(false);

		csize[q] = sum;
		assert(sum == grset[q].num_combined);
	}

	vector<combined_graph> cc;
	vector<bool> bb;
	for(int i = 0; i < grset.size(); i++)
	{
		if(grset[i].num_combined >= max_combined_num)
		{
			fixed.push_back(grset[i]);
		}
		else
		{
			cc.push_back(grset[i]);
			bb.push_back(merged[i]);
		}
	}

	grset = cc;
	merged = bb;
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

		//printf("final combine %d and %d with %lu and %lu vertices\n", x, k, grset[x].imap.size(), grset[k].imap.size());

		grset[x].combine(grset[k]);
		merged[k] = true;
	}
	return 0;
}

int incubator::write(const string &file, bool headers)
{
	ofstream fout(file.c_str());
	if(fout.fail()) exit(1);

	for(int k = 0; k < grset.size(); k++)
	{
		if(merged[k] == true) continue;
		grset[k].write(fout, k, headers);
	}

	for(int k = 0; k < fixed.size(); k++)
	{
		fixed[k].write(fout, k, headers);
	}

	fout.close();
	return 0;
}

int incubator::build_splice_map(MISI &mis)
{
	mis.clear();
	for(int k = 0; k < grset.size(); k++)
	{
		for(int i = 0; i < grset[k].splices.size(); i++)
		{
			int32_t p = grset[k].splices[i];
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

int incubator::print()
{
	for(int k = 0; k < grset.size(); k++) grset[k].print(k);
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
		sstr >> tmp >> gid >> chrm >> strand;

		combined_graph gr(line);
		gr.build(fin, chrm, strand[0]);
		vc.push_back(gr);
	}
	fin.close();
	return 0;
}

bool compare_graph_overlap(const PID &x, const PID &y)
{
	return x.second > y.second;
}
