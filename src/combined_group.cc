#include <algorithm>
#include "combined_group.h"
#include "boost/pending/disjoint_sets.hpp"

combined_group::combined_group(string c, char s)
{
	chrm = c;
	strand = s;
}

int combined_group::add_graph(const combined_graph &gr)
{
	assert(gr.chrm == chrm);
	gset.push_back(gr);
	return 0;
}

int combined_group::resolve()
{
	build_splice_map();
	build_similarity();
	combine_graphs();
	return 0;
}

int combined_group::build_splice_map()
{
	mis.clear();
	for(int k = 0; k < gset.size(); k++)
	{
		for(int i = 0; i < gset[k].splices.size(); i++)
		{
			int32_t p = gset[k].splices[i];
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

int combined_group::build_similarity()
{
	// maintain the similarity between any two graphs
	vector< set<int> > gmap(gset.size());		// success; into vpid
	vector< set<int> > fmap(gset.size());		// compared but failed

	for(MISI::iterator mi = mis.begin(); mi != mis.end(); mi++)
	{
		set<int> &ss = mi->second;
		vector<int> v(ss.begin(), ss.end());
		for(int xi = 0; xi < v.size(); xi++)
		{
			int i = v[xi];
			for(int xj = 0; xj < v.size(); xj++)
			{
				int j = v[xj];
				if(i >= j) continue;

				assert(gset[i].chrm == gset[j].chrm);
				assert(gset[i].strand == gset[j].strand);

				if(gmap[i].find(j) != gmap[i].end()) continue;
				if(fmap[i].find(j) != fmap[i].end()) continue;

				int c = gset[j].get_overlapped_splice_positions(gset[i].splices);
				double r = 2.0 * c / (gset[i].splices.size() + gset[j].splices.size());

				// TODO parameters
				bool b = true;
				if(c <= 1.5) b = false;
				if(r <= 0.8) b = false;

				//printf("r1 = %.3lf, r2 = %.3lf, r = %.3lf, size1 = %lu, size2 = %lu\n", r1, r2, r, gset[i].splices.size(), gset[j].splices.size());

				if(b == true)
				{
					gmap[i].insert(j);
					vpid.push_back(PID(PI(i, j), r));
				}
				else
				{
					fmap[i].insert(j);
				}
			}
		}
	}

	return 0;
}

int combined_group::combine_graphs()
{
	// maintain num_combined
	vector<int> csize(gset.size(), 0);
	for(int i = 0; i < gset.size(); i++)
	{
		csize[i] = gset[i].num_combined;
	}

	// disjoint map, maintain clusters
	vector<int> rank(gset.size(), -1);
	vector<int> parent(gset.size(), -1);

	disjoint_sets<int*, int*> ds(&rank[0], &parent[0]);
	for(int k = 0; k < gset.size(); k++)
	{
		ds.make_set(k);
	}

	sort(vpid.begin(), vpid.end(), compare_graph_similarity);

	for(int i = 0; i < vpid.size(); i++)
	{
		int x = vpid[i].first.first;
		int y = vpid[i].first.second;
		double r = vpid[i].second;
		assert(x < y);

		int px = ds.find_set(x);
		int py = ds.find_set(y);

		if(px == py) continue;

		int sum = csize[px] + csize[py]; 

		printf("combine graph %d (#splices = %lu) and %d (#splices = %lu) with score = %.3lf: %d + %d -> %d\n", 
				x, gset[x].splices.size(), y, gset[y].splices.size(), r, csize[px], csize[py], sum);

		ds.link(px, py);
		int q = ds.find_set(px);
		assert(q == ds.find_set(py));
		assert(q == px || q == py);
		csize[q] = sum;
	}

	vector<combined_graph> cc(gset.size());
	vector<int> pp(gset.size(), -1);
	for(int i = 0; i < gset.size(); i++)
	{
		int p = ds.find_set(i);
		pp[i] = p;
		if(p == i) cc[i] = gset[i];
	}

	for(int i = 0; i < gset.size(); i++)
	{
		if(pp[i] == i) continue;
		int p = pp[i];
		cc[p].combine(gset[i]);
	}

	mset.clear();
	for(int i = 0; i < gset.size(); i++)
	{
		if(i == pp[i])
		{
			mset.push_back(cc[i]);
		}
	}
	return 0;
}

bool compare_graph_similarity(const PID &x, const PID &y)
{
	return x.second > y.second;
}
