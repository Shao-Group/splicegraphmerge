#include "phasing_path.h"
#include <cassert>

int phasing_path::combine(const phasing_path &p)
{
	for(int k = 0; k < p.paths.size(); k++)
	{
		combine(p.paths[k], p.weights[k], p.counts[k]);
	}
	return 0;
}

int phasing_path::combine(const vector<int32_t> &v, double w, int c)
{
	assert(v.size() >= 2);

	map< int32_t, set<PI> >::iterator x1 = index.find(v.front()); 
	if(x1 == index.end())
	{
		add_new_path(v, w, c);
		return 0;
	}

	map< int32_t, set<PI> >::iterator x2 = index.find(v.back()); 
	if(x2 == index.end()) 
	{
		add_new_path(v, w, c);
		return 0;
	}

	set<PI> &s1 = x1->second;
	set<PI> &s2 = x2->second;

	vector<int> vv;
	vector<PI> pp;

	binary_merge(s1, s2, vv, pp);


	int f = -1;
	for(int i = 0; i < vv.size(); i++)
	{
		int p = vv[i];
		int s = pp[i].first;
		int t = pp[i].second;
		assert(s < t);
		assert(paths[p][s] == v.front());
		assert(paths[p][t] == v.back());
		if(t - s != v.size() - 1) continue;
		if(binary_verify(v, paths[p], s, 0, v.size()) == false) continue;
		f = p;
		break;
	}

	if(f == -1)
	{
		add_new_path(v, w, c);
		return 0;
	}
	else if(paths[f].size() == v.size())
	{
		weights[f] += w;
		counts[f] += c;
	}
	return 0;
}

int phasing_path::add_new_path(const vector<int32_t> &v, double w, int c)
{
	int n = paths.size();
	paths.push_back(v);
	weights.push_back(w);
	counts.push_back(c);

	for(int i = 0; i < v.size(); i++)
	{
		int32_t x = v[i];
		PI p(n, i);
		if(index.find(x) == index.end())
		{
			set<PI> s;
			s.insert(p);
			index.insert(pair<int32_t, set<PI> >(x, s));
		}
		else 
		{
			index[x].insert(p);
		}
	}
	return 0;
}

bool phasing_path::binary_verify(const vector<int32_t> &v, const vector<int32_t> &p, int o, int s, int t)
{
	if(s + 1 >= t) return true;
	int m = (s + t) / 2;
	if(v[m] != p[m + o]) return false;
	if(binary_verify(v, p, o, s, m) == false) return false;
	if(binary_verify(v, p, o, m, t) == false) return false;
	return true;
}

int phasing_path::binary_merge(const set<PI> &s1, const set<PI> &s2, vector<int> &v, vector<PI> &p)
{
	v.clear();
	p.clear();
	set<PI>::const_iterator x1 = s1.begin();
	set<PI>::const_iterator x2 = s2.begin();
	while(x1 != s1.end() && x2 != s2.end())
	{
		if(x1->first == x2->first)
		{
			v.push_back(x1->first);
			p.push_back(PI(x1->second, x2->second));
			x1++;
			x2++;
		}
		else if(x1->first < x2->first)
		{
			x1++;
		}
		else
		{
			x2++;
		}
	}
	return 0;
}
