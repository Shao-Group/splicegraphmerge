#ifndef __PHASING_PATH_H__
#define __PHASING_PATH_H__

#include <stdint.h>
#include <vector>
#include <map>
#include <set>

using namespace std;

typedef pair<int, int> PI;

class phasing_path
{
public:
	vector< vector<int32_t> > paths;
	vector<double> weights;
	vector<int> counts;
	map< int32_t, set<PI> > index;

public:
	int combine(const phasing_path &p);
	int combine(const vector<int32_t> &v, double w, int c);
	int size();

private:
	int add_new_path(const vector<int32_t> &v, double w, int c);
	int binary_merge(const set<PI> &s1, const set<PI> &s2, vector<int> &v, vector<PI> &p);
	bool binary_verify(const vector<int32_t> &v, const vector<int32_t> &p, int o, int s, int t);
};

#endif
