#include"stdafx.h"
#include <vector>
#include <stdio.h>
using namespace std;
//using std::vector;
//vector<int> vInts;
struct chromosome
{
	vector<double> s_chromosome;
	double s_value;
	double s_fit;
};
//基因
#define gene double
//种群
#define pool vector<struct chromosome>
