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
//����
#define gene double
//��Ⱥ
#define pool vector<struct chromosome>
