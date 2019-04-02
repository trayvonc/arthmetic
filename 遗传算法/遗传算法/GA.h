//参数：
//种群规模N
//交配概率Pc
//变异概率Pm
#include"stdafx.h"
#include"define.h"
#include<random>
using namespace std;
class GA
{
public:
	GA(int N, double pc, double pm);  //根据参数确定种群规模
	GA(int N);
	~GA(void);
	//初始化
	bool initial_pool();
	//开始迭代
	int run(int Maxnum, double d);
	double get_best_value();
	chromosome get_best_chro();
private:
	int c_n;                        //种群规模
	double c_Pc;                    //交配概率Pc
	double c_Pm;                     //变异概率Pm
	gene c_gene_tmp;                //临时基因
	chromosome c_chromosome_temp;   //一条临时染色体
	pool c_pool;                    //种群
									//适应值评价，最优染色体
	bool eval_pool();
	//选择
	bool selection();
	//交配
	bool crossover();
	//变异
	bool mutation();
};

//设定取值范围
inline double get_5_5_rand()
{
	random_device rd;
	double a = double((rd() * 6553) % 100000) / 100000;
	a = (a - .5) * 10;
	return a;
}

//得到0--3之间的随机数
inline int get_0_3_rand()
{
	random_device rd;
	int a = rd() % 4;
	return a;
}

//0--1
inline double get_0_1_rand()
{
	random_device rd;
	double a = double((rd() * 6553) % 100000) / 100000;
	return a;
}

inline double get_chro_value(chromosome chro)
{
	double a = 0, tmp = 0;
	for (int i = 0; i < chro.s_chromosome.size(); i++)
	{
		tmp = tmp + chro.s_chromosome[i] * chro.s_chromosome[i];
	}
	a = 1 / (tmp + 1);
	return a;
}
