//������
//��Ⱥ��ģN
//�������Pc
//�������Pm
#include"stdafx.h"
#include"define.h"
#include<random>
using namespace std;
class GA
{
public:
	GA(int N, double pc, double pm);  //���ݲ���ȷ����Ⱥ��ģ
	GA(int N);
	~GA(void);
	//��ʼ��
	bool initial_pool();
	//��ʼ����
	int run(int Maxnum, double d);
	double get_best_value();
	chromosome get_best_chro();
private:
	int c_n;                        //��Ⱥ��ģ
	double c_Pc;                    //�������Pc
	double c_Pm;                     //�������Pm
	gene c_gene_tmp;                //��ʱ����
	chromosome c_chromosome_temp;   //һ����ʱȾɫ��
	pool c_pool;                    //��Ⱥ
									//��Ӧֵ���ۣ�����Ⱦɫ��
	bool eval_pool();
	//ѡ��
	bool selection();
	//����
	bool crossover();
	//����
	bool mutation();
};

//�趨ȡֵ��Χ
inline double get_5_5_rand()
{
	random_device rd;
	double a = double((rd() * 6553) % 100000) / 100000;
	a = (a - .5) * 10;
	return a;
}

//�õ�0--3֮��������
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
