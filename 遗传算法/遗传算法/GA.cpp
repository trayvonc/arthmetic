#include"stdafx.h"
#include "GA.h"

#include"define.h"
#include <stdio.h>
#include <tchar.h>
#include <iostream>
using namespace std;

inline bool sort_value(const chromosome &aim1, const chromosome &aim2);

GA::GA(int N, double pc, double pm)
{
	c_n = N;
	c_Pc = pc;
	c_Pm = pm;
}
GA::GA(int N)
{
	c_n = N;
	c_Pc = 0.5;   //һ��ȡֵ0.4--0.99
	c_Pm = 0.1;   //һ��ȡֵ0.001-0.1
}

GA::~GA(void)
{
	cout << "����GA���ڴ�------------------" << endl;
}

bool GA::initial_pool()
{
	c_pool.clear();
	//��ʼ��Ⱦɫ��Ⱥ
	for (int j = 0; j < c_n; j++)
	{
		c_chromosome_temp.s_chromosome.clear();
		//�õ�һ��Ⱦɫ��
		for (int i = 0; i<4; i++)
		{
			c_chromosome_temp.s_chromosome.push_back(
				get_5_5_rand());
		}
		c_chromosome_temp.s_value = 0;
		c_pool.push_back(c_chromosome_temp);
	}
	cout << "��ʼֵ��" << endl;
	cout << "x1:" << get_best_chro().s_chromosome[0]
		<< "  x2:" << get_best_chro().s_chromosome[1]
		<< "  x3:" << get_best_chro().s_chromosome[2]
		<< "  x4:" << get_best_chro().s_chromosome[3] << endl;
	return true;
}

//��Ӧֵ����
bool GA::eval_pool()
{
	int i = 0;
	while (i != c_pool.size())
	{
		double temp = 0;
		for (int j = 0; j<c_pool[i].s_chromosome.size(); j++)
		{
			temp = temp +
				c_pool[i].s_chromosome[j] * c_pool[i].s_chromosome[j];
		}
		c_pool[i].s_value = 1 / (temp + 1);
		i++;
	}
	sort(c_pool.begin(), c_pool.end(), sort_value);

	return true;
}

inline bool sort_value(const chromosome &aim1, const chromosome &aim2)
{
	return (aim1.s_value>aim2.s_value);
}


//���ݶԴ���Ȼ����Ӧ������������Ȼѡ����ʤ��̭
bool GA::selection()
{
	double all_value = 0;
	int i = 0;
	while (true)
	{
		all_value += c_pool[i].s_value;
		i++;
		if (i >= c_pool.size())
		{
			i = 0;
			break;
		}
	}
	//������Ӧֵfit
	while (true)
	{
		c_pool[i].s_fit = c_pool[i].s_value / all_value;
		i++;
		if (i >= c_pool.size())
		{
			i = 0;
			break;
		}
	}
	//����ѡ��
	pool tmp_pool;
	for (int k = 0; k<c_pool.size(); k++)
	{
		double m = 0, r = get_0_1_rand();
		for (int j = 0; j<c_pool.size(); j++)
		{
			m = m + c_pool[j].s_fit;
			if (r <= m)
			{
				tmp_pool.push_back(c_pool[j]);
				break;
			}
		}
	}
	c_pool.clear();
	c_pool = tmp_pool;
	return true;
}

//Ⱦɫ�彻��
bool GA::crossover()
{
	vector<chromosome>::iterator it = c_pool.begin();
	pool cross_pool, no_cross_pool;
	//�õ�����غͷǽ����
	for (int i = 0; i<c_pool.size(); i++)
	{
		double a = get_0_1_rand();
		if (a<c_Pc)
			cross_pool.push_back(c_pool[i]);
		else
			no_cross_pool.push_back(c_pool[i]);
	}
	//���н���
	//����Ҫ����
	if (cross_pool.size() == 0)
	{
		c_pool.clear();
		c_pool = no_cross_pool;
		return true;
	}
	else if (cross_pool.size() == 1)
	{
		c_pool.clear();
		c_pool = no_cross_pool;
		c_pool.push_back(cross_pool[0]);
		return true;
	}
	else
	{
		chromosome cross_tmp1, cross_tmp2, cross_1, cross_2;
		for (int i = 0; i<4; i++)
		{
			cross_tmp1.s_chromosome.push_back(0);
			cross_tmp2.s_chromosome.push_back(0);
		}
		//��������λż��
		if (cross_pool.size() % 2 == 0)
		{
			for (int i = 0; i<cross_pool.size(); i += 2)
			{
				int lc = get_0_3_rand();
				//��ʼ����
				for (int j = 0; j<4; j++)
				{
					cross_1 = cross_pool[i];
					cross_2 = cross_pool[i + 1];
					if (lc<j)
					{
						cross_tmp1.s_chromosome[j] = cross_pool[i].s_chromosome[j];
						cross_tmp2.s_chromosome[j] = cross_pool[i + 1].s_chromosome[j];
					}
					else
					{
						cross_tmp2.s_chromosome[j] = cross_pool[i].s_chromosome[j];
						cross_tmp1.s_chromosome[j] = cross_pool[i + 1].s_chromosome[j];
					}

				}
				//������������µ�Ⱦɫ�����nocrosspool
				if (get_chro_value(cross_tmp1)>get_chro_value(cross_1))
					no_cross_pool.push_back(cross_tmp1);
				else
					no_cross_pool.push_back(cross_1);
				if (get_chro_value(cross_tmp2)>get_chro_value(cross_2))
					no_cross_pool.push_back(cross_tmp2);
				else
					no_cross_pool.push_back(cross_2);
			}
		}
		//��������λ���������������һ��
		else
		{
			no_cross_pool.push_back(cross_pool[0]);
			for (int i = 1; i<cross_pool.size(); i += 2)
			{
				int lc = get_0_3_rand();
				//��ʼ����
				for (int j = 0; j<4; j++)
				{
					cross_1 = cross_pool[i];
					cross_2 = cross_pool[i + 1];
					if (lc<j)
					{
						cross_tmp1.s_chromosome[j] = cross_pool[i].s_chromosome[j];
						cross_tmp2.s_chromosome[j] = cross_pool[i + 1].s_chromosome[j];
					}
					else
					{
						cross_tmp2.s_chromosome[j] = cross_pool[i].s_chromosome[j];
						cross_tmp1.s_chromosome[j] = cross_pool[i + 1].s_chromosome[j];
					}

				}
				//������������µ�Ⱦɫ�����nocrosspool
				//�������������һ�����ã��ͰѸ����滻�Ӵ�
				if (get_chro_value(cross_tmp1)>get_chro_value(cross_1))
					no_cross_pool.push_back(cross_tmp1);
				else
					no_cross_pool.push_back(cross_1);
				if (get_chro_value(cross_tmp2)>get_chro_value(cross_2))
					no_cross_pool.push_back(cross_tmp2);
				else
					no_cross_pool.push_back(cross_2);
			}
		}
		c_pool.clear();
		c_pool = no_cross_pool;
		return true;
	}
}

//���ݱ����ʣ�����Ⱦɫ�����
bool GA::mutation()
{
	int i = 0;
	while (true)
	{
		double rm = get_0_1_rand();
		chromosome bianyi = c_pool[i];
		if (rm<c_Pm)
		{
			int lc = get_0_3_rand();//����λ��
			bianyi.s_chromosome[lc] = get_5_5_rand();
			if (get_chro_value(bianyi)>get_chro_value(c_pool[i]))
				c_pool[i] = bianyi;
		}
		i++;
		if (i >= c_pool.size())
		{
			i = 0;
			break;
		}
	}
	return true;
}


int GA::run(int Maxnum, double d)
{
	int i = 0;
	double last_value = 0, current_value;
	while (true)
	{
		i++;
		eval_pool();
		current_value = c_pool[0].s_value;
		selection();
		crossover();
		mutation();
		if (i % 40 == 0)
			cout << "current_value:" << current_value << endl;
		if (abs(current_value - last_value)<d)
		{
			//cout<<"���ڵ������ﵽҪ��"<<endl;
			//break;
		}
		if (i >= Maxnum)
		{
			cout << "���������������㷨������" << endl;
			break;
		}
		last_value = current_value;
	}
	eval_pool();
	return i;
}

double GA::get_best_value()
{
	return c_pool[0].s_value;
}

chromosome GA::get_best_chro()
{
	return c_pool[0];
}