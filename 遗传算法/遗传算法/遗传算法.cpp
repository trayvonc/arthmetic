// 遗传算法.cpp: 定义控制台应用程序的入口点。
//

#include "stdafx.h"


#include "GA.h"
//#include "test.h"
#include <stdio.h>
#include <tchar.h>
#include <iostream>
using namespace std;
#define ERROR 0.00000001
#define MAXNUM 1000

int main()
{
	GA my_GA(10);
	my_GA.initial_pool();
	cout << "迭代次数：" << my_GA.run(MAXNUM, ERROR) << endl;
	cout << "最优解：" << my_GA.get_best_value() << endl;
	cout << "x1:" << my_GA.get_best_chro().s_chromosome[0]
		<< "  x2:" << my_GA.get_best_chro().s_chromosome[1]
		<< "  x3:" << my_GA.get_best_chro().s_chromosome[2]
		<< "  x4:" << my_GA.get_best_chro().s_chromosome[3] << endl;
	system("pause");
}
