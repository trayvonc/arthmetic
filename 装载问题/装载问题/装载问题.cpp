// 装载问题.cpp: 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
using namespace std;

template <class T>
class Loading
{
private:
	int 	n,	//集装箱数
		*x,		//当前解,前结点的路径
		*bestx;	//当前第一艘船的最优解
	T	c1,	//第一艘轮船的核定载重量
		c2,	//第二艘轮船的核定载重量
		*w,	//集装箱重量数组
		total,	//所有集装箱重量之和
		cw,		//当前第一艘船的载重量
		bestw,	//当前第一艘船的最优载重量
		r;		//剩余集装箱总重量
public:
	Loading(int n,T c1,T c2,T *w)		//构造函数
	{
		this->n = n;
		this->c1 = c1;
		this->c2 = c2;
		this->x = new int[n + 1];
		this->w = w;
		this->bestx = new int[n + 1];
		bestw = 0;
		cw = 0;
		r = 0;
		for (int i = 1; i <= n; i++)  
			r += w[i];
		total = 0;
		for (int i = 1; i <= n; i++)
			total += w[i];
    }

	~Loading()	//析构函数	
	{
		delete x;
		delete bestx;
	}
	void Backtrack(int i);	//找到最接近第一艘轮船载重c1的最佳装载方案，
										//最优载重值bestw，最优解数组bestx。
	void Show(int n, int bestw);	//输出整个装载方案
	int getbestw()
	{
		return bestw;
	}
	};

template <class T>
void Loading<T>::Backtrack(int i)
{	//搜索第i层结点
	if (i>n)
	{//到达叶节点
	 if (cw>bestw)
	 {
		for (int j = 1; j <= n; j++)	bestx[j] = x[j];
			bestw = cw;
	 }
		return;
	}
	//搜索子树
	r -= w[i];
	if (cw + w[i] <= c1)		//x[i]=1时的可行解约束条件
	{//搜索左子树
		x[i] = 1;
		cw += w[i];
		Backtrack(i + 1);
		cw -= w[i];
	}
	if (cw + r>bestw)		//x[i]=0时增加的约束函数，剪去不含最优解的分枝
	{//搜索右子树
		x[i] = 0;
		Backtrack(i + 1);
	}
	r += w[i];
}

template <class T>
void Loading<T>::Show(int n,int bestw)
{
	cout << "集装箱总重为：" << total << endl;
	cout << "所取物品：";
	for (int i = 1; i < n + 1; i++)
		cout << bestx[i];
		cout << endl;

		cout <<"放入第一艘船的集装箱总重为 ："<< bestw<<endl;

		if (total - bestw > c2)
			cout << "不能完成装载！" << endl;
		else
		{
			cout << "放入第二艘船的集装箱总重为 ：" << total-bestw << endl;
			cout << "能够完成装载！" << endl;
		}
				
}

int main()
{
	int w[6], c,n=5;
	cout << "请输入5个集装箱的重量" << endl;
	for (int i = 1; i < 6; i++)
		cin >> w[i];
	cout << "输入第一艘轮船的载重量" <<endl;
	cin >> c;
	cout << endl;
    Loading<int> ld(n,c,40,w);

    ld.Backtrack(1);

    ld.Show(n,ld.getbestw());
    system("pause");
	return 0;
}

