// 图的生成及欧拉回路的确定.cpp: 定义控制台应用程序的入口点。
//
/* 实验要求随机生成n个顶点的无向图。n由输入得到，声明一个整型的二维数组来表示无向图的关系矩阵，随机可使用c语言的srand()和rand()函数，
   对二维数组的值随机生成0或1，0表示两点间无直接路径，1表示两点间有直接的路径。

      判断一个无向图是否为欧拉图或者半欧拉图用的是书上定理7-4.1：
	  无向图G具有一条欧拉路，当且仅当G是连通的，且有零个或两个奇数度结点。所以问题主要分为三个部分。

1.判断无向图是否为连通图，这里采用深度优先搜索，在随机生成图的时候找一条边，任意取其一个端点做起点，对图进行深度优先遍历，
如果为连通图则从其中任何一个起点开始都能遍历整张图，反之，若从一点起遍历完还有未访问的点，则为非连通图。

2.判断奇数度结点的个数，这个很容易，直接计算每个点的度然后记录一下奇数度结点的个数就可以了。

3.若图是一个欧拉图或半欧拉图，输出其欧拉回路或欧拉路，这里采用弗罗莱算法：设G为一个无向欧拉图。

   （1）任取G中一顶点v0,令P0=v0

   （2）假设沿Pi=v0e1v1e2…eivi走到顶点vi,按下面方法从E(G)-{e1,e2,…ei}中选ei+1: 1) ei+1与vi关联 2)除非无别的边可供选择，
      否则ei+1不应该是GI=G-{e1,e2,…ei}中的桥。（3）当（2）不能再进行时算法停止。所得到的简单回路Pm=v0e1v1e2v2…emvm为G中的一条欧拉路

*/

#include "stdafx.h"
#include <iostream>  
#include <ctime>  
#include <cstdlib>  
using namespace std;
int map[55][55], n, m = 0; //图，节点数，边数，  
int ans[50], cnt = 0; //记录欧拉路的路径，路径数  
bool vis[50]; //标记点是否被访问  
int st; //判是否为连通图，搜索的起点  

struct stack
{
	int top, node[100];
}s; //定点的栈结构  

void show() //输出随机生成的无向图关系矩阵  
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			cout << map[i][j] << " ";
		cout << endl;
	}
}

void init() //初始化函数  
{
	bool flag = false;  //标记是否生成图  
	memset(map, 0, sizeof(map));
	memset(vis, false, sizeof(vis));
	memset(ans, 0, sizeof(ans));
	srand((unsigned)time(NULL));
	for (int i = 0; i < n; i++)
		for (int j = i + 1; j < n; j++)
			map[i][j] = 0 + rand() % 2; //随机生成无向图  
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			if (map[i][j])
			{
				flag = true;    //最少有一条边  
				m++;            //记录边的条数  
				map[j][i] = 1;  //无向图的定义  
				st = i;         //有一条边则其任意一个端点设为起点  
			}
	while (!flag)
		init();
}

void DFS2(int x)  //深度优先搜索  
{
	vis[x] = true;
	for (int i = 0; i < n; i++)
		if (!vis[i] && map[x][i])
			DFS2(i);
}

//从所设定的起点深度优先遍历图，若有一个点没被访问，则为非连通图  
bool judge()
{
	DFS2(st);
	for (int i = 0; i < n; i++)
		if (!vis[i])
			return false;
	return true;
}

void DFS(int x) //深度优先搜索  
{
	s.top++;
	s.node[s.top] = x;
	for (int i = 0; i < n; i++)
	{
		if (map[i][x] > 0)
		{
			map[i][x] = 0; //删边操作  
			map[x][i] = 0;
			DFS(i);
			break;
		}
	}
}

void Fleury(int x)  //Fleury算法  
{
	int b;
	s.top = 0;
	s.node[s.top] = x; //起点入栈  
	while (s.top >= 0)
	{
		b = 0;
		for (int i = 0; i < n; i++)
		{
			if (map[s.node[s.top]][i] > 0)
			{
				b = 1;
				break;
			}
		}
		if (b == 0) //如果没有可扩展的点，则记录下该点并将其出栈  
		{
			ans[cnt++] = s.node[s.top] + 1;
			s.top--;
		}
		else //如果有，则将其出栈并继续搜索  
		{
			s.top--;
			DFS(s.node[s.top + 1]);
		}
	}
	cout << endl;
}

void answer()  //输出答案  
{
	for (int i = 0; i < cnt; i++)
		cout << ans[i] << " ";
	cout << endl;
}

int main()
{
	int num = 0, start = 0, degree; // 奇度顶点个数, 欧拉路的起点, 每个顶点的度  
	cout << "请输入n (1 ~ 20):";
	cin >> n;
	init();
	cout << "生成的无向图为: " << endl;
	show();
	if (!judge())
	{
		cout << "非连通图" << endl;
		return 0;
	}
	//如果存在奇度顶点，则从奇度顶点出发，否则从0出发  
	for (int i = 0; i < n; i++)
	{
		degree = 0;
		for (int j = 0; j < n; j++)
			degree += map[i][j];
		if (degree % 2)
		{
			start = i;
			num++;
		}
	}
	//无向图具有一条欧拉路，当且仅当G是连通的，且有0个或2个奇数度结点  
	if (num == 0 || num == 2)
	{
		Fleury(start);
		//欧拉路径的头和尾相等，则说明欧拉路是回路  
		if (ans[0] == ans[cnt - 1])
			cout << "该图为欧拉图，欧拉回路为: ";
		else
			cout << "该图为半欧拉图，欧拉路为: ";
		answer();
	}
	else
		cout << "非欧拉图或半欧拉图" << endl;
	system("pause");
	return 0;
}
