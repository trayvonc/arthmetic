#include"stdafx.h"
#include<iostream>

using namespace std;
#define S 20
#define ERROR 0
#define OK 1
typedef int ElemType;
typedef int Status;
typedef struct
{
	int n;
	int e;
	ElemType **a;
	ElemType noEdge;
}mGraph;

Status Init(mGraph *mg, int nSize, ElemType noEdgeValue,int k)
{
	int i, j;
	mg->n = nSize;
	mg->e = 0;
	mg->noEdge = noEdgeValue;
	mg->a = (ElemType**)malloc(nSize * sizeof(ElemType*));
	if (!mg->a)
		return ERROR;
	for (i = 0; i<mg->n; i++)
	{
		mg->a[i] = (ElemType*)malloc(nSize * sizeof(ElemType));
		for (j = 0; j<mg->n; j++)mg->a[i][j] = mg->noEdge;
		mg->a[i][i] = k;
	}
	return OK;
}
void Output(mGraph *g, int A[])
{
	cout << "COV(A):" << endl;
	for (int i = 0; i<g->n; i++)
	{
		for (int j = 0; j<g->n; j++)
		{
			if (g->a[i][j] == 1)
			{
				cout << "<" << A[i] << "," << A[j] << ">" << "   ";
			}
		}
	}
	cout << endl;
}
int shangquejie(int *A, mGraph *g, int p, int q, int max)
{
	int temp1[S] = { 0 }, temp2[S] = { 0 };
	int tempa = 0, tempb = 0;
	int j1 = 0, j2 = 0;
	for (int i = 0; i < g->n; i++)
	{

		if (g->a[p][i] == 1)
		{
			tempa++;
			temp1[j1] = A[i];
			j1++;
		}
		if (g->a[q][i] == 1)
		{
			tempb++;
			temp2[j2] = A[i];
			j2++;
		}

	}

	for (int m = 0; m < j1; m++)
	{
		for (int n = 0; n < j2; n++)
		{
			if (temp1[m] == temp2[n])
				if (temp1[m] == max)
					return 1;
		}
	}
	return 0;
}
int xiaquejie(int *A, mGraph *g, int p, int q, int min)
{
	int temp3[S] = { 0 }, temp4[S] = { 0 };
	int tempc = 0, tempd = 0;
	int j1 = 0, j2 = 0;
	for (int i = 0; i < g->n; i++)
	{

		if (g->a[i][p] == 1)
		{
			tempc++;
			temp3[j1] = A[i];
			j1++;
		}
		if (g->a[i][q] == 1)
		{
			tempd++;
			temp4[j2] = A[i];
			j2++;
		}

	}

	for (int m = 0; m < j1; m++)
	{
		for (int n = 0; n < j2; n++)
		{
			if (temp3[m] == temp4[n])
				if (temp3[m] == min)
					return 1;
		}
	}
	return 0;
}
//判断盖住关系
mGraph* cover(int A[], int n)
{
	mGraph *g = (mGraph*)malloc(sizeof(mGraph));
	Init(g, n, 0,0);
	for (int i = 0; i<g->n; i++)
	{
		for (int j = 0; j<g->n; j++)
		{
			if (A[i] != A[j] && (A[j] % A[i]) == 0)
			{
				g->a[i][j] = 1;
			}
		}
	}
	for (int i = 0; i<g->n; i++)
	{
		for (int j = 0; j<g->n; j++)
		{
			if (g->a[i][j] == 1)
			{
				for (int u = 0; u<g->n; u++)
				{
					if (g->a[j][u] == 1)
					{
						g->a[i][u] = 0;
					}
				}
			}
		}
	}
	return g;
}

int JudgeGe(int *A, mGraph *g,int n)
{
	int flag1[S][S] = { 0 }, flag2[S][S] = {0};
	int temp1[S] = { 0 }, temp2[S] = { 0 }, temp3[S] = { 0 }, temp4[S] = { 0 };
	int tempa=0, tempb=0,tempc=0,tempd=0;
//	mGraph *g = (mGraph*)malloc(sizeof(mGraph));
	Init(g, n, 0, 1);
	for (int i = 0; i<g->n; i++)
	{
		for (int j = 0; j<g->n; j++)
		{
			if (A[i] != A[j] && (A[j] % A[i]) == 0)
			{
				g->a[i][j] = 1;
			}
		}
	}
	for (int p=0,q=0;A[p]!=A[n-1];q=p)
	{
       do{
		   q++;
		   /*
		   //如果直连则满足条件
		   if (g->a[p][q] == 1 || g->a[q][p] == 1)
		   {
			   flag1[p][q] = 1;
			   flag2[p][q] = 1;
			   
		   }
		   
		   //上确界
		   for (int i = 0; i < g->n; i++)
		   {
			   if (g->a[p][i] = 1)
			   {
				   if (g->a[q][i] = 1)
					   flag1[p][q] = 1;
			   }
		   }
		   //下确界
		   for (int i = 0; i < g->n; i++)
		   {
			   if (g->a[i][p] = 1)
			   {
				   if (g->a[i][q] = 1)
					   flag2[p][q] = 1;
			   }
		   }
		   */
		   //上界
		   int j1 = 0,j2=0;
		   for (int i = 0; i < g->n; i++)
		   {
		
			   if (g->a[p][i] == 1)
			   {
				   tempa++;
				   temp1[j1] = A[i];
				   j1++;
			   }
			   if (g->a[q][i] == 1)
			   {
				   tempb++;
				   temp2[j2] = A[i];
				   j2++;
			   }
			   
		   }
		   
		   for (int m = 0; m < j1; m++)
		   {
			   for(int n = 0; n < j2;n++ )
			   {
				   if (temp1[m] == temp2[n])
					   flag1[p][q] = 1;
			   }
		   }

		   //下界
		   int j3 = 0, j4 = 0;
		   for (int i = 0; i < g->n; i++)
		   {

			   if (g->a[i][p] == 1)
			   {
				   tempc++;
				   temp3[j3] = A[i];
				   j3++;
			   }
			   if (g->a[i][q] == 1)
			   {
				   tempd++;
				   temp4[j4] = A[i];
				   j4++;
			   }

		   }

		   for (int m = 0; m < j3; m++)
		   {
			   for (int n = 0; n < j4; n++)
			   {
				   if (temp1[m] == temp2[n])
					   flag2[p][q] = 1;
			   }
		   }

		     
				
	   } while (A[q] != A[n - 1]);
			p++;
	}
	
	for (int i = 0; i<g->n; i++)
	{
		for (int j = i + 1; j<g->n; j++)
		{
			if (flag1[i][j] == 0 || flag1[i][j] == 0)
			{
				return 0;
			}
		}
	}
	


	return 1;
}
int youbuge(int *A, mGraph *g,int n)
{
	Init(g, n, 0, 1);
	for (int i = 0; i<g->n; i++)
	{
		for (int j = 0; j<g->n; j++)
		{
			if (A[i] != A[j] && (A[j] % A[i]) == 0)
			{
				g->a[i][j] = 1;
			}
		}
	}
	int max=A[0],min=A[0];
	for (int i = 0; i < g->n; i++)
	{
		if (A[i] > max)
			max = A[i];
		if (A[i] < min)
			min = A[i];
	}
	for (int i = 0; i < g->n; i++)
	{
		if (!g->a[i][max] || !g->a[min][i])
			return 0;
	}
	for (int i = 0; i < g->n; i++)
	{
		for (int j = 0; j < g->n; j++)
		{
			if (shangquejie(A, g, i, j, max) && xiaquejie(A, g, i, j, min))
				return 1;
		}
	}
	return 0;
}


int main()
{
	int flag = 0;
	mGraph *g = (mGraph*)malloc(sizeof(mGraph));
	int number, A[S];
	cout << "请输入集合A元素的个数" << endl;
	cin >> number;
	cout << "请输入各元素" << endl;
	for (int i = 0; i<number; i++)
		cin >> A[i];
	g=cover(A, number);
	Output(g, A);
	flag=JudgeGe(A, g, number);
	if (flag)
	{
		cout << "该偏序集为格" << endl;
		//判断有补格
		int i = youbuge(A, g,number);
		if(i)
			cout << "该偏序集为有补格" << endl;
		else
			cout << "该偏序集不是有补格" << endl;
	}
	else
	{
		cout << "该偏序集不是格" << endl;
	}
	system("pause");
	return 0;
}

