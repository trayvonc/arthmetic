// NQueen.cpp: 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
using namespace std;
int cnt = 0;

bool Place(int k, int i, int *x)
{
	for (int j = 0; j < k; j++)
		if ((x[j] == i) || (abs(x[j] - i) == abs(j - k)))
			return false;
	return true;
}

void NQueens(int k, int n, int *x)
{
	for (int i = 0; i < n; i++)
	{
		if (Place(k, i, x))
		{
			x[k] = i;
			if (k == n - 1)
			{
				for (i = 0; i < n; i++)
					cout << x[i] << " ";
				cout << endl;
				cnt++;
			}
			else
			{
				NQueens(k + 1, n, x);
			}
		}
	}
}
void NQueens(int n, int *x)
{
	NQueens(0, n, x);
}

int main()
{
	int x[8];
	cnt = 0;
	for (int i = 0; i < 8; i++)
		x[i] = -1;
	NQueens(8, x);
	printf("Answers num: %d\n", cnt);
	system("pause");
	return 0;
}

