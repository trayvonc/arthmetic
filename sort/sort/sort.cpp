// sort.cpp: 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
//结构体
#define KeyTye int
#define DataType int
#define N 10000
#define maxSize 100000
typedef struct entry
{
	KeyTye key;
	DataType data;
}Entry;
typedef struct list {
	int n;
	Entry D[maxSize];
}List;
//简单排序
int FindMin(List list, int startIndex)
{
	int i, minIndex = startIndex;
	for (i = startIndex + 1; i<list.n; i++)
	{
		if (list.D[i].key<list.D[minIndex].key)
			minIndex = i;
	}
	return minIndex;
}
void Swap(Entry *D, int i, int j)
{
	if (i == j)return;
	Entry temp = *(D + i);
	*(D + i) = *(D + j);
	*(D + j) = temp;
}
void SelectSort(List *list)
{
	int minIndex, startIndex = 0;
	while (startIndex<list->n - 1)
	{
		minIndex = FindMin(*list, startIndex);
		Swap(list->D, startIndex, minIndex);
		startIndex++;
	}
}
//直接插入排序算法
void InsertSort(List *list)
{
	int i, j;
	for (i = 1; i<list->n; i++)
	{
		Entry insertItem = list->D[i];
		for (j = i - 1; j >= 0; j--)
		{
			if (insertItem.key<list->D[j].key)
				list->D[j + 1] = list->D[j];
			else break;
		}
		list->D[j + 1] = insertItem;
	}
}
//冒泡排序法
void BubbleSORT(List* list)
{
	int i, j;
	bool isSwap = false;
	for (i = list->n - 1; i>0; i--)
	{
		for (j = 0; j<i; j++)
		{
			if (list->D[j].key>list->D[j + 1].key)
			{
				Swap(list->D, j, j + 1);
				isSwap = true;
			}
		}
		if (!isSwap) break;
	}
}
//快速排序算法
int Partition(List *list, int low, int high)
{
	int i = low, j = high+1 ;
	Entry pivot = list->D[low];
	do {
		do i++; while (list->D[i].key<pivot.key);
		do j--; while (list->D[j].key>pivot.key);
		if (i<j) Swap(list->D, i, j);
	} while (i<j);
	Swap(list->D, low, j);
	return j;
}
void QuickSort(List *list, int low, int high)
{
	int k;
	if (low<high)
	{
		k = Partition(list, low, high);
		QuickSort(list, low, k - 1);
		QuickSort(list, k + 1, high);
	}
}
void QuickSort(List *list)
{
	QuickSort(list, 0, list->n - 1);
}
//两路合并排序算法
void Merge(List *list, Entry* temp, int low, int n1, int n2)
{
	int i = low, j = low + n1;
	while (i <= (low + n1 - 1 )&& j <= (low + n1 + n2 - 1))
	{
		if (list->D[i].key <= list->D[j].key)
			*temp++ = list->D[i++];
		else *temp++ = list->D[j++];
	}
	while (i <= (low + n1 - 1))
		*temp++ = list->D[i++];
	while (j <= (low + n1 + n2 - 1))
		*temp++ = list->D[j++];
}
void MergeSort(List *list)
{
	Entry temp[maxSize];
	int low, n1, n2, i, size = 1;
	while (size<list->n) {
		low = 0;
		while (low + size<list->n) {
			n1 = size;
			if (low + size * 2<list->n)n2 = size;
			else n2 = list->n - low - size;
			Merge(list, temp+low, low, n1, n2);
			low += n1 + n2;
		}
		for (i = 0; i<list->n; i++)
			list->D[i] = temp[i];
		size *= 2;
	}
}
void Output(list *list)
{
	for (int i = 0; i<list->n; i++) {
		printf("%d", list->D[i].key);
		printf("  ");
	}
}


int main()
{
	time_t start, finish;
	int k;
	double duration;//持续时间
	list *a=(List*)malloc(sizeof(List)), *b=(List*)malloc(sizeof(List));

	srand((unsigned)time(NULL));
	a->n = N;
	b->n = N;


	for (k = 0; k<N; k++)
	{
		a->D[k].key = rand() % 1001;
		b->D[k].key = a->D[k].key;
	}
	b->D[N].key = 9999;                   //针对解决快排出错
	start = clock();
	SelectSort(b);
	//Output(b); 
	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	printf("简单排序用时%lf秒\n", duration);
	

	for (k = 0; k<N; k++)
	{
		b->D[k].key = a->D[k].key;
	}
	start = clock();
	InsertSort(b);
	//Output(b);
	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	printf("插入排序用时%lf秒\n", duration);
	

	for (k = 0; k<N; k++)
	{
		b->D[k].key = a->D[k].key;
	}
	start = clock();
	BubbleSORT(b);
	//Output(b);
	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	printf("冒泡排序用时%lf秒\n", duration);
	
	
	for (k = 0; k<N; k++)
	{
	b->D[k].key = a->D[k].key;
	}
	start = clock();
	QuickSort(b);
	//Output(b);
	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	printf("快速排序用时%lf秒\n", duration);
	


	for (k = 0; k<N; k++)
	{
		b->D[k].key = a->D[k].key;
	}
	start = clock();
	MergeSort(b);
	//Output(b);
	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	printf("两路合并用时%lf秒\n", duration);
	

	system("pause");
	return 0;
}