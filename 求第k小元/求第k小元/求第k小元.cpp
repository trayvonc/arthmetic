#include<iostream>
#include"stdafx.h"
using namespace std;
enum ResultCode { OutOfBounds, Success };
class SortableList
{
public:
	SortableList(int mSize);
	~SortableList();
	void Input();
	void Output();
	ResultCode Select(int &x, int k);
private:
	int *l;
	int maxSize;
	int n;
	void Swap(int i, int j);
	void InsertSort(int left, int right);
	int Partition(int left, int right);
	int Select(int k, int left, int right, int r);
};
SortableList::SortableList(int mSize)
{
	maxSize = mSize;
	l = new int[maxSize];
	n = 0;
}
SortableList::~SortableList() { delete[]l; }
void SortableList::Input()
{
	cout << " 请输入带排序的数组 \n";
	for (int i = 0; i < maxSize; i++) {
		cin >> l[i]; 
			if (l[i] == -1)
				break;
		n++;
	}
}
void SortableList::Output()
{
	for (int i = 0; i < n; i++)
		cout << l[i] << " ";
}
void SortableList::Swap(int i, int j)
{
	int c = l[i];
	l[i] = l[j];
	l[j] = c;
}
int SortableList::Partition(int left, int right)
{
	int i = left, j = right + 1;
	do {
		do i++; while (l[i] < l[left]);
		do j--; while (l[j] > l[left]);
		if (i < j) Swap(i, j);
	} while (i < j);
	Swap(left, j);
	return j;
}
ResultCode SortableList::Select(int &x, int k)
{
	if (n <= 0 || k > n || k <= 0)return OutOfBounds;
	int j = Select(k, 0, n - 1, 5);
	x = l[j];
	return Success;
}
int SortableList::Select(int k, int left, int right, int r)
{
	int n = right - left + 1;
	if (n <= r) { // 若问题足够小，使用直接插入排序
		InsertSort(left, right); 
			return left + k - 1; // 取其中的第 k小元素，其下标为 left+k-1
	}
	for (int i = 1; i <= n / r; i++)
	{
		InsertSort(left + (i - 1)*r, left + i * r - 1); // 二次取中规则求每组的中间值
		Swap(left + i - 1, left + (i - 1)*r + (int)ceil((double)r / 2) - 1); // 将每组的中间值交
	}
	// 求二次中间值，其下标为 j
	int j = Select((int)ceil((double)n / r / 2), left, left + (n / r) - 1, r);
	Swap(left, j); // 二次中间值为枢纽元，并换至 left 处
	j = Partition(left, right); // 对表（子表）进行分划操作
	if (k == j - left + 1)
		return j; // 返回第 k小元素下标
	else if (k<j - left + 1)
		return Select(k, left, j - 1, r);// 在左子表求第 k小元素
	else return Select(k - (j - left + 1), j + 1, right, r); // 在右子表求第 k-(j-left+1) 小元素
}
void SortableList::InsertSort(int left, int right)
{
	for (int i = left + 1; i <= right; i++) {
		int j = i;
		int temp = l[i];
		while (j > left && temp < l[j - 1]) {
			l[j] = l[j - 1]; j--;
		}
		l[j] = temp;
	}
}
void main()
{
	int n = 10;
	int x = 4;
	SortableList myl(n);
	myl.Input();
	myl.Select(x, 4);
	myl.Output();
}