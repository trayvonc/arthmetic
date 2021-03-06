// 新密码算法.cpp: 定义控制台应用程序的入口点。
//
#include"stdafx.h"
#include<math.h>
#include<time.h>
#include<stdlib.h>
#include<iostream>

using namespace std;
int ext_euclid(int a, int b, int f, int e)
//因		e*0+φ*1=φ	(1)
//		e*1+φ*1=e	(2)
//则：e*(0-1*(φ/e))+.....=φ%e
//由于e和φ互质，因此一定有某一次运算后，等式右侧的φ%e==1。
//此时左侧等式中e所乘的系数就是所要求的d，即e^-1
//将(1)式e所乘的系数用a表示，(2)式中e所乘的系数用b表示，并且令：m=φ/e;  n=φ%e
{
	int m, n, t;
	if (e == 1) return b;
	m = f / e;	n = f % e;
	t = a - b * m;
	ext_euclid(b, t, e, n);
}

bool judge_prime(int a)//判断a是否为素数
{
	if (a <= 1) return false;
	if (a == 2) return true;
	int n = sqrt(a);
	for (int i = 2; i <= n; i++)
	{
		if (a%i == 0)
			return false;
	}
	return true;
}

int gcd(int m, int n)//辗转相除法求最大公约数
{
	int r = m % n;
	m = n;
	n = r;
	if (r == 0)
		return m;
	else
		return gcd(m, n);
}

bool judge_relativelyprime(int a, int b)//判断a,b是否互质
{
	if (gcd(a, b) == 1)
		return true;
	else
		return false;
}

int randomPrime(int max)
{
	srand((unsigned)time(NULL));
	int r = rand() % max;
	if (r < max / 2)//确保素数尽可能大
		r = max - r;
	if (judge_relativelyprime(r, max)) return r;
	for (int i = r; i < max; i++)
	{
		if (judge_relativelyprime(i, max))
			return i;
	}
	for (int i = r; i > 1; i--)
	{
		if (judge_relativelyprime(i, max))
			return i;
	}
}

int main()
{
	//输入质数p和q
	int p, q;
	do {
		cout << "输入一个质数p(如101):";
		cin >> p;
		cout << "输入一个质数q(如113):";
		cin >> q;
	} while (!(judge_prime(p) && judge_prime(q)));//p,q质数检测

												  //求得n=p*q的值
	int n = p * q;
	cout << "分组加密时，每个分组的大小不能超过n=p*q=";
	cout << n << endl;

	//求得φ(n)=(p-1)*(q-1)的值
	int f = (p - 1)*(q - 1);
	cout << "模φ(n)=(p-1)*(q-1)=";
	cout << f << endl << endl;

	//随机产生公钥e
	int e = randomPrime(f); //随机生成质数
	cout << "与φ(n)互质的公钥e为" << e << endl;

	//由e和φ(n)生成私钥d
	int d = ext_euclid(0, 1, f, e);
	while (d <= 0) d += f;
	cout << "通过调用扩展欧几里德算法，求得密钥d为：" << d << endl;

	//利用生成的公钥{e,n}对明文M进行加密
	int M, C;
	cout << "现在公钥{e,n}、私钥{d,n}均已生成完毕。\n\n请输入需要传输的明文内容进行加密(如9726)：";
	cin >> M;
	C = 1;
	for (int i = 1; i <= e; i++)
	{
		C = C * M%n;
	}
	cout << "明文M=" << M << "经加密后得到密文C=M^e(mod n)：" << C << endl;

	//利用生成的私钥私钥{e,n}对密文C进行解密
	M = 1;
	for (int i = 1; i <= d; i++)
	{
		M = M * C%n;
	}
	cout << "密文C=" << C << "经解密后得到明文M=C^d(mod n)：" << M << endl;
	
	system("pause");
	return 0;

}
