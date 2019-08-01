#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include<fstream>
#include<iomanip>

using namespace std;
class SVD_Decomposition
{//U*A*V即奇异值标准型
public:
	SVD_Decomposition(vector<vector<double>>	&A_);
	~SVD_Decomposition();

public:
	void First_Step();
	void Householder(vector<double> &x);//返回x.size()+1维数组,最后一位存储Beta
	void Givens(vector<double> &x);//形参为(a,b),返回为(c,s)
	void Set_Zeros();
	bool Searching_For_Diagonal_Matrix();
	void If_Diagonal_Zero(int k);//A[k][k]=0
	void Wilkinson_Step();
	void Decomposition();
	void Modify_Negatives();
	void Sort();
	void Print_To_File();
	//void Error_Estimation();//计算

public:
	vector<vector<double>>	A;
	vector<vector<double>>  A_Initial;
	vector<vector<double>>	U;
	vector<vector<double>>	V;

	int p, q;//意义参见课本算法7.6.3
};

