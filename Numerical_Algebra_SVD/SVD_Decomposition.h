#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include<fstream>
#include<iomanip>

using namespace std;
class SVD_Decomposition
{//U*A*V������ֵ��׼��
public:
	SVD_Decomposition(vector<vector<double>>	&A_);
	~SVD_Decomposition();

public:
	void First_Step();
	void Householder(vector<double> &x);//����x.size()+1ά����,���һλ�洢Beta
	void Givens(vector<double> &x);//�β�Ϊ(a,b),����Ϊ(c,s)
	void Set_Zeros();
	bool Searching_For_Diagonal_Matrix();
	void If_Diagonal_Zero(int k);//A[k][k]=0
	void Wilkinson_Step();
	void Decomposition();
	void Modify_Negatives();
	void Sort();
	void Print_To_File();
	//void Error_Estimation();//����

public:
	vector<vector<double>>	A;
	vector<vector<double>>  A_Initial;
	vector<vector<double>>	U;
	vector<vector<double>>	V;

	int p, q;//����μ��α��㷨7.6.3
};

