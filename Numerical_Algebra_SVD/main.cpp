#include<iostream>
#include<stdio.h>
#include<iomanip>
#include<stdlib.h>
#include<time.h>
using namespace std;
#include "SVD_Decomposition.h"

void print(vector<vector<double>> &A)
{
	cout << endl;
	for (size_t i = 0; i < A.size(); i++)
	{
		for (size_t j = 0; j < A[0].size(); j++)
		{
			cout << setprecision(5);
			cout.width(12);
			if (abs(A[i][j]) < 1.0e-6)
				cout << 0.0 << " ";
			else
				cout << A[i][j] << " ";
		}
		cout << endl;
	}
}

//110*100阶随机数矩阵
void Test2()
{

	int M = 110, N = 100;
	vector<vector<double>> A(M);
	for (size_t i = 0; i < M; i++)
	{
		A[i].resize(N);
		for (size_t j = 0; j < N; j++)
		{
			A[i][j] = rand()%100;
		}
	}
	clock_t start = clock();
	SVD_Decomposition* X = new SVD_Decomposition(A);
	X->Decomposition();
	clock_t end = clock();

	cout << "奇异值:" << endl;
	for (int i = 0; i < X->A[0].size(); i++)
	{
		cout << X->A[i][i] << " ";
		if (!((i + 1) % 10))cout << endl;
	}
	cout << "耗时" << (end - start) << "毫秒" << endl;
	cout << endl;

	cout << "U.Transpose[U]-I 的2-范数:";
	double Norm = 0;
	for (int i = 0; i < X->U.size(); i++)
	{
		for (int j = 0; j < X->U.size(); j++)
		{
			double temp = 0;
			for (int k = 0; k < X->U.size(); k++)
			{
				temp += X->U[i][k] * X->U[j][k];
			}
			if (i == j)
				Norm += (temp - 1)*(temp - 1);
			else
				Norm += temp * temp;
		}
	}
	cout << sqrt(Norm) << endl;

	cout << "V.Transpose[V]-I 的2-范数:";
	Norm = 0;
	for (int i = 0; i < X->V.size(); i++)
	{
		for (int j = 0; j < X->V.size(); j++)
		{
			double temp = 0;
			for (int k = 0; k < X->V.size(); k++)
			{
				temp += X->V[k][i] * X->V[k][j];
			}
			if (i == j)
				Norm += (temp - 1)*(temp - 1);
			else
				Norm += temp * temp;
		}
	}
	cout << sqrt(Norm) << endl;

	cout << "U.A.Transpose[V]-T (A原矩阵,T奇异值标准型)的2-范数:";
	vector<vector<double>> S(A.size());
	for (int i = 0; i < A.size(); i++)
	{
		S[i].resize(A[0].size());
		for (int j = 0; j < A[0].size(); j++)
		{
			S[i][j] = 0;
			for (int k = 0; k < A.size(); k++)
			{
				S[i][j] += X->U[i][k] * A[k][j];
			}
		}
	}
	vector<vector<double>> T(S);
	for (int i = 0; i < S.size(); i++)
	{
		for (int j = 0; j < S[0].size(); j++)
		{
			S[i][j] = 0;
			for (int k = 0; k < S[0].size(); k++)
			{
				S[i][j] += T[i][k] * X->V[j][k];
			}
		}
	}
	T.clear();
	Norm = 0;
	for (int i = 0; i < S.size(); i++)
	{
		for (int j = 0; j < S[0].size(); j++)
		{
			Norm += (S[i][j] - X->A[i][j])*(S[i][j] - X->A[i][j]);
		}
	}
	cout << sqrt(Norm) << endl;

	delete X;
	A.clear();
}

void Test1()
{
	int M = 28, N = 12;
	vector<vector<double>> A(M);
	double A_[28][12] = { {1, 4.9176, 1, 3.472, 0.998, 1, 7, 4, 42, 3, 1, 0}, {1, 5.0208, 1,
   3.531, 1.5, 2, 7, 4, 62, 1, 1, 0}, {1, 4.5429, 1, 2.275, 1.175, 1,
   6, 3, 40, 2, 1, 0}, {1, 4.5573, 1, 4.05, 1.232, 1, 6, 3, 54, 4, 1,
   0}, {1, 5.0597, 1, 4.455, 1.121, 1, 6, 3, 42, 3, 1, 0}, {1, 3.891,
   1, 4.455, 0.988, 1, 6, 3, 56, 2, 1, 0}, {1, 5.898, 1, 5.85, 1.24,
   1, 7, 3, 51, 2, 1, 1}, {1, 5.6039, 1, 9.52, 1.501, 0, 6, 3, 32, 1,
   1, 0}, {1, 15.4202, 2.5, 9.8, 3.42, 2, 10, 5, 42, 2, 1, 1}, {1,
   14.4598, 2.5, 12.8, 3, 2, 9, 5, 14, 4, 1, 1}, {1, 5.8282, 1, 6.435,
	1.225, 2, 6, 3, 32, 1, 1, 0}, {1, 5.3003, 1, 4.9883, 1.552, 1, 6,
   3, 30, 1, 2, 0}, {1, 6.2712, 1, 5.52, 0.975, 1, 5, 2, 30, 1, 2,
   0}, {1, 5.9592, 1, 6.666, 1.121, 2, 6, 3, 32, 2, 1, 0}, {1, 5.05,
   1, 5, 1.02, 0, 5, 2, 46, 4, 1, 1}, {1, 5.6039, 1, 9.52, 1.501, 0,
   6, 3, 32, 1, 1, 0}, {1, 8.2464, 1.5, 5.15, 1.664, 2, 8, 4, 50, 4,
   1, 0}, {1, 6.6969, 1.5, 6.092, 1.488, 1.5, 7, 3, 22, 1, 1, 1}, {1,
   7.7841, 1.5, 7.102, 1.376, 1, 6, 3, 17, 2, 1, 0}, {1, 9.0384, 1,
   7.8, 1.5, 1.5, 7, 3, 23, 3, 3, 0}, {1, 5.9894, 1, 5.52, 1.256, 2,
   6, 3, 40, 4, 1, 1}, {1, 7.5422, 1.5, 4, 1.69, 1, 6, 3, 22, 1, 1,
   0}, {1, 8.7951, 1.5, 9.89, 1.82, 2, 8, 4, 50, 1, 1, 1}, {1, 6.0931,
	1.5, 6.7265, 1.652, 1, 6, 3, 44, 4, 1, 0}, {1, 8.3607, 1.5, 9.15,
   1.777, 2., 8, 4, 48, 1, 1, 1}, {1, 8.14, 1, 8, 1.504, 2, 7, 3, 3,
   1, 3, 0}, {1, 9.1416, 1.5, 7.3262, 1.831, 1.5, 8, 4, 31, 4, 1,
   0}, {1, 12, 1.5, 5, 1.2, 2, 6, 3, 30, 3, 1, 1} };
	for (size_t i = 0; i < M; i++)
	{
		A[i].resize(N);
		for (size_t j = 0; j < N; j++)
		{
			A[i][j] = A_[i][j];
		}
	}
	//cout << "A:" << endl;
	//print(A);
	clock_t start = clock();
	SVD_Decomposition* X = new SVD_Decomposition(A);
	X->Decomposition();
	clock_t end = clock();
	//cout << "T:" << endl;
	//print(X->A);
	//cout << "U:" << endl;
	//print(X->U);
	//cout << "V:" << endl;
	//print(X->V);
	X->Print_To_File();

	cout << "奇异值:" << endl;
	for (int i = 0; i < X->A[0].size(); i++)
		cout << X->A[i][i] << endl;
	cout << "耗时" << (end - start) << "毫秒" << endl;
	cout << endl;

	cout << "U.Transpose[U]-I 的2-范数:";
	double Norm = 0;
	for (int i = 0; i < X->U.size(); i++)
	{
		for (int j = 0; j < X->U.size(); j++)
		{
			double temp = 0;
			for (int k = 0; k < X->U.size(); k++)
			{
				temp += X->U[i][k] * X->U[j][k];
			}
			if (i == j)
				Norm += (temp - 1)*(temp - 1);
			else
				Norm += temp * temp;
		}
	}
	cout << sqrt(Norm) << endl;

	cout << "V.Transpose[V]-I 的2-范数:";
	Norm = 0;
	for (int i = 0; i < X->V.size(); i++)
	{
		for (int j = 0; j < X->V.size(); j++)
		{
			double temp = 0;
			for (int k = 0; k < X->V.size(); k++)
			{
				temp += X->V[k][i] * X->V[k][j];
			}
			if (i == j)
				Norm += (temp - 1)*(temp - 1);
			else
				Norm += temp * temp;
		}
	}
	cout << sqrt(Norm) << endl;

	cout << "U.A.Transpose[V]-T (A原矩阵,T奇异值标准型)的2-范数:";
	vector<vector<double>> S(A.size());
	for (int i = 0; i < A.size(); i++)
	{
		S[i].resize(A[0].size());
		for (int j = 0; j < A[0].size(); j++)
		{
			S[i][j] = 0;
			for (int k = 0; k < A.size(); k++)
			{
				S[i][j] += X->U[i][k] * A[k][j];
			}
		}
	}
	vector<vector<double>> T(S);
	for (int i = 0; i < S.size(); i++)
	{
		for (int j = 0; j < S[0].size(); j++)
		{
			S[i][j] = 0;
			for (int k = 0; k < S[0].size(); k++)
			{
				S[i][j] += T[i][k] * X->V[j][k];
			}
		}
	}
	T.clear();
	Norm = 0;
	for (int i = 0; i < S.size(); i++)
	{
		for (int j = 0; j < S[0].size(); j++)
		{
			Norm += (S[i][j] - X->A[i][j])*(S[i][j] - X->A[i][j]);
		}
	}
	cout << sqrt(Norm) << endl;

	delete X;
	A.clear();
}

void main()
{
	Test1();
	//Test2();
}