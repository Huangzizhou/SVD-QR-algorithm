#include "SVD_Decomposition.h"



SVD_Decomposition::SVD_Decomposition(vector<vector<double>>	&A_)
{
	if (A_.size() < A_[0].size())
		cout << "行比列少！" << endl;
	A=A_;
	A_Initial = A;
	U.resize(A.size());
	V.resize(A[0].size());
	for (size_t i = 0; i < A.size(); i++)
		U[i].resize(A.size());
	for (size_t i = 0; i < A[0].size(); i++)
		V[i].resize(A[0].size());
	for (size_t i = 0; i < A.size(); i++)
		U[i][i] = 1;
	for (size_t i = 0; i < A[0].size(); i++)
		V[i][i] = 1;
	q = 0;
	p = A[0].size();
}

SVD_Decomposition::~SVD_Decomposition()
{
	A.clear();
	U.clear();
	V.clear();
}

void SVD_Decomposition::First_Step()
{
	int n = A[0].size(), m = A.size();
	for (size_t k = 1; k <= n; k++)
	{	
		//打第k列
		if (n < m || k < n)
		{//方阵不用打最后一列
			vector<double> x(m - k + 1);
			for (size_t i = k - 1; i < m; i++)
				x[i - k + 1] = A[i][k - 1];
			Householder(x);
			vector<vector<double>> T(m - k + 1);
			for (size_t i = 0; i < m - k + 1; i++)
				T[i].resize(n);
			for (size_t i = 0; i < m - k + 1; i++)
			{
				for (size_t j = 0; j < n; j++)
				{
					for (size_t t = 0; t < m - k + 1; t++)
					{
						T[i][j] += x[t] * A[t + k - 1][j];
					}
					T[i][j] *= x[i] * x[m - k + 1];
				}
			}
			for (size_t i = 0; i < m - k + 1; i++)
			{
				for (size_t j = 0; j < n; j++)
				{
					A[i + k - 1][j] -= T[i][j];
				}
			}
			T.clear();
			T.resize(m - k + 1);
			for (size_t i = 0; i < m - k + 1; i++)
				T[i].resize(m);
			for (size_t i = 0; i < m - k + 1; i++)
			{
				for (size_t j = 0; j < m; j++)
				{
					for (size_t t = 0; t < m - k + 1; t++)
					{
						T[i][j] += x[t] * U[t + k - 1][j];
					}
					T[i][j] *= x[i] * x[m - k + 1];
				}
			}
			for (size_t i = 0; i < m - k + 1; i++)
			{
				for (size_t j = 0; j < m; j++)
				{
					U[i + k - 1][j] -= T[i][j];
				}
			}
			x.clear();
		}
		//打第k行
		if (k < n - 1)
		{//不用打最后两行
			vector<double>	x(n - k);
			for (size_t j = k; j < n; j++)
				x[j-k] = A[k - 1][j];
			Householder(x);
			vector<vector<double>>	T(m);
			for (size_t i = 0; i < m; i++)
				T[i].resize(n - k);
			for (size_t i = 0; i < m; i++)
			{
				for (size_t j = 0; j < n - k; j++)
				{
					for (size_t t = 0; t < n - k; t++)
					{
						T[i][j] += A[i][t + k] * x[t];
					}
					T[i][j] *= x[j] * x[n - k];
				}
			}
			for (size_t i = 0; i < m; i++)
			{
				for (size_t j = 0; j < n - k; j++)
				{
					A[i][j + k] -= T[i][j];
				}
			}
			T.clear();
			T.resize(n);
			for (size_t i = 0; i < n; i++)
				T[i].resize(n - k);
			for (size_t i = 0; i < n; i++)
			{
				for (size_t j = 0; j < n - k; j++)
				{
					for (size_t t = 0; t < n - k; t++)
					{
						T[i][j] += V[t+k][i] * x[t];
					}
					T[i][j] *= x[j] * x[n - k];
				}
			}
			for (size_t i = 0; i < n; i++)
			{
				for (size_t j = 0; j < n - k; j++)
				{
					V[j+k][i] -= T[i][j];
				}
			}
			x.clear();
		}

	}
	for (size_t i = 0; i < m; i++)
	{
		for (size_t j = 0; j < n; j++)
		{
			if (i != j && i != j - 1)
				A[i][j] = 0;
		}
	}
}

void SVD_Decomposition::Householder(vector<double> &x)
{
	double Eta = 0, Sigma = 0;
	vector<double> v(x);
	for (size_t i = 0; i < x.size(); i++)
	{
		if (Eta < abs(x[i])) Eta = abs(x[i]);
	}
	for (size_t i = 0; i < x.size(); i++)
	{
		x[i] /= Eta;
		if (i)
		{
			Sigma += x[i] * x[i];
			v[i] = x[i];
		}
	}
	if (!Sigma)
	{
		if (x[0] < 0)
		{
			v[0] = 1;
			v.push_back(2);
		}
		else
		{
			v[0] = 0;
			v.push_back(0);
		}
	}
	else
	{
		double Alpha = sqrt(x[0] * x[0] + Sigma);
		if (x[0] <= 0)
		{
			v[0] = x[0] - Alpha;
		}
		else
		{
			v[0] = -Sigma / (x[0] + Alpha);
		}
		v.push_back(2 * v[0] * v[0] / (Sigma + v[0] * v[0]));
		double temp = v[0];
		for (size_t i = 0; i < x.size(); i++)
		{
			v[i] /= temp;
		}
	}
	x = v;
}

void SVD_Decomposition::Set_Zeros()
{
	double Infinity_Norm = 0;
	Infinity_Norm = abs(A[0][0]) > Infinity_Norm ? abs(A[0][0]) : Infinity_Norm;
	for (size_t i = 0; i < A[0].size() - 1; i++)
	{
		if (abs(A[i][i + 1]) < (abs(A[i][i])+abs(A[i+1][i+1]))*1.0e-20)
		{
			A[i][i + 1] = 0;
		}
		Infinity_Norm = abs(A[i][i + 1]) > Infinity_Norm ? abs(A[i][i + 1]) : Infinity_Norm;
		Infinity_Norm = abs(A[i + 1][i + 1]) > Infinity_Norm ? abs(A[i + 1][i + 1]) : Infinity_Norm;
	}
	for (size_t i = 0; i < A[0].size(); i++)
	{
		if (abs(A[i][i]) < Infinity_Norm*1.0e-20)
			A[i][i] = 0;
	}
}

bool SVD_Decomposition::Searching_For_Diagonal_Matrix()
{
	int n = A[0].size();
	int k;
	for (k = n - 1; k > 0; k--)
	{
		if (A[k - 1][k] != 0)
			break;
	}
	if (!k)
		return true;
	q = n - k - 1;
	for (; k > 0; k--)
	{
		if (A[k - 1][k] == 0)
			break;
	}
	p = k;
	return false;
}

void SVD_Decomposition::Givens(vector<double> &x)
{
	if (!x[1])
	{
		x[0] = 1;
	}
	else
	{
		if (abs(x[1]) > abs(x[0]))
		{
			double Tau = x[0] / x[1];
			x[1] = 1 / sqrt(1 + Tau * Tau);
			x[0] = x[1] * Tau;
		}
		else
		{
			double Tau = x[1] / x[0];
			x[0] = 1 / sqrt(1 + Tau * Tau);
			x[1] = x[0] * Tau;
		}
	}
}

void SVD_Decomposition::If_Diagonal_Zero(int k)
{
	int n = A[0].size(), m = A.size();
	int min = p, max = n - q - 1;//B22的边界行
	for (size_t t = k + 1; t <= max; t++)
	{
		vector<double> x(2);
		double a = A[k][t], b = A[t][t];
		x[0] = a;	x[1] = b;
		Givens(x);
		A[k][t] = x[0] * a + x[1] * b;
		A[t][t] = -x[1] * a + x[0] * b;
		if (t < max)
		{
			a = A[k][t + 1];	b = A[t][t + 1];
			A[k][t+1] = x[0] * a + x[1] * b;
			A[t][t + 1] = -x[1] * a + x[0] * b;
		}
		
		for (size_t j = 0; j < m; j++)
		{
			a = U[k][j];	b = U[t][j];
			U[k][j] = x[0] * a + x[1] * b;
			U[t][j] = -x[1] * a + x[0] * b;
		}
	}
}

void SVD_Decomposition::Wilkinson_Step()
{//希腊字母意义均与算法7.6.2一致
	int n = A[0].size(), m = A.size();
	int min = p, max = n - q - 1;//B22的边界行
	double Alpha = A[max][max] * A[max][max] + A[max - 1][max] * A[max - 1][max];
	double Delta;
	if (max - min > 1)
		Delta = A[max - 1][max - 1] * A[max - 1][max - 1] + A[max - 2][max - 1] * A[max - 2][max - 1];
	else if (max == min + 1)
		Delta = A[max - 1][max - 1] * A[max - 1][max - 1];
	else return;
		Delta = (Delta - Alpha) / 2;
		double Beta = A[max - 1][max - 1] * A[max - 1][max];
		double Miu;
		if (Delta >= 0)
			Miu = Delta + sqrt(Delta*Delta + Beta * Beta);
		else
			Miu = Delta - sqrt(Delta*Delta + Beta * Beta);
		Miu = Alpha - Beta * Beta / Miu;
		double y = A[min][min] * A[min][min] - Miu;
		double z = A[min][min] * A[min][min + 1];

		for (size_t k = 1; k < max-min+1; k++)
		{
			vector<double> G(2);
			G[0] = y; G[1] = z;
			Givens(G);
			double c = G[0], s = -G[1];
			if (k > 1)
				A[min + k - 2][min + k - 1] = c * y - s * z;

			y = c * A[min + k - 1][min + k - 1] - s * A[min + k - 1][min + k];
			z = -s * A[min + k][min + k];
			A[min + k - 1][min + k] = s * A[min + k - 1][min + k - 1] + c * A[min + k - 1][min + k];
			A[min + k][min + k] *= c;

			for (size_t i = 0; i < n; i++)
			{
				double a = V[min + k - 1][i], b = V[min + k][i];
				V[min + k - 1][i] = c * a - s * b;
				V[min + k][i] = s * a + c * b;
			}

			G[0] = y; G[1] = z;
			Givens(G);
			c = G[0], s = -G[1];
			A[min + k - 1][min + k - 1] = c * y - s * z;

			if (k < max - min)
			{
				y = c * A[min + k - 1][min + k] - s * A[min + k][min + k];
				z = -s * A[min + k][min + k + 1];
				A[min + k][min + k + 1] *= c;
				A[min + k][min + k] = s * A[min + k - 1][min + k] + c * A[min + k][min + k];
			}
			else
			{
				double a = A[min + k - 1][min + k], b = A[min + k][min + k];
				A[min + k - 1][min + k] = c * a - s * b;
				A[min + k][min + k] = s * a + c * b;
			}
			for (size_t j = 0; j < m; j++)
			{
				double a = U[min + k - 1][j], b = U[min + k][j];
				U[min + k - 1][j] = c * a - s * b;
				U[min + k][j] = s * a + c * b;
			}
		}
}

void SVD_Decomposition::Decomposition()
{
	First_Step();
	while (1)
	{
		Set_Zeros();
		if (Searching_For_Diagonal_Matrix())	break;
		int min = p, max = A[0].size() - q - 1;//B22的边界行
		bool flag = false;
		for (int k = max - 1; k >= min; k--)
		{
			if (!A[k][k])
			{
				flag = true;
				If_Diagonal_Zero(k);
				break;
			}
		}
		if (flag)	continue;
		
		vector<double> delta(A[0].size()), gamma(A[0].size() - 1);
		for (size_t i = 0; i < A[0].size(); i++)
		{
			delta[i] = A[i][i];
			if (i > 0)
				gamma[i-1] = A[i - 1][i];
		}

		Wilkinson_Step();
	}
	Modify_Negatives();
	Sort();
}

void SVD_Decomposition::Print_To_File()
{
	ofstream ofile;               //定义输出文件
	ofile.open("myfile.txt");     //作为输出文件打开
	ofile.precision(10);
	ofile.setf(ios::fixed);
	ofile << "A={";
	for (size_t i = 0; i < A_Initial.size(); i++)
	{
		ofile << "{";
		for (size_t j = 0; j < A_Initial[0].size(); j++)
		{
			ofile << A_Initial[i][j];
			if (j < A_Initial[0].size() - 1)
				ofile << ",";
		}
		ofile << "}";
		if (i < A_Initial.size() - 1)
			ofile << ",";
	}
	ofile << "};" << endl;
	ofile << "T={";
	for (size_t i = 0; i < A.size(); i++)
	{
		ofile << "{";
		for (size_t j = 0; j < A[0].size(); j++)
		{
			ofile << A[i][j];
			if (j < A[0].size() - 1)
				ofile << ",";
		}
		ofile << "}";
		if (i < A.size() - 1)
			ofile << ",";
	}
	ofile << "};" << endl;
	ofile << "U={";
	for (size_t i = 0; i < A.size(); i++)
	{
		ofile << "{";
		for (size_t j = 0; j < A.size(); j++)
		{
			ofile << U[i][j];
			if (j < A.size() - 1)
				ofile << ",";
		}
		ofile << "}";
		if (i < A.size() - 1)
			ofile << ",";
	}
	ofile << "};" << endl;
	ofile << "V={";
	for (size_t i = 0; i < A[0].size(); i++)
	{
		ofile << "{";
		for (size_t j = 0; j < A[0].size(); j++)
		{
			ofile << V[i][j];
			if (j < A[0].size() - 1)
				ofile << ",";
		}
		ofile << "}";
		if (i < A[0].size() - 1)
			ofile << ",";
	}
	ofile << "};";
	ofile.close();                //关闭文件
}

void SVD_Decomposition::Modify_Negatives()
{
	for (int i = 0; i < A[0].size(); i++)
	{
		if (A[i][i] < 0)
		{
			A[i][i] *= -1;
			for (int k = 0; k < A[0].size(); k++)
			{
				V[i][k] *= -1;
			}
		}
	}
}

void SVD_Decomposition::Sort()
{
	for (int i = 0; i < A[0].size() - 1; i++)
	{
		double max = i;
		for (int j = i + 1; j < A[0].size(); j++)
		{
			if (A[max][max] < A[j][j])
				max = j;
		}
		double temp = A[i][i];
		A[i][i] = A[max][max];
		A[max][max] = temp;

		vector<double> temp_ = U[i];
		U[i] = U[max];
		U[max] = temp_;

		temp_ = V[i];
		V[i] = V[max];
		V[max] = temp_;
	}
}