#include <iostream>
#include <fstream>
#define pi 3.14
using namespace std;

double func(double x) {
	return sin(4 * x * x);
}

double Chebyshev(int i, int n0, double a, double b) {
	return (b + a) / 2 + (b - a) / 2 * cos(pi * (2 * i + 1) / (2 * n0));
}

double Lagrange(int n, double X, double* x, double* y) {
	double sigma = 0;
	double proiz;
	for (int i = 0; i < n; i++)
	{
		proiz = 1;
		for (int j = 0; j < n; j++)
		{
			if (j == i) continue;
			proiz *= (X - x[j]) / (x[i] - x[j]);
		}
		sigma += proiz * y[i];
	}
	return sigma;
}
int main()
{
	const int size = 15;
	double x1[size];
	double y1[size];
	double a = 0, b = 1;

	for (int i = 0; i < size; i++) {
		double z = i * (b - a) / size;
		x1[i] = z;
		y1[i] = func(z);
	}

	cout << "\n1.1 LAGRANGE\n" << endl;

	for (int i = 0; i < size; i++) {
		double z = i * (b - a) / size;
		cout << func(z) << " : " << Lagrange(size, z, x1, y1) << endl;
	}

	double d = 10000;
	double delta, maxdelta = 0;

	ofstream fout("Lagrange (1.1).txt");
	for (int i = 0; i < d; i++) {
		double z = i * (b - a) /d;
		delta = abs(func(z) - Lagrange(size, z, x1, y1));
		fout << i << "\t" << delta << endl;
		if (delta > maxdelta) {
			maxdelta = delta;
		}
	}

	cout << "\n" << "maxDelta = " << maxdelta << endl;

	fout.close();

	cout << "\n1.2 CHEBYSHEV\n" << endl;

	const int size = 15;
	double x2C[size];
	double y2C[size];
	double x2L[size];
	double y2L[size];
	double a = 0, b = 1;

	for (int i = 0; i < size; i++) {
		double z = i * (b - a) / size;
		x2C[i] = func(z);
		y2C[i] = Chebyshev(i,size,a,b);
	}

	return 0;
}