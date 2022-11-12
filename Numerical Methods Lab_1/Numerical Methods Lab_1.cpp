#include <iostream>
#include <fstream>
#include <vector>
#define pi 3.14
using namespace std;

struct _XY_ {
	double x;
	double y;
};

double func(double x) {
	return sin(4 * x * x);
}

double Chebyshev(int i, int n0, double a, double b) {
	return (b + a) / 2 + (b - a) / 2 * cos(pi * (2 * i + 1) / (2 * n0));
}

double Lagrange(int n, double X, vector<_XY_> XY) {
	double sigma = 0;
	double proiz;
	for (int i = 0; i < n; i++)
	{
		proiz = 1;
		for (int j = 0; j < n; j++)
		{
			if (j == i) continue;
			proiz *= (X - XY[j].x) / (XY[i].x - XY[j].x);
		}
		sigma += proiz * XY[i].y;
	}
	return sigma;
}

double Newton(int n,double X, vector<_XY_> x, vector<vector<double>> V) {

	for (int i = 1; i < n; i++) {
		for (int j = 0; j < n - i; j++) {
			V[i][j] = (V[i - 1][j + 1] - V[i - 1][j]) / (x[j + i].y - x[j].y);
		}
	}

	double S = V[0][0];
	double P = 1;
	for (int i = 1; i < n; ++i) {
		P *= (X - x[i - 1].y);
		S += V[i][0] * P;
	}

	return S;
}

int main(){
	const int size = 15;
	vector<_XY_> XY(size);
	double a = 0, b = 1;

	for (int i = 0; i < size; i++) {
		double z = i * (b - a) / size;
		XY[i].x = z;
		XY[i].y = func(z);
	}

	cout << "\n1.1 LAGRANGE\n" << endl;

	for (int i = 0; i < size; i++) {
		double z = i * (b - a) / size;
		cout << func(z) << " : " << Lagrange(size, z, XY) << endl;
	}

	double d = 10000;
	double delta, maxdelta = 0;
	ofstream f("Lagrange (1.1).txt");
	for (int i = 0; i < d; i++) {
		double z = i * (b - a) /d;
		delta = abs(func(z) - Lagrange(size, z, XY));
		f << i << "\t" << delta << endl;
		if (delta > maxdelta) {
			maxdelta = delta;
		}
	}

	cout << "\n" << "maxDelta = " << maxdelta << endl;

	f.close();

	cout << "\n1.2 CHEBYSHEV\n" << endl;

	vector<_XY_> XY2C(size);
	vector<_XY_> XY2L(size);

	for (int i = 0; i < size; i++) {
		XY2C[i].x = Chebyshev(i, size, a, b);
		XY2C[i].y = func(Chebyshev(i, size, a, b));
	}

	for (int i = 0; i < size; i++) {
		double z = i * (b - a) / size;
		XY2L[i].x = z;
		XY2L[i].y = func(z);
	}

	double  maxdeltaC = 0, maxdeltaL = 0;

	ofstream f2("Chebyshev (1.2.1).txt");
	for (int i = 0; i < d; i++) {
		double z = i * (b - a) / d;
		delta = abs(func(z) - Lagrange(size, z, XY2C));
		f2 << i << "\t" << delta << endl;
		if (delta > maxdeltaC) {
			maxdeltaC = delta;
		}
	}
	f2.close();

	ofstream f3("Lagrange (1.2.2).txt");
	for (int i = 0; i < d; i++) {
		double z = i * (b - a) / d;
		delta = abs(func(z) - Lagrange(size, z, XY2L));
		f3 << i << "\t" << delta << endl;
		if (delta > maxdeltaL) {
			maxdeltaL = delta;
		}
	}
	f3.close();

	cout << "\n" << "maxDeltaChebyshev = " << maxdeltaC << endl;
	cout << "\n" << "maxDeltaLagrange = " << maxdeltaL << endl;

	cout << "\n1.3 NEWTON\n" << endl;

	vector<vector<double>> V(size);

	V[0].resize(size);
	for (int i = 0; i < size; i++)
		V[0][i] = XY[i].y;

	for (int i = 1; i < size; i++) {
		V[i].resize(size);
		for (int j = 0; j < size; j++) {
			V[i][j] = 0;
		}
	}

	double  maxdeltaN = 0;

	ofstream f4("Newton (1.3).txt");
	for (int i = 0; i < d; i++) {
		double z = i * (b - a) / d;
		delta = abs(func(z) - Newton(size, z, XY, V));
		f4 << i << "\t" << delta << endl;
		if (delta > maxdeltaN) {
			maxdeltaN = delta;
		}
	}
	f4.close();

	cout << "\n" << "maxDeltaNewton = " << maxdeltaN << endl;

	cout << "\n1.4 TRIGONOMETRY\n" << endl;

	return 0;
}