#include "stdafx.h"
#include "iostream"
#include "fstream"
#include "SolverPDE.h"

using namespace std;

SolverPDE::SolverPDE() : n(0), m(0), t(NULL), x(NULL), u(NULL), mtxKoef(NULL), fVec(NULL) {}

SolverPDE::~SolverPDE() {
	delete[] t;
	delete[] x;

	for (int i = 0; i < n; i++)
		delete[] u[i];
	delete[] u;

	for (int i = 0; i < m + 1; i++)
		delete[] mtxKoef[i];
	delete[] mtxKoef;

	delete[] fVec;
}

double SolverPDE::p(double t, double x) const {
	return 1.0;
}

double SolverPDE::q(double t, double x) const {
	return 1.0;
}

double SolverPDE::f(double t, double x) const {
	return 0.0;
}

double SolverPDE::b(double x) const {
	return 3 * x;
}

double SolverPDE::l(double t) const {
	return 0;
}

double SolverPDE::r(double t) const {
	return 0;
}

void SolverPDE::getFileDate(ifstream & fin) {
	double t0, T; // начальная и конечная точка по времени
	fin >> t0 >> T;

	double x0, X; // начальная и конечная точка по пространству
	fin >> x0 >> X;

	fin >> n >> m; // числа разбиений по времени и пространству
	fin.close();

	hT = (T - t0) / n;
	hX = (X - x0) / m;

	t = new double[n + 1];
	t[0] = t0;
	t[m] = T;
	for (int i = 1; i < n; i++)
		t[i] = t[i - 1] + hT;

	x = new double[m + 1];
	x[0] = x0;
	x[m] = X;
	for (int i = 1; i < m; i++)
		x[i] = x[i - 1] + hX;
}

void SolverPDE::solveByT() {
	u = new double*[n + 1];
	for (int i = 0; i < n + 1; i++)
		u[i] = new double[m + 1];

	mtxKoef = new double*[m + 1];
	for (int i = 0; i < m + 1; i++)
		mtxKoef[i] = new double[m + 1];

	// заполняем нулями элементы матрицы mtxKoef
	for (int i = 0; i < m + 1; i++){
		for (int j = 0; j < m + 1; j++) {
			mtxKoef[i][j] = 0.0;
		}
	}
	// верхняя и нижняя строка матрицы кожффмцентов будет неизменными
	mtxKoef[0][0] = mtxKoef[m][m] = 1.0;

	fVec = new double[m + 1];

	for (int j = 0; j < m + 1; j++) // используем начальные условия
		u[0][j] = b(x[j]);

	for (int i = 1; i < n + 1; i++) { // используем краевые условия
		u[i][0] = l(t[i]);
		u[i][m] = r(t[i]);
	}

	double pij, qij;

	for (int i = 1; i < n + 1; i++) { //на каждом шаге формируем систему для очередного временного слоя
		fVec[0] = u[i][0];
		for (int j = 1; j < m; j++){
			pij = p(t[i], x[j]);
			qij = q(t[i], x[j]);

			mtxKoef[j][j - 1] = (-1)*qij / (2 * hX);
			mtxKoef[j][j] = pij / hT;
			mtxKoef[j][j + 1] = qij / (2 * hX);
			fVec[j] = f(t[i], x[j]) + u[i - 1][j] / hT;
		}
		fVec[m] = u[i][m];

		printMtx(i);
		tridiagSolve(i);
	}
}

void SolverPDE::printMtx(int i1) const {
	cout << "Получаемая "  << i1 <<" система: \n";
	for (int i = 0; i < m + 1; i++) {
		for (int j = 0; j < m+ 1; j++) {
			cout << mtxKoef[i][j] << " ";
		}
		cout << " | " << fVec[i] << endl;
	}
}

void SolverPDE::tridiagSolve(int i1){
	double *a = new double[m + 1];
	a[0] = 0;

	for (int i = 1; i < m + 1; i++) 
		a[i] = mtxKoef[i][i - 1]; // вылетает здесь программа
		
	double *b = new double[m + 1];
	for (int i = 0; i < n + 1; i++)
		b[i] = mtxKoef[i][i];

	double *c = new double[m + 1];
	for (int i = 0; i < m; i++)
		c[i] = mtxKoef[i][i + 1];
	c[n] = 0;

	double *d = new double[m + 1];
	for (int i = 0; i < m + 1; i++)
		d[i] = fVec[i];
	
	double *A = new double[m + 1];
	double *B = new double[m + 1];

	A[0] = (-1)*c[0] / b[0];
	B[0] = d[0] / b[0];

	for (int i = 1; i < m + 1; i++) {
		A[i] = (-1)*c[i] / (b[i] + a[i] * A[i - 1]);
		B[i] = (d[i] - a[i] * B[i - 1]) / (b[i] + a[i] * A[i - 1]);
	}

	delete[] a;
	delete[] b;
	delete[] c;
	delete[] d;
	
	// Решение данной системы записываем в i-ый временнной слой
	u[i1][m] = B[m];

	for (int i = m - 1; i >= 0; i--)
		u[i1][i] = A[i] * u[i1][i + 1] + B[i];

	delete[] A;
	delete[] B;
}

void SolverPDE::showSolution() const {
	cout << "\nРешение сеточной задачи: \n";
	for (int i = 0; i < n + 1; i++) {
		cout << "Слой t = " << t[i] << endl;
		for (int j = 0; j < m + 1; j++) {
			cout << "  u(" << x[j] << ") = " << u[i][j] << endl;
		}
		cout << endl;
	}
}