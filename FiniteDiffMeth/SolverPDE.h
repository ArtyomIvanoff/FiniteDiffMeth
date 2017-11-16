#pragma once

#include "fstream"

using namespace std;

class SolverPDE {
public:
	SolverPDE();
	~SolverPDE(); //освобождение занятой памяти

	double p(double t, double x) const; // функция перед ЧП по t(время)
	double q(double t, double x) const;  // функция перед ЧП по x(пространство)
	double f(double t, double x) const;  // правая часть УСЧП

	double b(double x) const; // правая часть начального условия
	double l(double t) const; // правая часть левого граничного условия
	double r(double t) const; // правая часть правого граничного условия

	void getFileDate(ifstream & fin); // извлекаем из файла начальную информацию
	void solveByT(); // Т-образная разностная схема
	void printMtx(int i1) const;
	void tridiagSolve(int i1); //реализация метода прогонки для 3х-диагональной матрицы
	void showSolution() const; //вывести матрицу решения
private:
	int n; // число разбиений по t(время)
	int m; // число разбиений по x(пространство)
	double hT; // шаг по времени
	double hX; // шаг по пространству
	double* t; // разбиение отрезка времени
	double* x; // разбиение отрезка пространства

	double** u; // матрица решения на "сетке", т.е. решение сеточной задачи, размерность - (n+1)x(m+1)
	double** mtxKoef; // матрица коэффицентов для одного временного слоя, размерность - (m+1)x(m+1)
	double* fVec; // столбец свободных членов
};
