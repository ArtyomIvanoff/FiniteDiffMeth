#pragma once

#include "fstream"

using namespace std;

class SolverPDE {
public:
	SolverPDE();
	~SolverPDE(); //������������ ������� ������

	double p(double t, double x) const; // ������� ����� �� �� t(�����)
	double q(double t, double x) const;  // ������� ����� �� �� x(������������)
	double f(double t, double x) const;  // ������ ����� ����

	double b(double x) const; // ������ ����� ���������� �������
	double l(double t) const; // ������ ����� ������ ���������� �������
	double r(double t) const; // ������ ����� ������� ���������� �������

	void getFileDate(ifstream & fin); // ��������� �� ����� ��������� ����������
	void solveByT(); // �-�������� ���������� �����
	void printMtx(int i1) const;
	void tridiagSolve(int i1); //���������� ������ �������� ��� 3�-������������ �������
	void showSolution() const; //������� ������� �������
private:
	int n; // ����� ��������� �� t(�����)
	int m; // ����� ��������� �� x(������������)
	double hT; // ��� �� �������
	double hX; // ��� �� ������������
	double* t; // ��������� ������� �������
	double* x; // ��������� ������� ������������

	double** u; // ������� ������� �� "�����", �.�. ������� �������� ������, ����������� - (n+1)x(m+1)
	double** mtxKoef; // ������� ������������ ��� ������ ���������� ����, ����������� - (m+1)x(m+1)
	double* fVec; // ������� ��������� ������
};
