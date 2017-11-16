// FiniteDiffMeth.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include "iostream" 
#include "fstream"
#include  "SolverPDE.h"

using namespace std;


int _tmain(int argc, _TCHAR* argv[])
{
	setlocale(LC_ALL, "Russian");
	ifstream fin("input.txt");
	SolverPDE spde;
	spde.getFileDate(fin);
	
	spde.solveByT();
	spde.showSolution();
	return 0;
}

