#pragma once

class UINsInvterm;
class SolveMRhs 
{
public:
	SolveMRhs();
	~SolveMRhs();
public:
	void Init();
	int RANKNUMBER;
	int NUMBER;
	int COLNUMBER;
	double residual;
	double* TempA;    //< The linear operator matrix
	int* TempIA;      //< The row number of the matrix
	int* TempJA;      //< The column number of the matrix
	double** TempB;
	double** TempX;
	void Deallocate();

public:
	void BGMRES();
};
extern SolveMRhs Rank;
extern SolveMRhs bgx;

