#pragma once
#include<complex>
#include<cmath>
#include<vector>
#include"Tmath.h"

using namespace std;

using cd = complex<double>;

class Quantum{
private:
	static const double PI;
	int a, N;
	double noize;
	vector<cd> data;
	int logSize;
	vector<int> counter;

	bool checkBit(int mask, int bit);
	void Hadamar(int pos);
	void NOT(int pos);

	void CNOT(int pos, int m);
	void CCNOT(int pos1, int pos2, int m);

	void PhaseShift(int x, double alp);
	void CPhaseShift(int x, int y, double alp);
	void CCPhaseShift(int x, int y, int z, double alp);

	void QFFT(int l, int r);
	void ReverseQFFT(int l, int r);

	void ADD(int a, int l, int r, int znak);
	void CADD(int a, int l, int r, int m, int znak);
	void CCADD(int a, int l, int r, int m1, int m2, int znak);
	void CCADD_modN(int a, int l, int r, int m1, int m2);

	void CMULT(int a, int l, int r, int m);
	void ReverseCMULT(int a, int l, int r, int m);

	void CSWAP(int l, int r, int m);

	void UnitCMULT(int a, int l, int r, int m);

	void Shor(int l, int r);
public:
	Quantum(int a = 2, int N = 15, double noize = 0.0);
	vector<cd> getData();
	int getQubitNumber();
	vector<int> getCounter();
	int measure_qubit();
	void measures(int times);
	int getPeriod(int measuresTime = 1000);
};
