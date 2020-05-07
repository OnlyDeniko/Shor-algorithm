#pragma once

#include<iostream>
#include<vector>
#include<bitset>
#include<random>
#include<chrono>
#include<map>
#include<iomanip>
#include<cassert>
#include<omp.h>
#include"algoNumbers.h"
#include<complex>
#include<string>

#define deb(a) cout << #a << " = " << a << '\n';
using namespace std;
using cd = complex<double>;
extern const double eps;

template<class T>
bool operator <(const complex<T> & a, const complex<T> & b) {
	if (a.real() == b.real()) return a.imag() < b.imag();
	return a.real() < b.real();
}

class Quantum {
private:
	vector<cd> firstRegister;
	vector<cd> secondRegister;
	vector<cd> thirdRegister;
	cd fourthRegister;
	int a, N;
	int reverseNumber(int x, int len);
	string perevod(int x, int k);
	class BadGCD {};
	class BadAnswerReverse {};
	class BadAnswerShor {};
public:
	Quantum(int a, int N);
	vector<cd> QFFT(const vector<cd> & input, int cnt_qubits);
	vector<cd> ReverseQFFT(const vector<cd> & input, int cnt_qubits);

	vector<cd> Hadamar(const vector<cd > & a, int ind);
	vector<cd> Controlled_Z(const vector<cd > & a, int x, int y);
	vector<cd> CNOT(const vector<cd > & a, int x, int y);

	vector<cd> PhaseShift(const vector<cd> & a, int ind, int m);
	vector<cd> ReversePhaseShift(const vector<cd> & a, int ind, int m);

	vector<cd> Controlled_Rm(const vector<cd> & a, int x, int y, int m);
	vector<cd> ReverseControlled_Rm(const vector<cd> & a, int x, int y, int m);

	vector<cd> FADD(vector<cd> b, int a);
	vector<cd> ReverseFADD(vector<cd> b, int a);

	void FADD_modN(cd &_1, cd &_2, vector<cd> &b, cd &_4, int a, int N, int cnt_qub);
	void ReverseFADD_modN(cd &_1, cd &_2, vector<cd> &b, cd &_4, int a, int N, int cnt_qub);

	void CMULT_modN(cd &_1, vector<cd> & x, vector<cd> &b, cd &_4, int a, int N, int cnt_qub);
	void ReverseCMULT_modN(cd &_1, vector<cd> & x, vector<cd> &b, cd &_4, int a, int N, int cnt_qub);

	void UnitCMULT_modN(cd &_1, vector<cd> & x, vector<cd> & b, cd & _4, int a, int N, int cnt_qub);
	int Shor();

	vector<cd> DenseCoding(const vector<cd> & a, int cnt_qubits);
	vector<cd> OptimizedDenseCoding(const vector<cd> & a, int cnt_qubits);
	void TestDenseCoding(const vector<complex<double> > & input, int cnt_qubits, int cnt_tests);
	void TestOptimizedDenseCoding(const vector<complex<double> > & input, int cnt_qubits, int cnt_tests);

	int calc_whole_cubit(const vector<cd> & a);
	void test_for_whole_calc(const int &tests, const vector<complex<double> > & a);

	vector<cd> calc_one_cubit(const vector<cd> & a, int index);
	void test_for_one_calc(const int &tests, const vector<cd> & a, int index);
};