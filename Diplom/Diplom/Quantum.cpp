#include "Quantum.h"
#include<iostream>

const double Quantum::PI = acos(-1.0);

bool Quantum::checkBit(int mask, int bit){
	return mask & (1 << bit);
}

void Quantum::Hadamar(int pos){
	int shift = 1 << pos;
	int i = 0;
#pragma omp parallel for private(i)
	for (i = 0; i < (int)data.size(); i++)
	{
		if (!checkBit(i % (1 << (pos + 1)), pos)){
			cd q = data[i];
			cd w = data[i + shift];
			data[i] = (q + w) / sqrt(2);
			data[i + shift] = (q - w) / sqrt(2);
		}
	}
}

void Quantum::NOT(int pos) {
	int shift = 1 << pos;
	int i = 0;
#pragma omp parallel for private(i)
	for (i = 0; i < (int)data.size(); i++)
	{
		if (!checkBit(i % (1 << (pos + 1)), pos)) {
			swap(data[i], data[i + shift]);
		}
	}
}

void Quantum::CNOT(int pos, int m){
	int shift = 1 << m;
	int i = 0;
#pragma omp parallel for private(i)
	for (i = 0; i < (int)data.size(); i++) {
		if (checkBit(i % (1 << (pos + 1)), pos) && !checkBit(i % (1 << (m + 1)), m)) {
			swap(data[i], data[i + shift]);
		}
	}
}

void Quantum::CCNOT(int pos1, int pos2, int m){
	int shift = 1 << m;
	int i = 0;
#pragma omp parallel for private(i)
	for (i = 0; i < (int)data.size(); i++) {
		if (checkBit(i % (1 << (pos1 + 1)), pos1) && checkBit(i % (1 << (pos2 + 1)), pos2)  && !checkBit(i % (1 << (m + 1)), m)) {
			swap(data[i], data[i + shift]);
		}
	}
}

void Quantum::PhaseShift(int x, double alp) {
	alp *= (1 - noize);
	cd mult(cos(alp), sin(alp));
	int shift = 1 << x;
	int i = 0;
#pragma omp parallel for private(i)
	for (i = 0; i < (int)data.size(); i++) {
		if (!checkBit(i % (1 << (x + 1)), x)) {
			data[i + shift] *= mult;
		}
	}
}

void Quantum::CPhaseShift(int x, int y, double alp){
	alp *= (1 - noize);
	cd mult(cos(alp), sin(alp));
	int shift = 1 << y;
	int i = 0;
#pragma omp parallel for private(i)
	for (i = 0; i < (int)data.size(); i++) {
		if (checkBit(i % (1 << (x + 1)), x) && !checkBit(i % (1 << (y + 1)), y)) {
			data[i + shift] *= mult;
		}
	}
}

void Quantum::CCPhaseShift(int x, int y, int z, double alp) {
	alp *= (1 - noize);
	cd mult(cos(alp), sin(alp));
	int shift = 1 << z;
	int i = 0;
#pragma omp parallel for private(i)
	for (i = 0; i < (int)data.size(); i++) {
		if (checkBit(i % (1 << (x + 1)), x) && checkBit(i % (1 << (y + 1)), y) && !checkBit(i % (1 << (z + 1)), z)) {
			data[i + shift] *= mult;
		}
	}
}

void Quantum::QFFT(int l, int r){
	for (int i = r - 1; i >= l; i--) {
		Hadamar(i);
		for (int j = i - 1; j >= l; j--) {
			CPhaseShift(j, i, PI / (1 << (i - j)));
		}
	}
}

void Quantum::ReverseQFFT(int l, int r){
	for (int i = l; i < r; i++) {
		for (int j = l; j < i; j++) {
			CPhaseShift(j, i, -PI / (1 << (i - j)));
		}
		Hadamar(i);
	}
}

void Quantum::ADD(int a, int l, int r, int znak) {
	for (int j = l; j < r; j++) {
		for (int i = j; i < r; i++) {
			if (checkBit(a, j - l)) {
				PhaseShift(i, PI / (1 << (i - j)) * znak);
			}
		}
	}
}

void Quantum::CADD(int a, int l, int r, int m, int znak) {
	for (int j = l; j < r; j++) {
		for (int i = j; i < r; i++) {
			if (checkBit(a, j - l)) {
				CPhaseShift(m, i, PI / (1 << (i - j)) * znak);
			}
		}
	}
}

void Quantum::CCADD(int a, int l, int r, int m1, int m2, int znak) {
	for (int j = l; j < r; j++) {
		for (int i = j; i < r; i++) {
			if (checkBit(a, j - l)) {
				CCPhaseShift(m1, m2, i, PI / (1 << (i - j)) * znak);
			}
		}
	}
}

void Quantum::CCADD_modN(int a, int l, int r, int m1, int m2){
	CCADD(a, l, r - 1, m1, m2, 1);
	ADD(N, l, r - 1, -1);

	ReverseQFFT(l, r - 1);
	CNOT(r - 2, r - 1);
	QFFT(l, r - 1);

	CADD(N, l, r - 1, r - 1, 1);
	CCADD(a, l, r - 1, m1, m2, -1);

	ReverseQFFT(l, r - 1);
	NOT(r - 2);
	CNOT(r - 2, r - 1);
	NOT(r - 2);
	QFFT(l, r - 1);
	CCADD(a, l, r - 1, m1, m2, 1);
}

void Quantum::CMULT(int a, int l, int r, int m) {
	for (int i = l; i < (l + r) / 2; i++) {
		CCADD_modN(a, l + (r - l) / 2, r, i, m);
		(a *= 2) %= N;
	}
}

void Quantum::ReverseCMULT(int a, int l, int r, int m) {
	for (int i = l; i < (l + r) / 2; i++) {
		CCADD_modN(N - a, l + (r - l) / 2, r, i, m);
		(a *= 2) %= N;
	}
}

void Quantum::CSWAP(int l, int r, int m) {
	for (int i = l; i < l + (r - l) / 2; i++) {
		CNOT(i + (r - l) / 2, i);
		CCNOT(m, i, i + (r - l) / 2);
		CNOT(i + (r - l) / 2, i);
	}
}

void Quantum::UnitCMULT(int a, int l, int r, int m){
	CMULT(a, l, r, m);
	ReverseQFFT(l + (r - l) / 2, r - 1);
	CSWAP(l, r - 1, m);
	QFFT(l + (r - l) / 2, r - 1);
	ReverseCMULT(inverse_element(a, N), l, r, m);
}

void Quantum::Shor(int l, int r){
	data[1 << (logSize * 2)] = { 1, 0 };
	for (int i = 0; i < 2 * logSize; i++) {
		Hadamar(i);
	}
	QFFT(logSize * 3 + 1, r - 1);
	int _a = a;
	for (int i = 0; i < 2 * logSize; i++) {
		UnitCMULT(_a, logSize * 2, r, logSize * 2 - i - 1);
		(_a *= _a) %= N;
	}
	ReverseQFFT(logSize * 3 + 1, r - 1);
	ReverseQFFT(l, logSize * 2);
}

Quantum::Quantum(int a, int N, double noize):
	a(a), N(N), noize(noize),
	logSize(ceil(log2(N))) {
	data.resize(1 << (4 * logSize + 3));
}

vector<cd> Quantum::getData(){
	return data;
}

int Quantum::getQubitNumber() {
	return logSize;
}

int Quantum::measure_qubit(){
	double r = generate_probability();
	double sum = 0;
	int index = -1;
	while (sum <= r) {
		sum += norm(data[++index]);
	}
	return index;
}

void Quantum::measures(int times){
	counter.assign(1 << (logSize * 2), 0);
	for (int i = 0; i < times; i++) {
		counter[measure_qubit() % (1 << (logSize * 2))]++;
	}
}

vector<int> Quantum::getCounter() {
	return counter;
}

int Quantum::getPeriod(int measuresTime){
	Shor(0, 4 * logSize + 3);
	measures(measuresTime);
	int ans = 0;
	for (int i = 0; i < (int)counter.size(); i++) if (counter[i]) {
		ans = gcd(ans, i);
	}
	return (int)counter.size() / ans;
}