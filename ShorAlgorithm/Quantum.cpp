#include "Quantum.h"
#include<map>

int Quantum::reverseNumber(int x, int len) {
	string s;
	while (x) {
		s += char((x % 2) + '0');
		x /= 2;
	}
	reverse(s.begin(), s.end());
	while (s.size() < len) s = '0' + s;
	int ans = 0;
	int step2 = 1;
	for (auto i : s) {
		if (i == '1') {
			ans += step2;
		}
		step2 <<= 1;
	}
	return ans;
}

string Quantum::perevod(int x, int k) {
	string ans;
	while (x) {
		ans += char(x % k + '0');
		x /= k;
	}
	reverse(ans.begin(), ans.end());
	if (ans.size() == 0) ans = "0";
	return ans;
}

Quantum::Quantum(int _a, int _N) {
	a = _a;
	N = _N;
	int cnt_qub = (int)ceil(log2(N)) + 1;
	int sz = 1 << (2 * cnt_qub - 2);
	firstRegister.assign(sz, { 0, 0 });
	firstRegister[0] = { 1, 0 };
	secondRegister.assign(1 << (cnt_qub + 2), { 0, 0 });
	secondRegister[reverseNumber(1, cnt_qub + 1) << 1] = { 1, 0 };
	thirdRegister.assign(1 << (cnt_qub + 2), { 0, 0 });;
	thirdRegister[0] = { 1, 0 };
	fourthRegister = { 0, 0 };
}

vector<cd> Quantum::QFFT(const vector<cd>& input, int cnt_qubits) {
	auto a = input;
	for (int i = 0; i < cnt_qubits; i++) {
		a = Hadamar(a, i + 1);
		for (int j = i + 1; j < cnt_qubits; j++) {
			a = Controlled_Rm(a, i + 1, j + 1, j - i + 1);
		}
	}
	return a;
}

vector<cd> Quantum::ReverseQFFT(const vector<cd>& input, int cnt_qubits) {
	auto a = input;
	for (int i = cnt_qubits - 1; i >= 0; i--) {
		int gg = 0;
		for (int j = cnt_qubits - 1; j >= i + 1; j--) {
			a = ReverseControlled_Rm(a, i + 1, j + 1, j - i + 1);
		}
		a = Hadamar(a, i + 1);
	}
	return a;
}

vector<cd> Quantum::Hadamar(const vector<cd>& a, int ind) {
	int gg = (int)a.size();
	int step = 0;
	while (gg != 1) {
		step++;
		gg >>= 1;
	}

	//int index = step - ind;
	int index = ind - 1;
	int finish = 1 << (step - 1);
	vector<cd> matrix(a);

#pragma omp parallel
	{
#pragma omp for schedule(static)
		for (int mask = 0; mask < finish; mask++) {
			int ind1 = 0, ind2 = 0;
			for (int i = 0; i < index; i++) {
				ind1 |= mask & (1 << i);
			}
			for (int i = index + 1; i < step; i++) {
				bool gg = mask & (1 << (i - 1));
				ind1 |= (gg << i);
			}
			ind2 = ind1;
			ind2 |= (1 << index);
			cd q = a[ind1];
			cd w = a[ind2];
			matrix[ind1] = 1 / sqrt(2.) * (q + w);
			matrix[ind2] = 1 / sqrt(2.) * (q - w);
		}
	}
	return matrix;
}

vector<cd> Quantum::Controlled_Z(const vector<cd>& a, int x, int y) {
	int gg = (int)a.size();
	int step = 0;
	while (gg != 1) {
		step++;
		gg >>= 1;
	}
	/*x = step - x;
	y = step - y;*/
	x--;
	y--;
	int finish = 1 << (step - 2);
	vector<cd> res(a);
#pragma omp parallel
	{
#pragma omp for schedule(static)
		for (int mask = 0; mask < finish; mask++) {
			int ind1 = 0;
			for (int i = 0; i < x; i++) {
				ind1 |= mask & (1 << i);
			}
			for (int i = x + 1; i < y; i++) {
				bool gg = mask & (1 << (i - 1));
				ind1 |= (gg << i);
			}
			for (int i = y + 1; i < step; i++) {
				bool gg = mask & (1 << (i - 2));
				ind1 |= (gg << i);
			}
			int ind2 = ind1;
			ind2 |= (1 << y);
			int ind3 = ind1;
			ind3 |= (1 << x);
			int ind4 = ind3;
			ind4 |= (1 << y);
			res[ind1] = a[ind1];
			res[ind2] = a[ind2];
			res[ind3] = a[ind3];
			res[ind4] = -a[ind4];

		}
	}
	return res;
}

vector<cd> Quantum::CNOT(const vector<cd>& a, int x, int y) {
	int gg = (int)a.size();
	int step = 0;
	while (gg != 1) {
		step++;
		gg >>= 1;
	}
	/*x = step - x;
	y = step - y;*/
	x--;
	y--;
	int finish = 1 << (step - 2);
	vector<cd> res(a);
#pragma omp parallel
	{
#pragma omp for schedule(static)
		for (int mask = 0; mask < finish; mask++) {
			int ind1 = 0;
			for (int i = 0; i < x; i++) {
				ind1 |= mask & (1 << i);
			}
			for (int i = x + 1; i < y; i++) {
				bool gg = mask & (1 << (i - 1));
				ind1 |= (gg << i);
			}
			for (int i = y + 1; i < step; i++) {
				bool gg = mask & (1 << (i - 2));
				ind1 |= (gg << i);
			}
			int ind2 = ind1;
			ind2 |= (1 << y);
			int ind3 = ind1;
			ind3 |= (1 << x);
			int ind4 = ind3;
			ind4 |= (1 << y);
			res[ind1] = a[ind1];
			res[ind2] = a[ind2];
			res[ind4] = a[ind3];
			res[ind3] = a[ind4];
		}
	}
	return res;
}

vector<cd> Quantum::PhaseShift(const vector<cd>& a, int ind, int m) {
	auto ans = a;
	double alp = (double)2 * acos(-1.) / (1 << m);
	cd mult(cos(alp), sin(alp));
	int gg = (int)a.size();
	int step = 0;
	while (gg != 1) {
		step++;
		gg >>= 1;
	}
	int index = ind - 1;
	int finish = 1 << (step - 1);
	vector<cd> matrix(a);

#pragma omp parallel
	{
#pragma omp for schedule(static)
		for (int mask = 0; mask < finish; mask++) {
			int ind1 = 0, ind2 = 0;
			for (int i = 0; i < index; i++) {
				ind1 |= mask & (1 << i);
			}
			for (int i = index + 1; i < step; i++) {
				bool gg = mask & (1 << (i - 1));
				ind1 |= (gg << i);
			}
			ind2 = ind1;
			ind2 |= (1 << index);
			matrix[ind1] = a[ind1];
			matrix[ind2] = a[ind2] * mult;
		}
	}
	return matrix;
}

vector<cd> Quantum::ReversePhaseShift(const vector<cd>& a, int ind, int m) {
	auto ans = a;
	double alp = -(double)2 * acos(-1.) / (1 << m);
	cd mult(cos(alp), sin(alp));
	int gg = (int)a.size();
	int step = 0;
	while (gg != 1) {
		step++;
		gg >>= 1;
	}
	int index = ind - 1;
	int finish = 1 << (step - 1);
	vector<cd> matrix(a);

#pragma omp parallel
	{
#pragma omp for schedule(static)
		for (int mask = 0; mask < finish; mask++) {
			int ind1 = 0, ind2 = 0;
			for (int i = 0; i < index; i++) {
				ind1 |= mask & (1 << i);
			}
			for (int i = index + 1; i < step; i++) {
				bool gg = mask & (1 << (i - 1));
				ind1 |= (gg << i);
			}
			ind2 = ind1;
			ind2 |= (1 << index);
			matrix[ind1] = a[ind1];
			matrix[ind2] = a[ind2] * mult;
		}
	}
	return matrix;
}

vector<cd> Quantum::Controlled_Rm(const vector<cd>& a, int x, int y, int m) {
	double alp = (double)2 * acos(-1.) / (1 << m);
	cd mult(cos(alp), sin(alp));
	int gg = (int)a.size();
	int step = 0;
	while (gg != 1) {
		step++;
		gg >>= 1;
	}
	x--;
	y--;
	int finish = 1 << (step - 2);
	vector<cd> res(a.size());
#pragma omp parallel
	{
#pragma omp for
		for (int mask = 0; mask < finish; mask++) {
			int ind1 = 0;
			for (int i = 0; i < x; i++) {
				ind1 |= mask & (1 << i);
			}
			for (int i = x + 1; i < y; i++) {
				bool gg = mask & (1 << (i - 1));
				ind1 |= (gg << i);
			}
			for (int i = y + 1; i < step; i++) {
				bool gg = mask & (1 << (i - 2));
				ind1 |= (gg << i);
			}
			int ind2 = ind1;
			ind2 |= (1 << y);
			int ind3 = ind1;
			ind3 |= (1 << x);
			int ind4 = ind3;
			ind4 |= (1 << y);
			res[ind1] = a[ind1];
			res[ind2] = a[ind2];
			res[ind3] = a[ind3];
			res[ind4] = a[ind4] * mult;
		}
	}
	return res;
}

vector<cd> Quantum::ReverseControlled_Rm(const vector<cd>& a, int x, int y, int m) {
	double alp = -(double)2 * acos(-1.) / (1 << m);
	cd mult(cos(alp), sin(alp));
	int gg = (int)a.size();
	int step = 0;
	while (gg != 1) {
		step++;
		gg >>= 1;
	}
	x--;
	y--;
	int finish = 1 << (step - 2);
	vector<cd> res(a.size());
#pragma omp parallel
	{
#pragma omp for
		for (int mask = 0; mask < finish; mask++) {
			int ind1 = 0;
			for (int i = 0; i < x; i++) {
				ind1 |= mask & (1 << i);
			}
			for (int i = x + 1; i < y; i++) {
				bool gg = mask & (1 << (i - 1));
				ind1 |= (gg << i);
			}
			for (int i = y + 1; i < step; i++) {
				bool gg = mask & (1 << (i - 2));
				ind1 |= (gg << i);
			}
			int ind2 = ind1;
			ind2 |= (1 << y);
			int ind3 = ind1;
			ind3 |= (1 << x);
			int ind4 = ind3;
			ind4 |= (1 << y);
			res[ind1] = a[ind1];
			res[ind2] = a[ind2];
			res[ind3] = a[ind3];
			res[ind4] = a[ind4] * mult;
		}
	}
	return res;
}

vector<cd> Quantum::FADD(vector<cd> b, int a) {
	int gg = (int)b.size();
	int step = 0;
	while (gg != 1) {
		step++;
		gg >>= 1;
	}
	string to_2 = perevod(a, 2);
	while (to_2.size() < step) to_2 = '0' + to_2;

	for (int i = 0; i < step; i++) {
		for (int j = 1; j <= step - i; j++) {
			if (to_2[step - i - j] == '0') continue;
			b = PhaseShift(b, i + 1, j);
		}
	}
	return b;
}

vector<cd> Quantum::ReverseFADD(vector<cd> b, int a) {
	int gg = (int)b.size();
	int step = 0;
	while (gg != 1) {
		step++;
		gg >>= 1;
	}
	string to_2 = perevod(a, 2);
	while (to_2.size() < step) to_2 = '0' + to_2;

	for (int i = 0; i < step; i++) {
		for (int j = 1; j <= step - i; j++) {
			if (to_2[step - i - j] == '0') continue;
			b = ReversePhaseShift(b, i + 1, j);
		}
	}
	return b;
}

vector<cd> Quantum::Fredkin(const vector<cd>& a, int x, int y, int z) {
	int gg = (int)a.size();
	int step = 0;
	while (gg != 1) {
		step++;
		gg >>= 1;
	}
	/*x = step - x;
	y = step - y;*/
	x--;
	y--;
	z--;
	int finish = 1 << (step - 3);
	vector<cd> res(a);
	#pragma omp parallel
	{
		#pragma omp for schedule(static)
		for (int mask = 0; mask < finish; mask++) {
			int ind1 = 0;
			for (int i = 0; i < x; i++) {
				ind1 |= mask & (1 << i);
			}
			for (int i = x + 1; i < y; i++) {
				bool gg = mask & (1 << (i - 1));
				ind1 |= (gg << i);
			}
			for (int i = y + 1; i < z; i++) {
				bool gg = mask & (1 << (i - 2));
				ind1 |= (gg << i);
			}
			for (int i = z + 1; i < step; i++) {
				bool gg = mask & (1 << (i - 3));
				ind1 |= (gg << i);
			}
			int ind2(ind1), ind3(ind1), ind4(ind1), ind5(ind1), ind6(ind1), ind7(ind1), ind8(ind1), ind9(ind1);
			ind2 |= (1 << z);
			ind3 |= (1 << y);
			ind4 |= (1 << z); ind4 |= (1 << y);
			ind5 |= (1 << x);
			ind6 |= (1 << x);                 ind6 |= (1 << z);
			ind7 |= (1 << x); ind7 |= (1 << y);
			ind8 |= (1 << x); ind8 |= (1 << y); ind8 |= (1 << z);
			res[ind1] = a[ind1];
			res[ind2] = a[ind2];
			res[ind3] = a[ind3];
			res[ind4] = a[ind4];
			res[ind5] = a[ind5];
			res[ind6] = a[ind7];
			res[ind7] = a[ind6];
			res[ind8] = a[ind8];
		}
	}
	return res;
}

void Quantum::FADD_modN(cd & _1, cd & _2, vector<cd>& b, cd & _4, int a, int N, int cnt_qub) {
	if (norm(_1) > 0 && norm(_2) > 0) {
		b = FADD(b, a);
	}

	b = ReverseFADD(b, N);
	b = ReverseQFFT(b, cnt_qub);

	int ind = -1;
	for (int i = 0; i < b.size(); i++) {
		if (fabs(norm(b[i]) - 1) < eps) {
			ind = i;
			break;
		}
	}
	string tmp = perevod(ind, 2);
	if (tmp.back() == '0') {}
	else {
		_4 = { 1, 0 };
	}

	b = QFFT(b, cnt_qub);

	if (norm(_4) > 0) {
		b = FADD(b, N);
	}
	if (norm(_1) > 0 && norm(_2) > 0) {
		b = ReverseFADD(b, a);
	}
	b = ReverseQFFT(b, cnt_qub);
	ind = -1;
	for (int i = 0; i < b.size(); i++) {
		if (fabs(norm(b[i]) - 1) < eps) {
			ind = i;
			break;
		}
	}
	tmp = perevod(ind, 2);
	if (tmp.back() == '0')
		tmp.back() = '1';
	else tmp.back() = '0';

	if (tmp.back() == '1') {
		_4 = { 0, 0 };
	}
	b = QFFT(b, cnt_qub);
	if (norm(_1) > 0 && norm(_2) > 0) {
		b = FADD(b, a);
	}
}

void Quantum::ReverseFADD_modN(cd & _1, cd & _2, vector<cd>& b, cd & _4, int a, int N, int cnt_qub) {
	if (norm(_1) > 0 && norm(_2) > 0) {
		b = ReverseFADD(b, a);
	}

	b = ReverseQFFT(b, cnt_qub);

	int ind = -1;
	for (int i = 0; i < b.size(); i++) {
		if (fabs(norm(b[i]) - 1) < eps) {
			ind = i;
			break;
		}
	}
	//cout << "! " << _4 << '\n';
	if (ind == -1) throw 1;
	string tmp = perevod(ind, 2);
	if (tmp.back() == '1') _4 = { 1, 0 };

	//cout << "! " << _4 << '\n';
	b = QFFT(b, cnt_qub);

	if (norm(_4) > 0) {
		b = FADD(b, N);
	}
	if (norm(_1) > 0 && norm(_2) > 0) {
		b = FADD(b, a);
	}
	b = ReverseQFFT(b, cnt_qub);
	ind = -1;
	for (int i = 0; i < b.size(); i++) {
		if (fabs(norm(b[i]) - 1) < eps) {
			ind = i;
			break;
		}
	}
	if (ind == -1) throw 1;
	//cout << "! " << _4 << '\n';
	tmp = perevod(ind, 2);
	if (tmp.back() == '0') tmp.back() = '1';
	else tmp.back() = '0';
	if (tmp.back() == '1') {
		_4 = { 0, 0 };
	}
	//cout << "! " << _4 << '\n';
	b = QFFT(b, cnt_qub);
	if (norm(_1) > 0 && norm(_2) > 0) {
		b = ReverseFADD(b, a);
	}
}

void Quantum::CMULT_modN(cd & _1, vector<cd>& x, vector<cd>& b, cd & _4, int a, int N, int cnt_qub) {
	b = QFFT(b, cnt_qub);
	int step2 = 1;
	cd qwe = { 1, 0 };
	for (int i = 0; i < log2(x.size()); i++) {
		int gg = (step2 * a) % N;
		bool ok = 0;
		for (int j = 0; j < x.size(); j++) {
			string tmp = perevod(j, 2);
			while (tmp.size() < log2(x.size())) tmp = '0' + tmp;
			if (tmp[i] == '1' && norm(x[j]) > 0) {
				ok = 1;
				break;
			}
		}
		if (ok) FADD_modN(_1, qwe, b, _4, reverseNumber(gg, cnt_qub), reverseNumber(N, cnt_qub), cnt_qub);
		step2 *= 2;
		step2 %= N;
	}
	b = ReverseQFFT(b, cnt_qub);
}

void Quantum::ReverseCMULT_modN(cd & _1, vector<cd>& x, vector<cd>& b, cd & _4, int a, int N, int cnt_qub) {
	b = QFFT(b, cnt_qub);
	int step2 = 1;
	cd qwe = { 1, 0 };
	for (int i = 0; i < log2(x.size()); i++) {
		int gg = (step2 * a) % N;
		bool ok = 0;
		for (int j = 0; j < x.size(); j++) {
			string tmp = perevod(j, 2);
			while (tmp.size() < log2(x.size())) tmp = '0' + tmp;
			if (tmp[i] == '1' && norm(x[j]) > 0) {
				ok = 1;
				break;
			}
		}
		if (ok) ReverseFADD_modN(_1, qwe, b, _4, reverseNumber(gg, cnt_qub), reverseNumber(N, cnt_qub), cnt_qub);
		step2 *= 2;
		step2 %= N;
	}
	b = ReverseQFFT(b, cnt_qub);
}

void Quantum::UnitCMULT_modN(cd & _1, vector<cd>& x, vector<cd>& b, cd & _4, int a, int N, int cnt_qub) {
	if (norm(_1) > 0) {
		CMULT_modN(_1, x, b, _4, a, N, cnt_qub);
	}

	swap(x, b);// =)


	int ans, tmp;
	int gcd = gcdex(a, N, ans, tmp);

	try {
		if (gcd != 1) throw BadGCD();
		if (ans < 0) ans += N;
		if ((ans * a) % N != 1) throw BadAnswerReverse();
	}
	catch (Quantum::BadGCD) {
		cout << "BadGCD: UnitCMULT_modN\n";
		return;
	}
	catch (Quantum::BadAnswerReverse) {
		cout << "BadAnswerReverse: UnitCMULT_modN\n";
		return;
	}

	if (norm(_1) > 0) ReverseCMULT_modN(_1, x, b, _4, ans, N, cnt_qub);
}

int Quantum::Shor() {
	int cnt_qub = (int)ceil(log2(N)) + 1;
	int sz = 1 << (2 * cnt_qub - 2);
	vector<cd> kek(1 << (3 * cnt_qub - 3), 0);

	firstRegister = QFFT(firstRegister, log2(firstRegister.size()));


	for (int j = 0; j < firstRegister.size(); j++) {
		string tmp = perevod(j, 2);
		secondRegister.assign(1 << (cnt_qub + 2), { 0, 0 });
		secondRegister[reverseNumber(1, cnt_qub + 1) << 1] = { 1, 0 };
		thirdRegister.assign(1 << (cnt_qub + 2), { 0, 0 });
		thirdRegister[0] = { 1, 0 };
		while (tmp.size() < log2(firstRegister.size())) tmp = '0' + tmp;
		int gg = 1;
		long long step2 = 1;
		for (int k = 0; k < log2(firstRegister.size()); k++) {
			gg = bin_pow(a, step2, N);
			step2 *= 2;
			if (tmp[k] == '1') {
				cd qwe = { 1, 0 };
				UnitCMULT_modN(qwe, secondRegister, thirdRegister, fourthRegister, gg, N, cnt_qub + 2);
				thirdRegister.assign(1 << (cnt_qub + 2), { 0, 0 });
				thirdRegister[0] = { 1, 0 };
				#pragma omp parallel for
				for (int i = 0; i < secondRegister.size(); i++) {
					double q = secondRegister[i].real(), w = secondRegister[i].imag();
					if (fabs(secondRegister[i].real()) < eps) q = 0;
					if (fabs(secondRegister[i].imag()) < eps) w = 0;
					secondRegister[i] = { q, w };
				}
			}
		}

		int ans = -1;
		for (int k = 0; k < secondRegister.size(); k++) if (fabs(norm(secondRegister[k]) - 1) < eps) {
			ans = k;
			break;
		}
		//ans = reverseNumber(ans, cnt_qub + 2);
		ans >>= 3;
		kek[ans * (1 << (2 * cnt_qub - 2)) + j] = firstRegister[j];
	}
	for (int i = 0; i < cnt_qub - 1; i++) kek = measure_cubit(kek, log2(kek.size()));

	kek = QFFT(kek, log2(kek.size()));
	int ans = 0;

	for (int i = 0; i < kek.size(); i++) {
		if (norm(kek[i]) > eps) ans = __gcd(ans, i);
	}

	return kek.size() / ans;
}

vector<cd> Quantum::DenseCoding(const vector<cd>& a, int cnt_qubits) {
	auto res = a;
	int cnt_iter = (cnt_qubits - 2) >> 1;
	int gg = 3;
	for (int i = 0; i < cnt_iter; i++) {
		res = Hadamar(res, gg);
		res = CNOT(res, gg, gg + 1);
		res = CNOT(res, gg - 1, gg);
		res = Controlled_Z(res, gg - 2, gg);
		res = CNOT(res, gg, gg + 1);
		res = Hadamar(res, gg);
		gg += 2;
	}
	return res;
}

vector<cd> Quantum::OptimizedDenseCoding(const vector<cd>& a, int cnt_qubits) {
	auto res = a;
	int cnt_iter = (cnt_qubits - 2) >> 1;
	int gg = 2;
	for (int i = 0; i < cnt_iter; i++) {
		res = CNOT(res, gg, gg + 2);
		res = CNOT(res, gg - 1, gg + 1);
		gg += 2;
	}
	return res;
}

void Quantum::TestDenseCoding(const vector<complex<double>>& input, int cnt_qubits, int cnt_tests) {
	for (int test = 0; test < cnt_tests; test++) {
		cout << "ITERATION : " << test + 1 << '\n';

		double t1, t2;
		t1 = omp_get_wtime();

		vector<complex<double> > a(input);
		a = DenseCoding(a, cnt_qubits);
		/*for (int i = 0; i < (1 << cnt_qubits); i++) {
			cout << bitset<4>(i) << " : " << a[i] << '\n';
		}*/
		/*test_for_one_calc(100, a, cnt_qub - 1);
		test_for_one_calc(100, a, cnt_qub);*/

		t2 = omp_get_wtime();
		t2 -= t1;
		cout << "TIME ELAPSED = " << fixed << setprecision(20) << t2 << '\n';
	}
}

void Quantum::TestOptimizedDenseCoding(const vector<complex<double>>& input, int cnt_qubits, int cnt_tests) {
	for (int test = 0; test < cnt_tests; test++) {
		cout << "ITERATION : " << test + 1 << '\n';

		double t1, t2;
		t1 = omp_get_wtime();

		vector<complex<double> > a(input);
		a = OptimizedDenseCoding(a, cnt_qubits);
		//for (auto i : a) cout << i << '\n';
		/*test_for_one_calc(100, a, cnt_qub - 1);
		test_for_one_calc(100, a, cnt_qub);*/

		t2 = omp_get_wtime();
		t2 -= t1;
		cout << "TIME ELAPSED = " << fixed << setprecision(20) << t2 << '\n';
	}
}

int Quantum::calc_whole_cubit(const vector<cd>& a) {
	vector<double> tmp;
	for (auto i : a) tmp.push_back(norm(i));

	int gg = generate_random_number(2, INT_MAX);
	gg = abs(gg);
	double p = (double)gg / INT_MAX;

	double sum = 0;
	int index = -1;
	for (int i = 0; i < (int)a.size(); i++) {
		if (sum + tmp[i] >= p) {
			index = i;
			break;
		}
		sum += tmp[i];
	}
	return index;
}

void Quantum::test_for_whole_calc(const int & tests, const vector<complex<double>>& a) {
	map<int, int> mapik;

	for (int i = 0; i < tests; i++) {
		mapik[calc_whole_cubit(a)]++;
	}
	for (auto i : mapik) {
		cout << i.first << " : ";
		double pr = (double)i.second / tests;
		cout.width(20);
		cout << "PROBABILITY = " << setprecision(10) << pr << '\n';
	}
}

vector<cd> Quantum::calc_one_cubit(const vector<cd>& a, int index) {
	int gg = (int)a.size();
	int step = 0;
	while (gg != 1) {
		step++;
		gg >>= 1;
	}
	index--;
	double p0(0), p1(0);
	for (int i = 0; i < (int)a.size(); i++) {
		bitset<128> bit(i);
		if (bit[index]) p1 += norm(a[i]);
		else p0 += norm(a[i]);
	}

	int g = generate_random_number(2, INT_MAX);
	g = abs(g);
	double p = (double)g / INT_MAX;
	vector<cd> ans(a.begin(), a.end());
	if (p > p0) {
		for (int i = 0; i < (int)a.size(); i++) {
			bitset<128> bit(i);
			if (!bit[index]) {
				ans[i] = 0;
			}
			else {
				ans[i] /= sqrt(p1);
			}
		}
	}
	else {
		for (int i = 0; i < (int)a.size(); i++) {
			bitset<128> bit(i);
			if (bit[index]) {
				ans[i] = 0;
			}
			else {
				ans[i] /= sqrt(p0);
			}
		}
	}
	return ans;
}

vector<cd> Quantum::measure_cubit(const vector<cd>& a, int index) {
	int gg = (int)a.size();
	int step = 0;
	while (gg != 1) {
		step++;
		gg >>= 1;
	}
	index--;
	double p0(0), p1(0);
	for (int i = 0; i < (int)a.size(); i++) {
		bitset<128> bit(i);
		if (bit[index]) p1 += norm(a[i]);
		else p0 += norm(a[i]);
	}
	//cout << p1 << ' ' << p0 << '\n';
	int g = generate_random_number(2, INT_MAX);
	g = abs(g);
	double p = (double)g / INT_MAX;
	vector<cd> ans(a.size() / 2);
	if (p > p0) {
		for (int i = 0; i < (int)ans.size(); i++) {
			ans[i] = a[i + ans.size()];
			ans[i] /= sqrt(p1);
		}
	}
	else {
		for (int i = 0; i < (int)ans.size(); i++) {
			ans[i] = a[i];
			ans[i] /= sqrt(p0);
		}
	}
	return ans;
}

void Quantum::test_for_one_calc(const int & tests, const vector<cd>& a, int index) {
	int gg = (int)a.size();
	int step = 0;
	while (gg != 1) {
		gg >>= 1;
		step++;
	}

	map<vector<cd>, int> mapik;
	for (int i = 0; i < tests; i++) {
		auto res = calc_one_cubit(a, index);
		mapik[res]++;
	}
	for (auto i : mapik) {

		cd res = 0;
		for (int j = 0; j < i.first.size(); j++) {
			bitset<30> tmp = j;
			cout << tmp << " : " << i.first[j] << '\n';
		}
		/*for (auto j : i.first) {
			cout << j << ' ';
			res += norm(j);
		}*/
		double pr = (double)i.second / tests;
		cout << "PROBABILITY = " << setprecision(10) << pr << '\n' << '\n';

	}
}

