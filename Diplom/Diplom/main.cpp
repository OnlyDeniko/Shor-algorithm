#include<iostream>
#include<iomanip>
#include<fstream>
#include<omp.h>
#include"Quantum.h"
#include"Tmath.h"

using namespace std;
ofstream out;
void solve(int n, double noize) {
	if (is_prime_rabin_miller(n)) {
		out << "PRIME\n";
		return;
	}
	int log = ceil(log2(n));
	for (int i = 2; i <= log; i++) {
		int l = 2, r = ceil(pow(n, 1. / i)) + 1;

		while (r - l > 1) {
			int mid = (l + r) >> 1ll;

			long long res = bin_pow(mid, i, LLONG_MAX);
			if (res == -1 || res > n) {
				r = mid - 1;
			}
			else l = mid;
		}
		if (bin_pow(l, i, LLONG_MAX) == n) {
			out << l << '\n';
			// cout << l << '^' << i << '\n';
			return;
		}
		if (bin_pow(r, i, LLONG_MAX) == n) {
			out << r << '\n';
			// cout << r << '^' << i << '\n';
			return;
		}
	}
	int step = 0;
	for (int x = 2; x < n; x++) {
		if (gcd(x, n) != 1) continue;
		step++;
		Quantum magic(x, n, noize);
		int R = magic.getPeriod();
		int logSize = magic.getQubitNumber();
		if (R % 2 == 0) {
			int bp = bin_pow(x, R / 2, n);
			int first = gcd(bp + 1, n);
			int second = gcd(bp - 1, n);
			bool nice = false;
			if (n % first == 0 && first != n && first != 1) {
				nice = true;
				out << first << '\n';
				return;
			}
			if (n % second == 0 && second != n && second != 1) {
				nice = true;
				out << second << '\n';
				return;
			}
		}
	}
}

int main() {
	/*double t1 = omp_get_wtime();
	Quantum a(2, 33, 0.0);
	int R = a.getPeriod();
	t1 = omp_get_wtime() - t1;
	cout << R << ' ' << t1 << '\n';*/
	
	out.open("output.txt");
	for (int n = 2; n < 64; n++) {
		out << n << ' ';
		cout << n << '\n';
		double t1 = omp_get_wtime();
		solve(n, 0.0);
		t1 = omp_get_wtime() - t1;
		out << fixed << setprecision(8) << t1 << " sec\n";
	}
	out.close();
	/*ofstream out;
	out.open("output.txt");
	vector<int> a = { 2 };
	vector<double> noizes = { 0, 0.01, 0.05, 0.1 };
	vector<int> measurements = { 5000, 10000, 50000, 100000 };
	for(int i = 0; i < noizes.size();i++) {
		cout << noizes[i] << '\n';
		for (int j : a) {
			Quantum s(j, 15, noizes[i]);
			s.getPeriod(measurements[i]);
			out << fixed << setprecision(20) << j << ' ' << 15 << ' ' << noizes[i] << ' ' << (1 << (s.getQubitNumber() * 2)) << '\n';
			auto counter = s.getCounter();
			int cnt = 0;
			for (int k = 0; k < counter.size(); k++) if (counter[k]) cnt++;
			out << cnt << '\n';
			for (int k = 0; k < counter.size(); k++) if (counter[k]) out << k << ' ' << counter[k] << '\n';
		}
	}
	out.close();*/
}