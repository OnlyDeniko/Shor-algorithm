#include "solve.h"
#include<set>
#include"algoNumbers.h"
#include"Quantum.h"

void solve(int n, double noize) {
	cout << "Noize = " << fixed << setprecision(12) << noize << '\n';
	if (!(n & 1ll)) {
		cout << 2 << '\n';
		return;
	}

	if (is_prime_rabin_miller(n)) {
		cout << "PRIME\n";
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
			cout << "DEGREE ";
			cout << l << '^' << i << '\n';
			return;
		}
		if (bin_pow(r, i, LLONG_MAX) == n) {
			cout << "DEGREE ";
			cout << r << '^' << i << '\n';
			return;
		}
	}
	int step = 0;
	for (int i = 2; i < n; i++) {
		int x = i;
		if (gcd(x, n) != 1) continue;
		step++;
		Quantum magic(x, n, noize);
		int R = magic.Shor();
		cout << x << ' ' << R / 2 << ' ' << n << '\n';
		if (!(R & 1ll)) {
			int bp = bin_pow(x, R / 2, n);
			int first = gcd(bp + 1, n);
			int second = gcd(bp - 1, n);
			if (n % first == 0 && first != n && first != 1) {
				cout << "Found at step " << step << ", x = " << x << '\n';
				cout << first << '\n';
			}
			else if (n % second == 0 && second != n && second != 1) {
				cout << "Found at step " << step << ", x = " << x << '\n';
				cout << second << '\n';
			}
		}
	}
	cout << "Not found\n";
}
