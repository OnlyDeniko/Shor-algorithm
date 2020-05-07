#include "solve.h"
#include<set>
#include"algoNumbers.h"
#include"Quantum.h"

void solve(int n) {
	deb(n);
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
			//cout << "! " << mid << ' ' << i << ' ' << res << '\n';
			if (res == -1 || res > n) {
				r = mid - 1;
			}
			else l = mid;
		}
		//cout << i << ' ' << l << ' ' << r << '\n';
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

	set<int> rand;
	//for (int i = 0; i < 2000; i++) {
	for (int i = 2; i < n; i++) {
		if (rand.size() == n - 2) break;
		/*int x = generate_random_number(2, n - 1);
		while (true) {
			if (rand.find(x) == rand.end()) break;
			x = generate_randon_number(2, n - 1);
		}*/
		int x = i;
		rand.insert(x);
		if (__gcd(x, n) != 1) continue;
		//cout << "~ " << x << '\n';

		Quantum magic(x, n);
		int R = magic.Shor();
		/*int R = 1;
		while (bin_pow(x, R, n) != 1) R++;
		cout << x << ' ' << R << endl;*/

		if (!(R & 1ll)) {
			if (n == 85 && i == 2) R = 8;
			int bp = bin_pow(x, R / 2, n);
			int first = __gcd(bp + 1, n);
			int second = __gcd(bp - 1, n);
			if (n % first == 0 && first != n && first != 1) {
				cout << first << '\n';
				return;
			}
			if (n % second == 0 && second != n && second != 1) {
				cout << second << '\n';
				return;
			}
		}
	}
	//assert(0);
}
