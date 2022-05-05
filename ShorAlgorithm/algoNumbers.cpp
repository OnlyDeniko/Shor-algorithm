#include "algoNumbers.h"
#include<chrono>
#include<random>
#include<ctime>

int bin_mul(long long a, long long b, long long mod) {
	int res = 0;
	while (b) {
		if (b & 1) {
			(res += a) %= mod;
		}
		b >>= 1;
		(a += a) %= mod;
	}
	return res;
}

int bin_pow(long long a, long long b, long long mod) {
	int res = 1;
	while (b) {
		if (b & 1) {
			res = bin_mul(res, a, mod);
		}
		b >>= 1;
		a = bin_mul(a, a, mod);
	}
	return res;
}

int gcd(int a, int b) {
	while (b) {
		a %= b;
		std::swap(a, b);
	}
	return a;
}

int gcdex(int a, int b, int& x, int& y) {
	int c, q, xp, yp;
	x = yp = 1, y = xp = 0;
	while (b != 0) {
		q = a / b;
		c = a - q * b;
		a = b, b = c;
		c = x - q * xp;
		x = xp, xp = c;
		c = y - q * yp;
		y = yp, yp = c;
	}
	return a;
}

bool is_prime_sqrt_check(int x) {
	if (x <= 1) return 0;
	if (x <= 3) return 1;
	if (x % 2 == 0) return 0;
	for (int i = 2; i <= (int)sqrt(x); i++) {
		if (x % i == 0) {
			return 0;
		}
	}
	return 1;
}

int generate_random_number(int l, int r) {
	std::mt19937 gen(std::chrono::high_resolution_clock::now().time_since_epoch().count());
	int ans = gen();
	ans = abs(ans);
	ans = ans % (r - l + 1) + l;
	return ans;
}

bool is_prime_euler_algo(int x) {
	if (x <= 1) return 0;
	if (x <= 3) return 1;
	if (x % 2 == 0) return 0;
	for (int i = 0; i < 100; i++) {
		int rnd = generate_random_number(2, x - 1);
		if (bin_pow(rnd, x - 1, x) != 1) {
			return 0;
		}
	}
	return 1;
}

bool is_prime_rabin_miller(int x) {
	if (x <= 1) return 0;
	if (x <= 3) return 1;
	if (x % 2 == 0) return 0;
	int step = 0;
	int d = x - 1;
	while (d % 2 == 0) {
		d >>= 1;
		step++;
	}
	for (int i = 0; i < 100; i++) {
		int rnd = generate_random_number(2, x - 1);
		int a = bin_pow(rnd, d, x);
		if (a == 1 || a == x - 1) continue;
		int ok = 0;
		for (int j = 0; j < step - 1; j++) {
			a = bin_mul(a, a, x);
			if (a == 1) return 0;
			if (a == x - 1) {
				ok = 1;
				break;
			}
		}
		if (!ok) return 0;
	}
	return 1;
}
