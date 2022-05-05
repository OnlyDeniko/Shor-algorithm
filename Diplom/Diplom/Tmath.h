#pragma once

long long bin_mul(long long a, long long b, long long mod);
long long bin_pow(long long a, long long b, long long mod);

long long gcd(long long a, long long b);
long long gcdex(long long a, long long b, long long& x, long long& y);
long long inverse_element(long long a, long long m);

int generate_random_number(int l, int r);
double generate_probability();

bool is_prime_sqrt_check(int x);
bool is_prime_euler_algo(int x);
bool is_prime_rabin_miller(int x);
