#pragma once
#include<chrono>
#include<random>
#include<ctime>

int bin_mul(long long a, long long b, long long mod);
int bin_pow(long long a, long long b, long long mod);
int generate_random_number(int l, int r);

int __gcd(int a, int b);
int gcdex(int a, int b, int &x, int &y);

bool is_prime_sqrt_check(int x);
bool is_prime_euler_algo(int x);
bool is_prime_rabin_miller(int x);
