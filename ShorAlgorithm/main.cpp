#include<iostream>
#include<iomanip>
#include"solve.h"
#include<omp.h>
#pragma comment(linker, "/STACK:5000000000")

#define deb(a) cout << #a << " = " << a << '\n';

static const double eps = 1e-9;

using namespace std;

int main() {
	//ios_base::sync_with_stdio(0); cin.tie(0); cout.tie(0);
	int x;
	cin >> x;
	double qwe = omp_get_wtime();
	solve(x);
	qwe = omp_get_wtime() - qwe;
	cout << fixed << setprecision(20) << "TIME ELAPSED: " << qwe << '\n';
	return 0;
}