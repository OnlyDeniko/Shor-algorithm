#include<iostream>
#include<vector>
#include<iomanip>
#include<omp.h>
#include"solve.h"

#define deb(a) cout << #a << " = " << a << '\n';

static const double eps = 1e-9;

using namespace std;

int main() {
	int x;
	cin >> x;
	vector<double> noizes = { 0 };
	for (double noize : noizes) {
		double t1 = omp_get_wtime();
		solve(x, noize);
		t1 = omp_get_wtime() - t1;
		cout << fixed << setprecision(20) << "TIME ELAPSED: " << t1 << '\n';
	}
	return 0;
}