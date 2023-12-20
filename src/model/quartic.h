#ifndef TOSN_QUARTIC_H
#define TOSN_QUARTIC_H

#include <iostream>
#include <complex>
#include <iostream>

using namespace std;

// The solve_quartic routine solves the generic quartic equation:
//
//     a * x^4 + b * x^3 + c * x^2 + d * x + e == 0
//
// Usage:
//
//     solve_quartic({a, b, c, d, e}, roots).

void solve_quartic(complex<double> coefficients[5], complex<double> roots[4]);

#endif // TOSN_QUARTIC_H
