#include "quartic.h"

static complex<double> complex_sqrt(const complex<double> &z) {
    return pow(z, 1. / 2.);
}

static complex<double> complex_cbrt(const complex<double> &z) {
    return pow(z, 1. / 3.);
}

void solve_quartic(complex<double> coefficients[5], complex<double> roots[4]) {
    const int size = 5;
    std::complex<double> temp;

    for (int i = 0; i < size / 2; ++i) {
        temp = coefficients[i];
        coefficients[i] = coefficients[size - 1 - i];
        coefficients[size - 1 - i] = temp;
    }

    // The algorithm below was derived by solving the quartic in Mathematica, and simplifying the resulting expression by hand.

    const complex<double> a = coefficients[4];
    const complex<double> b = coefficients[3] / a;
    const complex<double> c = coefficients[2] / a;
    const complex<double> d = coefficients[1] / a;
    const complex<double> e = coefficients[0] / a;

    const complex<double> Q1 = c * c - 3. * b * d + 12. * e;
    const complex<double> Q2 = 2. * c * c * c - 9. * b * c * d + 27. * d * d + 27. * b * b * e - 72. * c * e;
    const complex<double> Q3 = 8. * b * c - 16. * d - 2. * b * b * b;
    const complex<double> Q4 = 3. * b * b - 8. * c;

    const complex<double> Q5 = complex_cbrt(Q2 / 2. + complex_sqrt(Q2 * Q2 / 4. - Q1 * Q1 * Q1));
    const complex<double> Q6 = (Q1 / Q5 + Q5) / 3.;
    const complex<double> Q7 = 2. * complex_sqrt(Q4 / 12. + Q6);

    roots[0] = (-b - Q7 - complex_sqrt(4. * Q4 / 6. - 4. * Q6 - Q3 / Q7)) / 4.;
    roots[1] = (-b - Q7 + complex_sqrt(4. * Q4 / 6. - 4. * Q6 - Q3 / Q7)) / 4.;
    roots[2] = (-b + Q7 - complex_sqrt(4. * Q4 / 6. - 4. * Q6 + Q3 / Q7)) / 4.;
    roots[3] = (-b + Q7 + complex_sqrt(4. * Q4 / 6. - 4. * Q6 + Q3 / Q7)) / 4.;
}