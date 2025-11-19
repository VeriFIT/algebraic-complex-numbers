#include <cassert>
#include <cmath>
#include <cstdlib>
#include <ostream>
#include <iostream>

#include "arith.hpp"

using namespace AlgebraicComplexNumbers;

std::ostream& operator<<(std::ostream& os, const ComplexNumber& number) {
    const char* sign = number.im > 0 ? "+" : "-";

    os << number.real << " "
       << sign
       << std::abs(number.im) << " i";

    return os;
}

std::ostream& AlgebraicComplexNumbers::operator<<(std::ostream& os, const AlgebraicComplexNumber4& number) {
    os << "("
       << mpz_get_si(number.a) << ", "
       << mpz_get_si(number.b) << ", "
       << mpz_get_si(number.c) << ", "
       << mpz_get_si(number.d) << ", "
       << mpz_get_si(number.k) << ")";
    return os;
}

std::ostream& AlgebraicComplexNumbers::operator<<(std::ostream& os, const FixedPrecisionComplexNumber& number) {
    os << "("
       << number.a << ", "
       << number.b << ", "
       << number.c << ", "
       << number.d << ", "
       << number.k << ")";
    return os;
}

AlgebraicComplexNumber<DenseNumberStore> AlgebraicComplexNumbers::from_fp_vector(std::vector<s64> coefs, s64 scale) {
    assert (std::popcount(coefs.size()) == 1);

    AlgebraicComplexNumber<DenseNumberStore> num (coefs.size());
    for (s64 idx = 0; idx < coefs.size(); idx++) {
        num.coefficients.set_fp(idx, coefs[idx]);
    }

    return num;
}

