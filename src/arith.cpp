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

    AlgebraicComplexNumber<DenseNumberStore> num (coefs.size(), scale);
    for (s64 idx = 0; idx < coefs.size(); idx++) {
        num.coefficients.set_fp(idx, coefs[idx]);
    }

    return num;
}

std::ostream& AlgebraicComplexNumbers::operator<<(std::ostream& out, const AlgebraicComplexNumber<DenseNumberStore>& val) {
    for (s64 i = 0; i < val.coefficients.width; i++) {
        s64 coef_val = mpz_get_si(val.coefficients.numbers[i]);
        out << "(" << coef_val << "e^" << i << "/" << val.coefficients.width << ")";
        if (i + 1 < val.coefficients.width) {
            out << " + ";
        }
    }

    out << " sf=" << mpz_get_si(val.scaling_factor);

    return out;
}

