#pragma once

#include "matrix.hpp"

using namespace AlgebraicComplexNumbers;

namespace AlgebraicComplexNumbers {

template <typename T>
ComplexMatrix<T> square_acn_matrix_from_ints(const std::vector<s64>& ints, T *zero) {
    u64 dim = static_cast<s64>(std::sqrt(ints.size()));
    assert(dim*dim == ints.size());

    ComplexMatrix<T> result(dim, dim, zero);

    for (u64 elem_idx = 0; elem_idx < ints.size(); elem_idx++) {
        s64 elem_int_value = ints[elem_idx];
        auto& target_elem = result.at_linear(elem_idx);
        mpz_set_si(target_elem.at(0), elem_int_value);
    }

    return result;
}

template <typename T>
ComplexMatrix<T> row_from_ints(const std::vector<s64>& row_data, T *zero) {
    ComplexMatrix<T> result(1, row_data.size(), zero);

    for (u64 elem_idx = 0; elem_idx < row_data.size(); elem_idx++) {
        s64 elem_int_value = row_data[elem_idx];
        auto& target_elem = result.at(0, elem_idx);
        mpz_set_si(target_elem.at(0), elem_int_value);
    }

    return result;
}

template <typename T>
ComplexMatrix<T> column_from_ints(const std::vector<s64>& column_data, T *zero) {
    ComplexMatrix<T> result(column_data.size(), 1, zero);

    for (u64 elem_idx = 0; elem_idx < column_data.size(); elem_idx++) {
        s64 elem_int_value = column_data[elem_idx];
        auto& target_elem = result.at(elem_idx, 0);
        mpz_set_si(target_elem.at(0), elem_int_value);
    }

    return result;
}


template <typename T>
struct MatrixBaker {
    T *zero;

    ComplexMatrix<T> square_acn_matrix_from_ints(const std::vector<s64>& ints) {
        return AlgebraicComplexNumbers::square_acn_matrix_from_ints<T>(ints, zero);
    }

    ComplexMatrix<T> row_from_ints(const std::vector<s64>& ints) {
        return AlgebraicComplexNumbers::row_from_ints<T>(ints, zero);
    }

    ComplexMatrix<T> column_from_ints(const std::vector<s64>& ints) {
        return AlgebraicComplexNumbers::column_from_ints<T>(ints, zero);
    }
};

}; // Namespace end
