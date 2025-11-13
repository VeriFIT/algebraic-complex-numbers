#include <cassert>
#include <cmath>
#include <cstdlib>
#include <new>
#include <ostream>
#include <iostream>
#include <vector>

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

s64 AlgebraicComplexNumbers::add_row_to_row_echelon_matrix_no_copy(ComplexMatrix& matrix, ComplexMatrix& row) {
    assert(matrix.width == row.width);

    if (row.contains_only_zeros()) return -1;

    u64 inserted_row_pivot_idx = row.find_nonzero_elem_in_row(0);

    s64 row_inserted_at = -1; // Not inserted

    for (u64 row_idx = 0; row_idx < matrix.height; row_idx++) {
        u64 row_pivot_idx = matrix.find_nonzero_elem_in_row(row_idx);

        if (row_pivot_idx < inserted_row_pivot_idx) {
            continue;
        }

        if (row_pivot_idx > inserted_row_pivot_idx) {
            matrix.insert_row_at(row, row_idx);
            row_inserted_at = row_idx;
            return row_inserted_at;
        }

        // Subtract the current matrix row and continue;
        auto& matrix_pivot = matrix.at(row_idx, row_pivot_idx);
        auto& row_pivot    = row.at(0, inserted_row_pivot_idx);

        row.subtract_from_ith_row(
            0,             // Which row of the `row` matrix to subtract from
            matrix_pivot,  // What coefficient should multiply the row before multiplication
            matrix,        // Matrix containing rows that we can subtract from
            row_idx,       // Which of the matrix rows are we subtracting
            row_pivot);    // What coeffcient should multiply the row

        inserted_row_pivot_idx = row.find_nonzero_elem_in_row(0);

        if (inserted_row_pivot_idx >= row.width) {
            return -1; // Not inserted, Gaussian elimination reduced the row to 0
        }
    }

    assert(false);  // Should be unreachable - the row is indepentent and we insert it, or our matrix has full rank and hence we do not insert it
}

s64 AlgebraicComplexNumbers::add_row_to_row_echelon_matrix(ComplexMatrix& matrix, const ComplexMatrix& row) {
    ComplexMatrix row_copy (1, row.width); // Make a local copy, since we will be modifying it
    for (auto row_elem_idx = 0; row_elem_idx < row.width; row_elem_idx++) {
        auto& elem_value = row.at(0, row_elem_idx);
        row_copy.set(0, row_elem_idx, elem_value);
    }
    return add_row_to_row_echelon_matrix_no_copy(matrix, row_copy);
}

u64 compute_square_matrix_dim_from_1d_repr(const std::vector<s64>& repr) {
    u64 dimension = static_cast<s64>(std::sqrt(repr.size()));

    assert(dimension*dimension == repr.size());

    return dimension;
}

ComplexMatrix AlgebraicComplexNumbers::square_acn_matrix_from_ints(const std::vector<s64>& ints) {
    u64 dim = compute_square_matrix_dim_from_1d_repr(ints);

    ComplexMatrix result(dim, dim);
    AlgebraicComplexNumber4* matrix_slots = new AlgebraicComplexNumber4[dim*dim];

    for (u64 elem_idx = 0; elem_idx < ints.size(); elem_idx++) {
        s64 elem_int_value = ints[elem_idx];
        mpz_set_si(matrix_slots[elem_idx].a, elem_int_value);
    }

    result.data = matrix_slots;
    return result;
}

ComplexMatrix AlgebraicComplexNumbers::row_from_ints(const std::vector<s64>& row_data) {
    AlgebraicComplexNumber4* matrix_slots = new AlgebraicComplexNumber4[row_data.size()];

    for (u64 elem_idx = 0; elem_idx < row_data.size(); elem_idx++) {
        s64 elem_int_value = row_data[elem_idx];
        mpz_set_si(matrix_slots[elem_idx].a, elem_int_value);
    }

    ComplexMatrix result(1, row_data.size(), matrix_slots);
    return result;
}

ComplexMatrix AlgebraicComplexNumbers::column_from_ints(const std::vector<s64>& column_data) {
    AlgebraicComplexNumber4* matrix_slots = new AlgebraicComplexNumber4[column_data.size()];

    for (u64 elem_idx = 0; elem_idx < column_data.size(); elem_idx++) {
        s64 elem_int_value = column_data[elem_idx];
        mpz_set_si(matrix_slots[elem_idx].a, elem_int_value);
    }

    ComplexMatrix result(column_data.size(), 1, matrix_slots);
    return result;
}

std::ostream& AlgebraicComplexNumbers::operator<<(std::ostream& os, const ComplexMatrix& matrix) {
    os << "[";
    if (matrix.height > 1) os << "\n";
    for (u64 row_idx = 0; row_idx < matrix.height; row_idx++) {
        for (u64 col_idx = 0; col_idx < matrix.width; col_idx++) {
            os << matrix.at(row_idx, col_idx) << ", ";
        }

        if (matrix.height > 1) os << "\n";
    }
    os << "]";

    return os;
}

