#pragma once

#include <cassert>
#include <cmath>
#include <gmp-x86_64.h>
#include <iostream>
#include <ostream>
#include <vector>

#include "gmp.h"
#include "basics.hpp"

namespace AlgebraicComplexNumbers {
    struct ComplexNumber {
        float real, im;
    };

    struct FixedPrecisionComplexNumber {
        s64 a, b, c, d, k;

        FixedPrecisionComplexNumber(s64 a, s64 b, s64 c, s64 d, s64 k) :
            a(a), b(b), c(c), d(d), k(k) {}

        bool operator==(const FixedPrecisionComplexNumber& other) const {
            return (a == other.a) && (b == other.b) && (c == other.c) && (d == other.d) && (k == other.k);
        }
    };

    std::ostream& operator<<(std::ostream& os, const ComplexNumber& number);
    std::ostream& operator<<(std::ostream& os, const FixedPrecisionComplexNumber& number);

    /**
     * Algebraic representation of a complex number with arbitrary precision. Number C is represented
     * as a 5-tuple (a, b, c, d, k) such that: C = (1/sqrt(2))^k (a + bO + cO^2 + dO^3) where
     * O (Omega) is e^(PI*i/4).
     */
    struct AlgebraicComplexNumber4 {
        mpz_t a, b, c, d, k;

        AlgebraicComplexNumber4() {
            mpz_inits(a, b, c, d, k, 0);
        }

        AlgebraicComplexNumber4(s64 a, s64 b, s64 c, s64 d, s64 k = 1) {
            mpz_init_set_si(this->a, a);
            mpz_init_set_si(this->b, b);
            mpz_init_set_si(this->c, c);
            mpz_init_set_si(this->d, d);
            mpz_init_set_si(this->k, k);
        }

        AlgebraicComplexNumber4(const AlgebraicComplexNumber4& other) {
            mpz_init_set(this->a, other.a);
            mpz_init_set(this->b, other.b);
            mpz_init_set(this->c, other.c);
            mpz_init_set(this->d, other.d);
            mpz_init_set(this->k, other.k);
        }

        AlgebraicComplexNumber4(AlgebraicComplexNumber4&& other) {
            mpz_init_set(this->a, other.a);
            mpz_init_set(this->b, other.b);
            mpz_init_set(this->c, other.c);
            mpz_init_set(this->d, other.d);
            mpz_init_set(this->k, other.k);
        };

        ~AlgebraicComplexNumber4() {
            mpz_clears(a, b, c, d, k, 0);
        }

        // Useful constants:
        static AlgebraicComplexNumber4 ONE_OVER_SQRT2() {
            return AlgebraicComplexNumber4(1, 0, 0, 0, 1);
        }

        static AlgebraicComplexNumber4 ONE() {
            return AlgebraicComplexNumber4(1, 0, 0, 0, 0);
        }

        static AlgebraicComplexNumber4 ZERO() {
            return AlgebraicComplexNumber4(0, 0, 0, 0, 0);
        }

        AlgebraicComplexNumber4 operator*(const AlgebraicComplexNumber4& other) const {

            AlgebraicComplexNumber4 result;

            mpz_t immediate;
            mpz_init(immediate);

            mpz_add(result.k, this->k, other.k);

            { // (this->a*other.a) - (this->b*other.d) - (this->c*other.c) - (this->d*other.b)
                mpz_mul(immediate, this->a, other.a);
                mpz_add(result.a, result.a, immediate);

                mpz_mul(immediate, this->b, other.d);
                mpz_sub(result.a, result.a, immediate);

                mpz_mul(immediate, this->c, other.c);
                mpz_sub(result.a, result.a, immediate);

                mpz_mul(immediate, this->d, other.b);
                mpz_sub(result.a, result.a, immediate);
            }

            { // (this->a*other.b) + (this->b*other.a) - (this->c*other.d) - (this->d*other.c)
                mpz_mul(immediate, this->a, other.b);
                mpz_add(result.b, result.b, immediate);

                mpz_mul(immediate, this->b, other.a);
                mpz_add(result.b, result.b, immediate);

                mpz_mul(immediate, this->c, other.d);
                mpz_sub(result.b, result.b, immediate);

                mpz_mul(immediate, this->d, other.c);
                mpz_sub(result.b, result.b, immediate);
            }

            { // (this->a*other.c) + (this->b*other.b) + (this->c*other.a) - (this->d*other.d)
                mpz_mul(immediate, this->a, other.c);
                mpz_add(result.c, result.c, immediate);

                mpz_mul(immediate, this->b, other.b);
                mpz_add(result.c, result.c, immediate);

                mpz_mul(immediate, this->c, other.a);
                mpz_add(result.c, result.c, immediate);

                mpz_mul(immediate, this->d, other.d);
                mpz_sub(result.c, result.c, immediate);
            }

            { // result_d = (this->a*other.d) + (this->b*other.c) + (this->c*other.b) + (this->d*other.a);
                mpz_mul(immediate, this->a, other.d);
                mpz_add(result.d, result.d, immediate);

                mpz_mul(immediate, this->b, other.c);
                mpz_add(result.d, result.d, immediate);

                mpz_mul(immediate, this->c, other.b);
                mpz_add(result.d, result.d, immediate);

                mpz_mul(immediate, this->d, other.a);
                mpz_add(result.d, result.d, immediate);
            }

            mpz_clear(immediate);

            if (result.is_zero()) {
                mpz_set_si(result.k, 0);
            }

            return result;
        }

        AlgebraicComplexNumber4 rescale(const mpz_t larger_k) const {
            assert(mpz_cmp(larger_k, this->k) >= 0);

            s64 scale_difference_int = mpz_get_si(larger_k) - mpz_get_si(this->k);
            s64 half_scale_diff = scale_difference_int / 2;

            AlgebraicComplexNumber4 rescaled;

            mpz_mul_2exp(rescaled.a, this->a, half_scale_diff);
            mpz_mul_2exp(rescaled.b, this->b, half_scale_diff);
            mpz_mul_2exp(rescaled.c, this->c, half_scale_diff);
            mpz_mul_2exp(rescaled.d, this->d, half_scale_diff);

            if (scale_difference_int % 2) { // Multiply by sqrt(2) if needed
                mpz_t imm;
                mpz_init(imm);

                AlgebraicComplexNumber4 multiplied_by_sqrt2;

                // Use scale_difference as immediate value
                // COMPUTE: s64 new_a = -rescaled.b - rescaled.d;
                mpz_set_ui(imm, 0);
                mpz_sub(imm, imm, rescaled.b);
                mpz_sub(multiplied_by_sqrt2.a, imm, rescaled.d);

                // COMPUTE: s64 new_b = rescaled.a + rescaled.c;
                mpz_add(multiplied_by_sqrt2.b, rescaled.a, rescaled.c);

                // COMPUTE: s64 new_c = rescaled.b + rescaled.d;
                mpz_add(multiplied_by_sqrt2.c, rescaled.b, rescaled.d);

                // COMPUTE: s64 new_d = rescaled.c - rescaled.a;
                mpz_sub(multiplied_by_sqrt2.d, rescaled.c, rescaled.a);

                mpz_clear(imm);

                return multiplied_by_sqrt2;
            }

            return rescaled;
        }

        AlgebraicComplexNumber4 operator-() const {
            AlgebraicComplexNumber4 result;

            mpz_neg(result.a, this->a);
            mpz_neg(result.b, this->b);
            mpz_neg(result.c, this->c);
            mpz_neg(result.d, this->d);
            mpz_set(result.k, this->k);

            return result;
        }

        AlgebraicComplexNumber4 operator+(const AlgebraicComplexNumber4& other) const {
            const AlgebraicComplexNumber4* smaller = this;
            const AlgebraicComplexNumber4* larger  = &other;
            if (mpz_cmp(larger->k, smaller->k) >= 1) {
                std::swap(smaller, larger);
            }

            AlgebraicComplexNumber4 larger_rescaled = larger->rescale(smaller->k);

            AlgebraicComplexNumber4 result;

            mpz_add(result.a, smaller->a, larger_rescaled.a);
            mpz_add(result.b, smaller->b, larger_rescaled.b);
            mpz_add(result.c, smaller->c, larger_rescaled.c);
            mpz_add(result.d, smaller->d, larger_rescaled.d);
            mpz_set(result.k, smaller->k);

            if (result.is_zero()) {
                mpz_set_si(result.k, 0);
            }

            return result;
        }

        AlgebraicComplexNumber4 operator-(const AlgebraicComplexNumber4& other) const {
            auto other_negated = -other;
            return *this + other_negated;
        }

        void operator=(const AlgebraicComplexNumber4& other) {
            mpz_set(this->a, other.a);
            mpz_set(this->b, other.b);
            mpz_set(this->c, other.c);
            mpz_set(this->d, other.d);
            mpz_set(this->k, other.k);
        }

        void operator+=(const AlgebraicComplexNumber4& other) {
            auto addition_result = *this + other;

            mpz_set(this->a, addition_result.a);
            mpz_set(this->b, addition_result.b);
            mpz_set(this->c, addition_result.c);
            mpz_set(this->d, addition_result.d);
            mpz_set(this->k, addition_result.k);
        }

        bool operator==(const AlgebraicComplexNumber4& other) const {
            return (mpz_cmp(this->a, other.a) == 0) &&
                   (mpz_cmp(this->b, other.b) == 0) &&
                   (mpz_cmp(this->c, other.c) == 0) &&
                   (mpz_cmp(this->d, other.d) == 0) &&
                   (mpz_cmp(this->k, other.k) == 0);
        }

        bool operator!=(const AlgebraicComplexNumber4& other) const {
            return !(*this == other);
        }

        bool is_zero() const {
            return (mpz_cmp_ui(this->a, 0) == 0) &&
                   (mpz_cmp_ui(this->b, 0) == 0) &&
                   (mpz_cmp_ui(this->c, 0) == 0) &&
                   (mpz_cmp_ui(this->d, 0) == 0);
        }

        bool is_integer() const {
            bool is_real = (mpz_cmp_ui(this->b, 0) == 0) &&
                           (mpz_cmp_ui(this->c, 0) == 0) &&
                           (mpz_cmp_ui(this->d, 0) == 0);

            if (!is_real) return false;

            s64 k = mpz_get_si(this->k);
            if (k > 0) {
                if (k & static_cast<s64>(1)) return false;
                s64 dividing_power_of_two = k >> 1;

                s64 a = mpz_get_si(this->a);
                if (a == 0) return true;

                s64 a_two_power = 0;

                while (!(a & 1)) {
                    a_two_power += 1;
                    a >>= 1;
                }

                return (a_two_power >= dividing_power_of_two);
            } else {
                bool multiplying_by_sqrt2 = std::abs(k) & 1;
                return !multiplying_by_sqrt2;
            }
        }

        void swap(AlgebraicComplexNumber4& other) noexcept {
            mpz_swap(this->a, other.a);
            mpz_swap(this->b, other.b);
            mpz_swap(this->c, other.c);
            mpz_swap(this->d, other.d);
            mpz_swap(this->k, other.k);
        }

        ComplexNumber into_approx() const {
            float one_over_sqrt2 = 1.0f / std::sqrt(2);

            float real = static_cast<float>(mpz_get_si(this->a)) + one_over_sqrt2*(static_cast<float>(mpz_get_si(b)) - static_cast<float>(mpz_get_si(this->d)));
            float im   = static_cast<float>(mpz_get_si(this->c)) + one_over_sqrt2*(static_cast<float>(mpz_get_si(b)) + static_cast<float>(mpz_get_si(this->d)));

            real *= std::pow(one_over_sqrt2, mpz_get_si(this->k));
            im   *= std::pow(one_over_sqrt2, mpz_get_si(this->k));

            return {.real = real, .im = im};
        }

        FixedPrecisionComplexNumber into_fixed_precision() const {
            return FixedPrecisionComplexNumber (
                mpz_get_si(this->a),
                mpz_get_si(this->b),
                mpz_get_si(this->c),
                mpz_get_si(this->d),
                mpz_get_si(this->k)
            );
        }
        /**
         *  Return current number C = a + b(e^{\pi * i / 4}) + c(e^{2 \pi * i / 4}) + d(e^{3 *\pi * i / 4})
         *  as (e, f, g, h, k) such that C = e + f*(1/sqrt(2)) + g*i + h*(1/sqrt(2))*i. 
         */
        FixedPrecisionComplexNumber into_sqrt2_repr() const {
            auto fp = this->into_fixed_precision();

            s64 e = fp.a;
            s64 f = fp.b - fp.d;
            s64 g = fp.c;
            s64 h = fp.b + fp.d;
            return FixedPrecisionComplexNumber(e, f, g, h, fp.k);
        }

        void normalize() {
            if (this->is_zero()) {
                mpz_set_ui(this->k, 0);
            }

            u64 a_2pow = mpz_scan1(this->a, 0);
            u64 b_2pow = mpz_scan1(this->b, 0);
            u64 c_2pow = mpz_scan1(this->c, 0);
            u64 d_2pow = mpz_scan1(this->d, 0);

            u64 spare_2pow;
            {
                u64 spare_2pow_a = std::min(a_2pow, b_2pow);
                u64 spare_2pow_b = std::min(c_2pow, d_2pow);
                spare_2pow = std::min(spare_2pow_a, spare_2pow_b);
            }

            if (spare_2pow == 0) return; // There is no common power of 2 dividing all of the components (multiplying omega)

            s64 available_k = mpz_get_ui(this->k) / 2;
            if (available_k < 0) return;

            u64 normalization_2pow = std::min(spare_2pow, static_cast<u64>(available_k));

            // Divide (a, b, c, d) by the agreed power of 2
            mpz_div_2exp(this->a, this->a, normalization_2pow);
            mpz_div_2exp(this->b, this->b, normalization_2pow);
            mpz_div_2exp(this->c, this->c, normalization_2pow);
            mpz_div_2exp(this->d, this->d, normalization_2pow);

            // Subtract from k what has been used in in the process
            mpz_sub_ui(this->k, this->k, 2*normalization_2pow);
        }
    };

    std::ostream& operator<<(std::ostream& os, const AlgebraicComplexNumber4& number);

    struct ComplexMatrix {
        u64 height, width;
        AlgebraicComplexNumber4* data = nullptr;

        ComplexMatrix(u64 height, u64 width, AlgebraicComplexNumber4* data_ptr = nullptr) :
            height(height), width(width), data(data_ptr)
        {
            if (this->data == nullptr) {
                this->data = new AlgebraicComplexNumber4[this->width*this->height]; // Zero initialized
            }
        }

        ComplexMatrix(const ComplexMatrix& other) {
            this->width  = other.width;
            this->height = other.height;

            this->data   = new AlgebraicComplexNumber4[this->width*this->height];

            for (u64 idx = 0; idx < other.width*other.height; idx++) {
                this->data[idx] = other.data[idx];
            }
        }

        ComplexMatrix(ComplexMatrix&& other) {
            this->width  = other.width;
            this->height = other.height;

            this->data = other.data;
            other.data = nullptr;
        }

        ComplexMatrix operator*(const ComplexMatrix& other) const {
            assert(this->width == other.height);

            u64 result_height = this->height;
            u64 result_width  = other.width;

            AlgebraicComplexNumber4* result_data = new AlgebraicComplexNumber4[result_height*result_width];
            u64 target_cell_idx = 0;

            for (u64 row_idx = 0; row_idx < this->height; row_idx++) {
                for (u64 col_idx = 0; col_idx < other.width; col_idx++) {

                    AlgebraicComplexNumber4 dot_product;
                    for (u64 elem_idx = 0; elem_idx < this->width; elem_idx++) {
                        AlgebraicComplexNumber4 fragment = this->at(row_idx, elem_idx) * other.at(elem_idx, col_idx);
                        dot_product += fragment;
                    }

                    result_data[target_cell_idx] = dot_product;
                    target_cell_idx += 1;
                }
            }

            return ComplexMatrix(result_height, result_width, result_data);
        }

        AlgebraicComplexNumber4& at(u64 row_idx, u64 col_idx) const {
            return this->data[row_idx*this->width + col_idx];
        }

        void set(u64 row_idx, u64 col_idx, const AlgebraicComplexNumber4& value) {
            this->data[row_idx*this->width + col_idx] = value;
        }

        u64 find_nonzero_elem_in_row(u64 row_idx) const {
            assert (row_idx < this->height);

            for (u64 elem_idx = 0; elem_idx < this->width; elem_idx++) {
                AlgebraicComplexNumber4& elem = this->at(row_idx, elem_idx);
                if (!elem.is_zero()) {
                    return elem_idx;
                }
            }

            return this->width;
        }

        ComplexMatrix extract_row(u64 row_idx) const {
            ComplexMatrix row (1, this->width);
            for (u64 column_idx = 0; column_idx < this->width; column_idx++) {
                auto& elem_value = this->at(row_idx, column_idx);
                row.set(0, column_idx, elem_value);
            }
            return row;
        }


        void insert_row_at(const ComplexMatrix& row, u64 row_idx, bool skip_shifting_subsequent_rows = false) {
            assert(row.width == this->width);

            if (!skip_shifting_subsequent_rows) {
                s64 last_elem_idx  = (this->height - 1) * this->width - 1;
                s64 first_elem_idx = row_idx * this->width;

                for (s64 elem_idx = last_elem_idx; elem_idx >= first_elem_idx; elem_idx--) {
                    this->data[elem_idx + this->width] = this->data[elem_idx];
                }
            }

            for (u64 elem_idx = 0; elem_idx < row.width; elem_idx++) {
                this->data[row_idx*this->width + elem_idx] = row.at(0, elem_idx);
            }
        }
        /**
         * @Note: We are copying the row_to_subtract_coef in the function invocation, becasue if we were storing just a reference to some external
         * number (e.g. matrix cell) then the subtraction we are doing might in fact modify the coefficient during the subtraction for-loop.
         */
        void subtract_from_ith_row(u64 row_idx, AlgebraicComplexNumber4& row_coef, ComplexMatrix& rows_to_subtract, u64 row_to_subtract_idx, AlgebraicComplexNumber4 row_to_subtract_coef) const {
            for (u64 elem_idx = 0; elem_idx < this->width; elem_idx++) {
                auto subtractee_weighted = this->at(row_idx, elem_idx) * row_coef;
                auto subtractor_weighted = rows_to_subtract.at(row_to_subtract_idx, elem_idx) * row_to_subtract_coef;

                // std::cout << "Subtractee: " << this->at(row_idx, elem_idx) << "*" << row_coef << " = " << subtractee_weighted << "\n";
                // std::cout << "Subtractor: " << rows_to_subtract.at(row_to_subtract_idx, elem_idx) << "*" << row_to_subtract_coef << " = " << subtractor_weighted << "\n";

                auto result_elem = subtractee_weighted - subtractor_weighted;

                this->data[row_idx*this->width + elem_idx] = result_elem;
            }
        }

        bool contains_only_zeros() const {
            for (u64 idx = 0; idx < this->width*this->height; idx++) {
                if (!this->data[idx].is_zero()) return false;
            }
            return true;
        }

        bool operator==(const ComplexMatrix& other) const {
            if (width != other.width || height != other.height) return false;

            for (u64 idx = 0; idx < height*width; idx++) {
                if (data[idx] != other.data[idx]) return false;
            }

            return true;
        }

        ~ComplexMatrix() {
            if (this->data != nullptr) delete[] this->data;
        }
    };

    std::ostream& operator<<(std::ostream& os, const ComplexMatrix& matrix);

    /**
     *  Add a given row into the matrix in row-echelon reduced form.
     */
    s64 add_row_to_row_echelon_matrix(ComplexMatrix& matrix, const ComplexMatrix& row);

    /**
     *  Add a given row into the matrix in row-echelon reduced form. The row is modified in process
     *  into the vector that is eventually inserted into the matrix.
     */
    s64 add_row_to_row_echelon_matrix_no_copy(ComplexMatrix& matrix, ComplexMatrix& row);

    ComplexMatrix square_acn_matrix_from_ints(const std::vector<s64>& ints);
    ComplexMatrix row_from_ints(const std::vector<s64>& row_data);
    ComplexMatrix column_from_ints(const std::vector<s64>& column_data);
}
