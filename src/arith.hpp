#pragma once

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <gmp-x86_64.h>
#include <iostream>
#include <ostream>
#include <utility>
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

        mpz_t& at(s64 idx) {
            switch (idx) {
                case 0: return this->a;
                case 1: return this->b;
                case 2: return this->c;
                case 3: return this->d;
            }
            assert(false);
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

    struct DenseNumberStore {
        s64 width;
        mpz_t* numbers;

        struct Iterator {
            using iterator_category = std::forward_iterator_tag;
            using difference_type   = std::ptrdiff_t;
            using value_type        = std::pair<s64, mpz_t&>;
            using pointer           = mpz_t*;
            using reference         = mpz_t&;

            s64    idx;
            mpz_t* numbers;

            Iterator(s64 idx, mpz_t* nums) : idx(idx), numbers(nums) {}

            value_type operator*() const { return std::pair<s64, mpz_t&>(this->idx, numbers[idx]); }
            pointer   operator->() { return &numbers[idx]; }
            Iterator& operator++() { idx++; return *this; }
            Iterator  operator++(int) { Iterator tmp = *this; ++(*this); return tmp; }

            friend bool operator== (const Iterator& a, const Iterator& b) { return a.numbers == b.numbers && a.idx == b.idx; };
            friend bool operator!= (const Iterator& a, const Iterator& b) { return a.numbers != b.numbers || a.idx != b.idx; }
        };

        struct ConstIterator {
            using iterator_category = std::forward_iterator_tag;
            using difference_type   = std::ptrdiff_t;
            using value_type        = std::pair<s64, const mpz_t&>;
            using pointer           = const mpz_t*;
            using reference         = mpz_t&;

            s64    idx;
            const mpz_t* numbers;

            ConstIterator(s64 idx, mpz_t* nums) : idx(idx), numbers(nums) {}

            value_type     operator*() const { return std::pair<s64, const mpz_t&>(this->idx, numbers[idx]); }
            pointer        operator->() { return &numbers[idx]; }
            ConstIterator& operator++() { idx++; return *this; }
            ConstIterator  operator++(int) { ConstIterator tmp = *this; ++(*this); return tmp; }

            friend bool operator== (const ConstIterator& a, const ConstIterator& b) { return a.numbers == b.numbers && a.idx == b.idx; };
            friend bool operator!= (const ConstIterator& a, const ConstIterator& b) { return a.numbers != b.numbers || a.idx != b.idx; }
        };

        DenseNumberStore(s64 w) : width(w) {
            this->numbers = new mpz_t[w];

            for (s64 i = 0; i < this->width; i++) {
                mpz_init(this->numbers[i]);
            }
        }

        DenseNumberStore(const DenseNumberStore& other) : width(other.width) {
            this->numbers = new mpz_t[other.width];

            for (s64 i = 0; i < this->width; i++) {
                mpz_init(this->numbers[i]);
                mpz_set(this->numbers[i], other.numbers[i]);
            }
        }

        DenseNumberStore(DenseNumberStore&& other) : width(other.width) {
            this->numbers = other.numbers;
            other.numbers = nullptr;
        }

        ~DenseNumberStore() {
            if (this->numbers == nullptr) return;

            for (s64 i = 0; i < this->width; i++) {
                mpz_clear(this->numbers[i]);
            }

            delete[] this->numbers;
            this->numbers = nullptr;
        }

        mpz_t& at(s64 index) {
            return this->numbers[index];
        }

        void set_fp(s64 idx, s64 value) {
            mpz_set_si(this->numbers[idx], value);
        }

        Iterator begin() const {
            return Iterator(0, this->numbers);
        }

        ConstIterator cbegin() const {
            return ConstIterator(0, this->numbers);
        }

        Iterator end() const {
            return Iterator(this->width, this->numbers);
        }

        ConstIterator cend() const {
            return ConstIterator(this->width, this->numbers);
        }

        bool operator==(const DenseNumberStore& other) const {
            for (s64 i = 0; i < this->width; i++) {
                if (mpz_cmp(this->numbers[i], other.numbers[i])) {
                    return false;
                }
            }
            return true;
        }

        DenseNumberStore& operator=(const DenseNumberStore& other) {
            for (s64 idx = 0; idx < this->width; idx++) {
                mpz_clear(this->numbers[idx]);
            }
            delete[] this->numbers;

            this->width = other.width;
            this->numbers = new mpz_t[this->width];

            for (s64 idx = 0; idx < this->width; idx++) {
                mpz_init(this->numbers[idx]);
                mpz_set(this->numbers[idx], other.numbers[idx]);
            }

            return *this;
        }
    };

    template <typename NumberStorage>
    struct AlgebraicComplexNumber {
        s64 n;  // Into how many parts we divide the unit half-circle; needs to be a power of 2

        NumberStorage coefficients;
        mpz_t scaling_factor;

        AlgebraicComplexNumber(s64 width) : n(width), coefficients(n) {
            mpz_init(this->scaling_factor);
        }

        AlgebraicComplexNumber(s64 width, s64 scale) : n(width), coefficients(n) {
            mpz_init(this->scaling_factor);
            mpz_set_si(this->scaling_factor, scale);
        }

        AlgebraicComplexNumber(s64 width, const mpz_t scale) : n(width), coefficients(n) {
            mpz_init(this->scaling_factor);
            mpz_set(this->scaling_factor, scale);
        }

        AlgebraicComplexNumber(const AlgebraicComplexNumber& other) : n(other.n), coefficients(other.coefficients) {
            mpz_init(this->scaling_factor);
            mpz_set(this->scaling_factor, other.scaling_factor);
        }

        AlgebraicComplexNumber(AlgebraicComplexNumber&& other) : n(other.n), coefficients(other.coefficients) {
            mpz_init(this->scaling_factor);
            mpz_set(this->scaling_factor, other.scaling_factor);
        }

        ~AlgebraicComplexNumber() {
            mpz_clear(this->scaling_factor);
        }

        std::pair<bool, s64> calc_neighbour(s64 current_idx, s64 offset) {
            offset = offset % (2*this->n); // Normalize, so we do not spin too much

            bool overflowed = false;
            s64 neighbour = (current_idx + offset) % (2*this->n);

            if (neighbour >= this->n) {
                overflowed = true;
                neighbour = (neighbour % this->n);
            }

            if (neighbour < 0) {
                overflowed = true;
                neighbour += (this->n);
            }

            return std::make_pair(overflowed, neighbour);
        }

        /**
         * Rescale the complex number so that this->scaling_factor = desired_scaling_factor
         */
        AlgebraicComplexNumber<NumberStorage> rescale(const mpz_t& desired_scaling_factor) const {
            assert(mpz_cmp(desired_scaling_factor, this->scaling_factor) >= 0);

            if (DEBUG) {
                mpz_t scale_diff;
                mpz_init(scale_diff);
                mpz_sub(scale_diff, desired_scaling_factor, scaling_factor);
                assert(mpz_fits_slong_p(scale_diff));
            }

            s64 scale_diff = mpz_get_si(desired_scaling_factor) - mpz_get_si(this->scaling_factor);
            s64 half_scale_diff = scale_diff / 2;

            AlgebraicComplexNumber rescaled (this->n, desired_scaling_factor);
            for (const auto& [coef_idx, coef] : this->coefficients) {
                mpz_mul_2exp(rescaled.coefficients.at(coef_idx), coef, half_scale_diff);
            }

            if (scale_diff % 2) { // Multiply by sqrt(2) if needed
                // sqrt(2) = e^(2pi*i/n  * n/4) + e^(2pi*i/n  * 3n/4)

                AlgebraicComplexNumber multiplied_by_sqrt2(this->n, desired_scaling_factor);

                auto add_to_friend = [&rescaled, &multiplied_by_sqrt2](s64 idx, s64 friend_offset) {
                    auto [overflowed, friend_idx] = rescaled.calc_neighbour(idx, friend_offset);
                    mpz_t& friend_val = multiplied_by_sqrt2.coefficients.at(friend_idx);
                    if (overflowed) {
                        mpz_sub(friend_val, friend_val, rescaled.coefficients.at(idx));
                    } else {
                        mpz_add(friend_val, friend_val, rescaled.coefficients.at(idx));
                    }
                };

                for (auto [idx, coef] : rescaled.coefficients) {
                    add_to_friend(idx,   n/4);
                    add_to_friend(idx, 3*n/4);
                }

                return multiplied_by_sqrt2;
            }

            return rescaled;
        }

        AlgebraicComplexNumber<NumberStorage> rescale(s64 desired_scale_fp) {
            mpz_t desired_scale;
            mpz_init(desired_scale);
            mpz_set_si(desired_scale, desired_scale_fp);
            auto result = this->rescale(desired_scale);
            mpz_clear(desired_scale);
            return result;
        }

        AlgebraicComplexNumber<NumberStorage> operator+(const AlgebraicComplexNumber<NumberStorage>& other) const {
            const AlgebraicComplexNumber<NumberStorage>* num_with_larger_scale  = this;
            const AlgebraicComplexNumber<NumberStorage>* num_with_smaller_scale = &other;

            if (mpz_cmp(num_with_larger_scale->scaling_factor, num_with_smaller_scale->scaling_factor) < 0) { // e.g. (1/sqrt(2))^-4 vs (1/sqrt(2))^-2
                std::swap(num_with_larger_scale, num_with_smaller_scale);
            }

            AlgebraicComplexNumber<NumberStorage> smaller_rescaled = num_with_smaller_scale->rescale(num_with_larger_scale->scaling_factor);
            for (const auto& [idx, val] : num_with_larger_scale->coefficients) {
                mpz_add(smaller_rescaled.coefficients.at(idx), smaller_rescaled.coefficients.at(idx), val);
            }

            return smaller_rescaled;
        }

        AlgebraicComplexNumber<NumberStorage> operator-(const AlgebraicComplexNumber<NumberStorage>& other) const {
            if (mpz_cmp(this->scaling_factor, other.scaling_factor) < 0) {
                AlgebraicComplexNumber<NumberStorage> this_rescaled = this->rescale(other.scaling_factor);
                for (const auto& [idx, val] : other.coefficients) {
                    mpz_sub(this_rescaled.coefficients.at(idx), this_rescaled.coefficients.at(idx), val);
                }
                return this_rescaled;
            } else {
                AlgebraicComplexNumber<NumberStorage> other_rescaled = other.rescale(this->scaling_factor);
                AlgebraicComplexNumber<NumberStorage> result (*this);

                for (const auto& [idx, val] : other_rescaled.coefficients) {
                    mpz_sub(result.coefficients.at(idx), result.coefficients.at(idx), val);
                }
                return result;
            }
        }

        AlgebraicComplexNumber<NumberStorage> operator*(const AlgebraicComplexNumber<NumberStorage>& other) const {
            AlgebraicComplexNumber<NumberStorage> result (this->n);

            mpz_mul(result.scaling_factor, this->scaling_factor, other.scaling_factor);

            mpz_t imm;
            mpz_init(imm);

            for (const auto& [this_idx, this_val] : this->coefficients) {
                for (const auto& [other_idx, other_val] : other.coefficients) {
                    mpz_mul(imm, this_val, other_val);

                    s64 dst_idx = this_idx + other_idx;
                    if (dst_idx >= this->n) {
                        mpz_neg(imm, imm);
                        dst_idx = dst_idx - n;
                    }

                    mpz_t& target = result.coefficients.at(dst_idx);
                    mpz_add(target, target, imm);
                }
            }

            mpz_clear(imm);

            return result;
        }

        AlgebraicComplexNumber<NumberStorage> operator-() {
            AlgebraicComplexNumber<NumberStorage> result(*this);
            for (auto [idx, coef] : result.coefficients) {
                mpz_neg(coef, coef);
            }
            return result;
        }

        bool operator==(const AlgebraicComplexNumber<NumberStorage>& other) const {
            if (this->n != other.n) {
                return false;
            }

            bool scale_eq = (mpz_cmp(this->scaling_factor, other.scaling_factor) == 0);
            if (!scale_eq) return false;
            bool coefs_eq = (this->coefficients == other.coefficients);

            return coefs_eq;
        }

        bool operator!=(const AlgebraicComplexNumber<NumberStorage>& other) const {
            return !(*this == other);
        }

        bool is_zero() const {
            bool result = true;

            for (auto [idx, coef] : this->coefficients) {
                if (mpz_sgn(coef) != 0) {
                    result  = false;
                    break;
                }
            }

            return result;
        }

        void operator+=(const AlgebraicComplexNumber<NumberStorage>& other) {
            // @Optimize: This is a silly implementation (needless memory allocations)
            auto sum = *this + other;
            *this = sum;
        }

        AlgebraicComplexNumber<NumberStorage> operator=(const AlgebraicComplexNumber<NumberStorage>& other) {
            this->coefficients = other.coefficients;
            this->n = other.n;
            mpz_set(this->scaling_factor, other.scaling_factor);
            return *this;
        }

        mpz_t& at(s64 idx) {
            return this->coefficients.at(idx);
        }
    };

    std::ostream& operator<<(std::ostream& out, const AlgebraicComplexNumber<DenseNumberStore>& val);
    AlgebraicComplexNumber<DenseNumberStore> from_fp_vector(std::vector<s64> coefs, s64 scale);
}

template class AlgebraicComplexNumbers::AlgebraicComplexNumber<AlgebraicComplexNumbers::DenseNumberStore>;
