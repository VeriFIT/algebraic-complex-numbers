#include "arith.hpp"
#include "matrix.hpp"

#define CATCH_CONFIG_MAIN

#include <catch2/catch.hpp>

using AlgebraicComplexNumbers::AlgebraicComplexNumber4;
using AlgebraicComplexNumbers::FixedPrecisionComplexNumber;
using AlgebraicComplexNumbers::ComplexMatrix;

using AlgebraicComplexNumbers::row_from_ints;
using AlgebraicComplexNumbers::square_acn_matrix_from_ints;

TEST_CASE( "Simple multiplication", "[Algebraic complex numbers]" ) {
    AlgebraicComplexNumber4 acn(0, 1, 0, 0, 0);
    auto result = acn * acn;
    auto fp_result = result.into_fixed_precision();

    REQUIRE(fp_result.a == 0);
    REQUIRE(fp_result.b == 0);
    REQUIRE(fp_result.c == 1);
    REQUIRE(fp_result.d == 0);
    REQUIRE(fp_result.k == 0);
}


TEST_CASE( "Multiplication comutativity", "[Algebraic complex numbers]" ) {
    AlgebraicComplexNumber4 left  (1, 2, 3, 4, 0);
    AlgebraicComplexNumber4 right (0, 1, 2, 3, 1);
    auto left_imm  = left*right;
    auto right_imm = right*left;
    auto result = left_imm - right_imm;

    REQUIRE(result.is_zero());
}

TEST_CASE( "Scaling during addition", "[Algebraic complex numbers]" ) {
    {
        AlgebraicComplexNumber4 left (1, 2, 0, 0, 0);
        AlgebraicComplexNumber4 right (-2, -5, 0, 0, -2);

        AlgebraicComplexNumber4 result = left + right;
        FixedPrecisionComplexNumber direct_result = result.into_fixed_precision();

        // Result should be  -9 - 6/sqrt(2) - 10i + 2i/sqrt(2)
        FixedPrecisionComplexNumber expected (-3, -8, 0, 0, 0);
        REQUIRE(direct_result == expected);
    }

    {
        AlgebraicComplexNumber4 left (1, 2, 0, 0, 0);
        AlgebraicComplexNumber4 right (-2, -5, 0, 0, -3);

        AlgebraicComplexNumber4 result = left + right;
        auto fp_result = result.into_fixed_precision();

        FixedPrecisionComplexNumber expected_result (11, -2, -10, 4, 0);
        REQUIRE(expected_result == fp_result);
    }

    {
        // Represents: 1 + 2/sqrt(2) + 2*i/sqrt(2) + 1*i
        AlgebraicComplexNumber4 left(1, 2, 1, 0, 0);

        // Represents: -7 -4/sqrt(2) -7i/sqrt(2)
        AlgebraicComplexNumber4 right(-2, -5, 0, 2, -1);

        auto result = left + right;
        auto fp_result = result.into_fixed_precision();

        FixedPrecisionComplexNumber expected_result (4, 0, -2, 2, 0);

        REQUIRE(expected_result == fp_result);
    }
}

TEST_CASE( "Conversion into direct representation", "[Algebraic complex numbers]") {
    // Represents: 1 + -2/sqrt(2) + 3i + 6i/sqrt(2)
    {
        AlgebraicComplexNumber4 number (1, 2, 3, 4, 0);
        FixedPrecisionComplexNumber direct_repr = number.into_sqrt2_repr();

        FixedPrecisionComplexNumber expected_result = {1, -2, 3, 6, 0};
        REQUIRE(expected_result == direct_repr);
    }

    {
        AlgebraicComplexNumber4 number (0, 1, 2, 3, 1);
        FixedPrecisionComplexNumber direct_repr = number.into_sqrt2_repr();

        FixedPrecisionComplexNumber expected_result = {0, -2, 2, 4, 1};
        REQUIRE(expected_result == direct_repr);
    }

    {
        AlgebraicComplexNumber4 number (-2, -5, 2, -3, -3);
        FixedPrecisionComplexNumber direct_repr = number.into_sqrt2_repr();

        FixedPrecisionComplexNumber expected_result = {-2, -2, 2, -8, -3};
        REQUIRE(expected_result == direct_repr);
    }
}

TEST_CASE( "Add row to a row-echelon-form matrix", "[ACN Matrix]") {
    {
        ComplexMatrix matrix = square_acn_matrix_from_ints({
            0, 0, 0,
            0, 0, 0,
            0, 0, 0,
        });
        ComplexMatrix row = row_from_ints({1, 0, 0});

        s64 row_slot_idx        = add_row_to_row_echelon_matrix(matrix, row);
        s64 subsequent_slot_idx = add_row_to_row_echelon_matrix(matrix, row);

        REQUIRE(row_slot_idx == 0);
        REQUIRE(subsequent_slot_idx == -1);
    }

    {
        ComplexMatrix matrix = square_acn_matrix_from_ints({
            1, 0, 0,
            0, 1, 0,
            0, 0, 1,
        });
        ComplexMatrix row = row_from_ints({1, 2, 3});

        s64 row_slot_idx = add_row_to_row_echelon_matrix(matrix, row);

        REQUIRE(row_slot_idx == -1);
    }

    {
        ComplexMatrix matrix = square_acn_matrix_from_ints({
            1, 0, 0,
            0, 3, 0,
            0, 0, 8,
        });
        ComplexMatrix row = row_from_ints({1, 2, 0});

        s64 row_slot_idx = add_row_to_row_echelon_matrix(matrix, row);

        REQUIRE(row_slot_idx == -1);
    }

    {
        ComplexMatrix matrix = square_acn_matrix_from_ints({
            1, -1, 0,
            0,  0, 0,
            0,  0, 0,
        });
        ComplexMatrix row = row_from_ints({1, -1, 0});

        s64 row_slot_idx = add_row_to_row_echelon_matrix(matrix, row);

        REQUIRE(row_slot_idx == -1);
    }
}

TEST_CASE( "MMUL", "[ACN Matrix]") {
    {
       ComplexMatrix row = row_from_ints({1, -2});
       ComplexMatrix mat = square_acn_matrix_from_ints({
           1, 2,
           3, 4
       });
       auto result = row * mat;

       auto expected = row_from_ints({-5, -6});
       REQUIRE(result == expected);
    }
    {
       ComplexMatrix row = row_from_ints({0, -2, 1, 0});
       ComplexMatrix mat = square_acn_matrix_from_ints({
           0, -2, 1, 0,
           0,  1, 0, 1,
           0,  0, 1, 2,
           0,  0, 0, 0,
       });
       auto result = row * mat;

       auto expected = row_from_ints({0, -2, 1, 0});
       REQUIRE(result == expected);
    }
}
