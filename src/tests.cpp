#include "arith.hpp"
#include "matrix.hpp"
#include "matrix_constructs.hpp"

#define CATCH_CONFIG_MAIN

#include <catch2/catch.hpp>

using namespace AlgebraicComplexNumbers;

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

TEST_CASE ( "Rescaling with 4", "[Algebraic complex numbers]" ) {
    auto num = from_fp_vector({1, 2, 3, 4}, 0);
    auto rescaled = num.rescale(1);
    auto expected_num = from_fp_vector({-6, -2, -2, 4}, 1);

    REQUIRE(rescaled == expected_num);
}

TEST_CASE ( "Rescaling with width 8", "[Algebraic complex numbers]" ) {
    auto num = from_fp_vector({1, 2, 3, 4, 1, 2, 3, 4}, 0);
    auto rescaled = num.rescale(1);
    auto expected_num = from_fp_vector({-6, -8, 0, 0, 0, 0, 2, 4}, 1);

    REQUIRE(rescaled == expected_num);
}

TEST_CASE ( "Rescaling with width 8, factor=3", "[Algebraic complex numbers]" ) {
    auto num = from_fp_vector({1, 2, 3, 4, 1, 2, 3, 4}, 0);
    auto rescaled = num.rescale(3);
    auto expected_num = from_fp_vector({-12, -16, 0, 0, 0, 0, 4, 8}, 3);

    REQUIRE(rescaled == expected_num);
}

TEST_CASE ( "Rescaling with width 4, factor=-3 > -2", "[Algebraic complex numbers]" ) {
    auto num = from_fp_vector({1, 2, 3, 4}, -3);
    auto rescaled = num.rescale(-2);
    auto expected_num = from_fp_vector({-6, -2, -2, 4}, -2);

    REQUIRE(rescaled == expected_num);
}

TEST_CASE ( "Rescaling with width 8, factor=-3 > -1", "[Algebraic complex numbers]" ) {
    auto num = from_fp_vector({1, 2, 3, 4, 1, 2, 3, 4}, -3);
    auto rescaled = num.rescale(-1);
    auto expected_num = from_fp_vector({2, 4, 6, 8, 2, 4, 6, 8}, -1);

    REQUIRE(rescaled == expected_num);
}

TEST_CASE ("Addition without rescaling", "[Algebraic complex numbers]") {
    auto a = from_fp_vector({1, 2, 3, 4, 1, 2, 3, 4}, -1);
    auto b = from_fp_vector({1, 2, 3, 4, 1, 2, 3, 4}, -1);
    auto r = a + b;

    auto expected = from_fp_vector({2, 4, 6, 8, 2, 4, 6, 8}, -1);
    REQUIRE(expected == r);
}

TEST_CASE ("Addition with rescaling", "[Algebraic complex numbers]") {
    auto a = from_fp_vector({1, 0, 0, 0, 0, 0, 0, 0}, -1);
    auto b = from_fp_vector({0, 1, 0, 0, 0, 0, 0, 0}, 0);
    auto r = a + b;

    auto expected = from_fp_vector({0, 1, 1, 0, 0, 0, 1, 0}, 0);
    REQUIRE(expected == r);
}

TEST_CASE ("Subtraction without rescaling", "[Algebraic complex numbers]") {
    auto a = from_fp_vector({1, 0, 2, 0, 0, 0, 0, 0}, -1);
    auto b = from_fp_vector({0, 0, 3, 4, 0, 0, 0, 0}, -1);
    auto r = a - b;

    auto expected = from_fp_vector({1, 0, -1, -4, 0, 0, 0, 0}, -1);
    REQUIRE(expected == r);
}

TEST_CASE ("Subtraction with rescaling (left operand has smaller sf)", "[Algebraic complex numbers]") {
    auto a = from_fp_vector({0, 1, 0, 0, 0, 0, 0, 0}, -2); // After scaling ==> from_fp_vector({0, 0, 0, 1, 0, 0, 0, 1}, -1);
    auto b = from_fp_vector({0, 0, 3, 4, 0, 0, 0, 0}, -1);
    auto r = a - b;

    auto expected = from_fp_vector({0, 0, -3, -3, 0, 0, 0, 1}, -1);
    REQUIRE(expected == r);
}

TEST_CASE ("Subtraction with rescaling (right operand has smaller sf)", "[Algebraic complex numbers]") {
    auto a = from_fp_vector({0, 0, 3, 4, 0, 0, 0, 0}, -1);
    auto b = from_fp_vector({0, 1, 0, 0, 0, 0, 0, 0}, -2); // After scaling ==> from_fp_vector({0, 0, 0, 1, 0, 0, 0, 1}, -1);
    auto r = a - b;

    auto expected = from_fp_vector({0, 0, 3, 3, 0, 0, 0, -1}, -1);
    REQUIRE(expected == r);
}

TEST_CASE ("Multiplication 0", "[Algebraic complex numbers]") {
    auto a = from_fp_vector({2, 0, 0, 0}, 2);
    auto b = from_fp_vector({3, 0, 0, 0}, -1);
    auto r = a * b;

    auto expected = from_fp_vector({6, 0, 0, 0}, -2);
    REQUIRE(r == expected);
}

TEST_CASE ("Multiplication 1", "[Algebraic complex numbers]") {
    auto a = from_fp_vector({1, 2, 3, 0}, 2);
    auto b = from_fp_vector({2, 3, 4, 5}, -1);
    auto r = a * b;
    // {  2,   3, 4, 5}
    // {-10,   4, 6, 8}
    // {-12, -15, 6, 9}
    auto expected = from_fp_vector({-20, -8, 16, 22}, -2);
    REQUIRE(r == expected);
}

TEST_CASE ("Multiplication 2", "[Algebraic complex numbers]") {
    auto a = from_fp_vector({1, 2, 0, 0, 0, 0, 1, 0}, 2);
    auto b = from_fp_vector({2, 3, 4, 5, 0, 1, 0, 0}, -1);
    auto r = a * b;
    auto expected = from_fp_vector({-2, 2, 10, 12, 10, 1, 4, 3}, -2);
    REQUIRE(r == expected);
}

TEST_CASE ("Negation 0", "[Algebraic complex numbers]") {
    auto a = from_fp_vector({1, 2, 0, 0, 0, 0, 1, 0}, 2);
    auto r = -a;
    auto expected = from_fp_vector({-1, -2, 0, 0, 0, 0, -1, 0}, 2);
    REQUIRE(r == expected);
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
    using ACN = AlgebraicComplexNumber<DenseNumberStore>;
    ACN zero(4);

    MatrixBaker<ACN> baker { &zero };

    {
        ComplexMatrix<ACN> matrix = baker.square_acn_matrix_from_ints({
            0, 0, 0,
            0, 0, 0,
            0, 0, 0,
        });
        ComplexMatrix<ACN> row = baker.row_from_ints({1, 0, 0});

        s64 row_slot_idx        = add_row_to_row_echelon_matrix<ACN>(matrix, row);
        s64 subsequent_slot_idx = add_row_to_row_echelon_matrix<ACN>(matrix, row);

        REQUIRE(row_slot_idx == 0);
        REQUIRE(subsequent_slot_idx == -1);
    }

    {
        ComplexMatrix<ACN> matrix = baker.square_acn_matrix_from_ints({
            1, 0, 0,
            0, 1, 0,
            0, 0, 1,
        });
        ComplexMatrix<ACN> row = baker.row_from_ints({1, 2, 3});

        s64 row_slot_idx = add_row_to_row_echelon_matrix<ACN>(matrix, row);

        REQUIRE(row_slot_idx == -1);
    }

    {
        ComplexMatrix<ACN> matrix = baker.square_acn_matrix_from_ints({
            1, 0, 0,
            0, 3, 0,
            0, 0, 8,
        });
        ComplexMatrix<ACN> row = baker.row_from_ints({1, 2, 0});

        s64 row_slot_idx = add_row_to_row_echelon_matrix<ACN>(matrix, row);

        REQUIRE(row_slot_idx == -1);
    }

    {
        ComplexMatrix<ACN> matrix = baker.square_acn_matrix_from_ints({
            1, -1, 0,
            0,  0, 0,
            0,  0, 0,
        });
        ComplexMatrix<ACN> row = baker.row_from_ints({1, -1, 0});

        s64 row_slot_idx = add_row_to_row_echelon_matrix<ACN>(matrix, row);

        REQUIRE(row_slot_idx == -1);
    }
}

TEST_CASE( "MMUL", "[ACN Matrix]") {
    using ACN = AlgebraicComplexNumber<DenseNumberStore>;
    ACN zero (4);
    MatrixBaker<ACN> baker { &zero };

    {
       auto row = baker.row_from_ints({1, -2});
       auto mat = baker.square_acn_matrix_from_ints({
           1, 2,
           3, 4
       });
       auto result = row * mat;

       auto expected = baker.row_from_ints({-5, -6});
       REQUIRE(result == expected);
    }

    {
       auto row = baker.row_from_ints({0, -2, 1, 0});
       auto mat = baker.square_acn_matrix_from_ints({
           0, -2, 1, 0,
           0,  1, 0, 1,
           0,  0, 1, 2,
           0,  0, 0, 0,
       });
       auto result = row * mat;

       auto expected = baker.row_from_ints({0, -2, 1, 0});
       REQUIRE(result == expected);
    }
}
