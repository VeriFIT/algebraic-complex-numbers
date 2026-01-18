#pragma once

#include "arith.hpp"
#include <vector>

using namespace AlgebraicComplexNumbers;

namespace AlgebraicComplexNumbers {
    template <typename Number>
    struct ComplexMatrix {
        u64 height, width;
        Number* data = nullptr;

        ComplexMatrix(u64 height, u64 width, Number* data_ptr = nullptr) :
            height(height), width(width), data(data_ptr)
        {
            if (this->data == nullptr) {
                this->data = new Number[this->width*this->height]; // Zero initialized
            }
        }

        ComplexMatrix(const ComplexMatrix<Number>& other) {
            this->width  = other.width;
            this->height = other.height;

            this->data   = new Number[this->width*this->height];

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

            Number *result_data = new Number[result_height*result_width];
            u64 target_cell_idx = 0;

            for (u64 row_idx = 0; row_idx < this->height; row_idx++) {
                for (u64 col_idx = 0; col_idx < other.width; col_idx++) {

                    Number dot_product;
                    for (u64 elem_idx = 0; elem_idx < this->width; elem_idx++) {
                        Number fragment = this->at(row_idx, elem_idx) * other.at(elem_idx, col_idx);
                        dot_product += fragment;
                    }

                    result_data[target_cell_idx] = dot_product;
                    target_cell_idx += 1;
                }
            }

            return ComplexMatrix(result_height, result_width, result_data);
        }

        Number& at(u64 row_idx, u64 col_idx) const {
            return this->data[row_idx*this->width + col_idx];
        }

        void set(u64 row_idx, u64 col_idx, const Number &value) {
            this->data[row_idx*this->width + col_idx] = value;
        }

        u64 find_nonzero_elem_in_row(u64 row_idx) const {
            assert (row_idx < this->height);

            for (u64 elem_idx = 0; elem_idx < this->width; elem_idx++) {
                Number& elem = this->at(row_idx, elem_idx);
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
        void subtract_from_ith_row(u64 row_idx, Number& row_coef, ComplexMatrix& rows_to_subtract, u64 row_to_subtract_idx, Number row_to_subtract_coef) const {
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

    template <typename T>
    std::ostream& operator<<(std::ostream& os, const ComplexMatrix<T>& matrix) {
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

    /**
     *  Add a given row into the matrix in row-echelon reduced form.
     */
    template <typename T>
    s64 add_row_to_row_echelon_matrix(ComplexMatrix<T>& matrix, const ComplexMatrix<T>& row) {
        ComplexMatrix<T> row_copy (1, row.width); // Make a local copy, since we will be modifying it
        for (auto row_elem_idx = 0; row_elem_idx < row.width; row_elem_idx++) {
            auto& elem_value = row.at(0, row_elem_idx);
            row_copy.set(0, row_elem_idx, elem_value);
        }
        return add_row_to_row_echelon_matrix_no_copy(matrix, row_copy);
    }

    /**
     *  Add a given row into the matrix in row-echelon reduced form. The row is modified in process
     *  into the vector that is eventually inserted into the matrix.
     */
    template <typename T>
    s64 add_row_to_row_echelon_matrix_no_copy(ComplexMatrix<T>& matrix, ComplexMatrix<T>& row) {
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
}
