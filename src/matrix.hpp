#pragma once

#include "arith.hpp"

using namespace AlgebraicComplexNumbers;

namespace AlgebraicComplexNumbers {
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
