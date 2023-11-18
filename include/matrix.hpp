#pragma once

#include <iostream>
#include <utility>
#include <iomanip>
#include "exceptions.hpp"
#include "double_funcs.hpp"


namespace matrix {
    template <typename T = int> class matrix_t;

    namespace detail {
        void add_submatrix_koeff_row(matrix_t<double>& matrix, int row, int add_row, double koeff = 1.0);
        int find_max_elem_submatrix_row(matrix_t<double>& matrix, int start);
    }

    template <typename T> class matrix_t final {
        T** data_;
        T* buffer_;
        int cols_, rows_;

        class proxy_row_t {
            T* row_ = nullptr;
            public:
                proxy_row_t(T* row) : row_(row) {}

                const T& operator[](int col) const { return row_[col]; }
                      T& operator[](int col)       { return row_[col]; }
        };

    public:
        matrix_t(int cols = 1, int rows = 1, T val = T{}) : cols_(cols), rows_(rows),
            data_(new T*[rows]), buffer_(new T[rows * cols]) {
            for (int i = 0; i < rows; ++i) {
                data_[i] = &buffer_[i * cols];

                for (int j = 0; j < cols; ++j) data_[i][j] = val;
            }

        }

        template<typename It>
        matrix_t(int cols, int rows, It start, It fin) : cols_(cols), rows_(rows),
            data_(new T*[rows]), buffer_(new T[rows * cols]) {

            for (int i = 0; i < rows; ++i) {
                data_[i] = &buffer_[i * cols];

                for (int j = 0; j < cols; ++j) {
                    if (start != fin) {
                        data_[i][j] = *start;
                        start++;
                    }
                    else data_[i][j] = 0;
                }
            }
        }

        matrix_t(const matrix_t& matrix) : cols_(matrix.cols_), rows_(matrix.rows_),
            data_(new T*[matrix.rows_]), buffer_(new T[matrix.rows_ * matrix.cols_]) {
            for (int i = 0; i < rows_; ++i) {
                data_[i] = &buffer_[i * cols_];

                for (int j = 0; j < cols_; ++j) data_[i][j] = matrix.data_[i][j];
            }
        }

        matrix_t(matrix_t&& matrix) : cols_(matrix.cols_), rows_(matrix.rows_),
            data_(matrix.data_), buffer_(matrix.buffer_) {
                matrix.data_ = nullptr;
                matrix.buffer_ = nullptr;
            }

        ~matrix_t() { delete [] buffer_; delete [] data_; }

        matrix_t& operator=(const matrix_t& matrix) {
            if (&matrix == this) return *this;

            delete[] data_; delete[] buffer_;

            cols_ = matrix.cols_;
            rows_ = matrix.rows_;
            data_ = new T*[rows_];
            buffer_ = new T[rows_ * cols_];

            for (int i = 0; i < rows_; ++i) {
                data_[i] = &buffer_[i * cols_];

                for (int j = 0; j < cols_; ++j) data_[i][j] = matrix.data_[i][j];
            }
            return *this;
        }

        matrix_t& operator=(matrix_t&& matrix) {
            if (&matrix == this) return *this;

            cols_ = matrix.cols_; rows_ = matrix.rows_;
            std::swap(data_, matrix.data_);
            std::swap(buffer_, matrix.buffer_);

            return *this;
        }

        static matrix_t eye(int n = 1) {
            matrix_t mtx{n, n, 0};

            for (int i = 0; i < n; ++i) {
                mtx.data_[i][i] = 1;
            }
            return mtx;
        }

        matrix_t& negate() & {
            for (int i = 0; i < rows_; ++i) {
                for (int j = 0; j < cols_; ++j) {
                    data_[i][j] *= -1;
                }
            }
            return *this;
        }

        matrix_t& transpose() & {
            T** new_data = new T*[cols_];
            T*  new_buffer = new T[cols_ * rows_];

            for (int i = 0; i < cols_; ++i) new_data[i] = &new_buffer[i * rows_];

            for (int i = 0; i < rows_; ++i) {
                for (int j = 0; j < cols_; ++j) new_data[j][i] = data_[i][j];
            }
            delete[] buffer_;
            delete[] data_;

            data_ = new_data;
            buffer_ = new_buffer;

            int tmp = cols_;
            cols_ = rows_;
            rows_ = tmp;

            return *this;
        }

        T trace() const {
            if (cols_ != rows_) throw matrix_exceptions::MatrixIsNotSquare();

            T ans = data_[0][0];
            for (int i = 1; i < cols_; ++i) ans += data_[i][i];

            return ans;
        }

        double determinant() const {
            if (cols_ != rows_) throw matrix_exceptions::MatrixIsNotSquare();

            return gauss_algorithm();
        }

        void fill_double_matrix(matrix_t<double>& matrix) const {
            for (int i = 0; i < rows_; ++i) {
                for (int j = 0; j < cols_; ++j) matrix[i][j] = data_[i][j];
            }
        }

        void swap_rows(int row1, int row2) {
            T* tmp = data_[row1];
            data_[row1] = data_[row2];
            data_[row2] = tmp;
        }

    private:

        double gauss_algorithm() const {
            matrix_t<double> matrix{cols_, rows_, NAN};

            fill_double_matrix(matrix);

            double det = 1.0;

            for (int i = 0; i < cols_ - 1; ++i) {
                int max_elem_row = detail::find_max_elem_submatrix_row(matrix, i);

                if (double_funcs::equal(matrix[max_elem_row][i], 0)) return 0;
                else if (max_elem_row != i) {
                    matrix.swap_rows(max_elem_row, i);
                    det *= -1;
                }

                for (int row = i + 1; row < rows_; ++row)
                    detail::add_submatrix_koeff_row(matrix, row, i, -1 * matrix[row][i] / matrix[i][i]);
            }

            for (int i = 0; i < rows_; ++i)  det *= matrix[i][i];

            return det;
        }

    public:
        int ncols() const { return cols_; }
        int nrows() const { return rows_; }

        matrix_t& operator*=(T num) {
            for (int i = 0; i < rows_; ++i)
                for (int j = 0; j < cols_; ++j) data_[i][j] *= num;
            return *this;
        }

        matrix_t& operator+=(const matrix_t& matrix) {
            for (int i = 0; i < rows_; ++i)
                for (int j = 0; j < cols_; ++j) data_[i][j] += matrix.data_[i][j];
            return *this;
        }

        bool operator==(const matrix_t& matrix) const {
            for (int i = 0; i < rows_; ++i)
                for (int j = 0; j < cols_; ++j)
                    if (!double_funcs::equal(data_[i][j], matrix.data_[i][j])) return false;

            return true;
        }

        proxy_row_t operator[](int row) const { return proxy_row_t{data_[row]}; }

        void dump(std::ostream& os) const {
            for (int i = 0; i < rows_; ++i) {
                for (int j = 0; j < cols_; ++j) {
                    os << std::setw(9) << data_[i][j] << " ";
                }
                os << std::endl;
            }
            os << std::endl;
        }
    };

    template <typename T>
    matrix_t<T> operator*(const matrix_t<T>& matrix, T num) {
        matrix_t<T> tmp{matrix}; tmp *= num;
        return tmp;
    }

    template <typename T>
    matrix_t<T> operator*(T num, const matrix_t<T>& matrix) {
        matrix_t<T> tmp{matrix}; tmp *= num;
        return tmp;
    }

    template <typename T>
    matrix_t<T> operator+(const matrix_t<T>& mtx1, const matrix_t<T>& mtx2) {
        matrix_t<T> tmp{mtx1}; tmp += mtx2;
        return tmp;
    }

namespace detail {

    void add_submatrix_koeff_row(matrix_t<double>& matrix, int row, int add_row, double koeff) {
        for (int i = add_row + 1, iend = matrix.ncols(); i < iend; ++i)
            matrix[row][i] += koeff * matrix[add_row][i];

        matrix[row][add_row] = 0;
    }

    int find_max_elem_submatrix_row(matrix_t<double>& matrix, int start) {
        int i_max = start;
        double max_elem = std::abs(matrix[start][start]), curr = NAN;

        for (int i = start + 1, iend = matrix.nrows(); i < iend; ++i) {
            curr = std::abs(matrix[i][start]);
            if (max_elem < curr) {
                i_max = i;
                max_elem = curr;
            }
        }

        return i_max;
    }
}

}
