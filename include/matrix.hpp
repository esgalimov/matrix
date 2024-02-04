#pragma once

#include <iostream>
#include <utility>
#include <iomanip>
#include <type_traits>
#include <memory>
#include <algorithm>
#include "exceptions.hpp"
#include "double_funcs.hpp"


namespace matrix {

    template <typename T> class matrix_buffer_t {
    protected:
        std::unique_ptr<T[]>  buffer_;
        std::unique_ptr<T*[]> data_;
        size_t cols_, rows_;

        matrix_buffer_t(const matrix_buffer_t&)            = delete;
        matrix_buffer_t& operator=(const matrix_buffer_t&) = delete;
        ~matrix_buffer_t()                                 = default;

        matrix_buffer_t(matrix_buffer_t&& rhs) noexcept :
            buffer_(std::move(rhs.buffer_)), data_(std::move(rhs.data_)), cols_(rhs.cols_), rows_(rhs.rows_) {
                rhs.cols_ = 0;
                rhs.rows_ = 0;
            }

        matrix_buffer_t& operator=(matrix_buffer_t&& rhs) noexcept {
            if (this == &rhs) return *this;

            std::swap(buffer_, rhs.buffer_);
            std::swap(data_, rhs.data_);
            std::swap(cols_, rhs.cols_);
            std::swap(rows_, rhs.rows_);

            return *this;
        }

        matrix_buffer_t(size_t cols = 0, size_t rows = 0) : cols_(cols), rows_(rows) {
            if (cols_ * rows_ == 0) throw matrix_exceptions::MatrixZeroColsOrRows();

            buffer_ = std::make_unique_for_overwrite<T[]>(cols * rows);
            data_   = std::make_unique_for_overwrite<T*[]>(rows);

            for (size_t i = 0; i < rows_; ++i)
                data_[i] = &buffer_[i * cols_];
        }
    };


    template <typename T>
        requires std::is_arithmetic_v<T>
    class matrix_t final : private matrix_buffer_t<T>  {
        using matrix_buffer_t<T>::buffer_;
        using matrix_buffer_t<T>::data_;
        using matrix_buffer_t<T>::cols_;
        using matrix_buffer_t<T>::rows_;

        class proxy_row_t {
            T* row_ = nullptr;
            size_t cols_;
            public:
                proxy_row_t(T* row, size_t cols) : row_(row), cols_(cols) {}

                const T& operator[](size_t col) const {
                    if (col >= cols_) throw matrix_exceptions::MatrixOutOfRange();
                    return row_[col];
                }
                T& operator[](size_t col) {
                    if (col >= cols_) throw matrix_exceptions::MatrixOutOfRange();
                    return row_[col];
                }
        };

    public:
        matrix_t(size_t cols = 1, size_t rows = 1, T val = T{}) : matrix_buffer_t<T>(cols, rows) {
            std::fill(buffer_.get(), buffer_.get() + rows_ * cols_, val);
        }

        template<typename It>
        matrix_t(size_t cols, size_t rows, It start, It fin) : matrix_buffer_t<T>(cols, rows) {
            if (std::distance(start, fin) != cols_ * rows_)
                throw matrix_exceptions::MatrixCtorBadElemCnt();

            std::copy(start, fin, buffer_.get());
        }

        matrix_t(const matrix_t& matrix) : matrix_buffer_t<T>(matrix.cols_, matrix.rows_) {
            std::copy(matrix.buffer_.get(), matrix.buffer_.get() + cols_ * rows_, buffer_.get());
        }

        template<typename U>
        matrix_t(const matrix_t<U>& matrix) : matrix_buffer_t<T>(matrix.ncols(), matrix.nrows()) {
            std::copy(&matrix[0][0], &matrix[rows_ - 1][cols_ - 1] + 1, buffer_.get());
        }

        matrix_t(matrix_t&& matrix)            = default;
        matrix_t& operator=(matrix_t&& matrix) = default;
        ~matrix_t()                            = default;


        matrix_t& operator=(const matrix_t& matrix) {
            matrix_t tmp(matrix);
            std::swap(*this, tmp);
            return *this;
        }


        static matrix_t eye(size_t n = 1) {
            matrix_t mtx{n, n, 0};

            for (size_t i = 0; i < n; ++i) mtx.data_[i][i] = 1;

            return mtx;
        }

        matrix_t& negate() & {
            std::for_each(buffer_.get(), buffer_.get() + cols_ * rows_, [](T& elem) { elem *= -1; });
            return *this;
        }

        matrix_t& transpose() & {
            matrix_t tmp(rows_, cols_);

            for (size_t i = 0; i < cols_; ++i)
                for (size_t j = 0; j < rows_; ++j)
                    tmp.data_[i][j] = data_[j][i];

            std::swap(*this, tmp);
            return *this;
        }

        matrix_t& multiply(const matrix_t& matrix) & {
            if (cols_ != matrix.nrows())
                throw matrix_exceptions::MatrixesCannotBeMultiplied();

            matrix_t tmp(matrix.ncols(), rows_);

            for (size_t i = 0, iend = tmp.nrows(); i < iend; ++i)
                for (size_t j = 0, jend = tmp.ncols(); j < jend; ++j) {
                    T elem = 0;

                    for (size_t k = 0; k < cols_; ++k)
                        elem += data_[i][k] * matrix.data_[k][j];

                    tmp[i][j] = elem;
                }

            std::swap(*this, tmp);

            return *this;
        }

        T trace() const {
            if (cols_ != rows_) throw matrix_exceptions::MatrixIsNotSquare();

            T ans = data_[0][0];
            for (size_t i = 1; i < cols_; ++i) ans += data_[i][i];

            return ans;
        }

        double determinant() const {
            if (cols_ != rows_) throw matrix_exceptions::MatrixIsNotSquare();

            return gauss_algorithm();
        }

        void swap_rows(size_t row1, size_t row2) {
            T* tmp = data_[row1];
            data_[row1] = data_[row2];
            data_[row2] = tmp;
        }

    private:

        static void add_submatrix_koeff_row(matrix_t<double>& matrix, size_t row, size_t add_row, double koeff) {
            for (size_t i = add_row + 1, iend = matrix.ncols(); i < iend; ++i)
                matrix[row][i] += koeff * matrix[add_row][i];

            matrix[row][add_row] = 0;
        }

        static size_t find_max_elem_submatrix_row(const matrix_t<double>& matrix, size_t start) {
            size_t i_max = start;
            double max_elem = std::abs(matrix[start][start]), curr = NAN;

            for (size_t i = start + 1, iend = matrix.nrows(); i < iend; ++i) {
                curr = std::abs(matrix[i][start]);
                if (max_elem < curr) {
                    i_max = i;
                    max_elem = curr;
                }
            }
            return i_max;
        }

        double gauss_algorithm() const {
            matrix_t<double> matrix{*this};

            double det = 1.0;

            for (size_t i = 0; i < cols_ - 1; ++i) {
                size_t max_elem_row = find_max_elem_submatrix_row(matrix, i);

                if (double_funcs::equal(matrix[max_elem_row][i], 0)) return 0;
                else if (max_elem_row != i) {
                    matrix.swap_rows(max_elem_row, i);
                    det *= -1;
                }

                for (size_t row = i + 1; row < rows_; ++row)
                    add_submatrix_koeff_row(matrix, row, i, -1 * matrix[row][i] / matrix[i][i]);

                det *= matrix[i][i];
            }

            return det * matrix[rows_ - 1][cols_ - 1];
        }

    public:
        size_t ncols() const { return cols_; }
        size_t nrows() const { return rows_; }

        matrix_t& operator*=(T num) {
            for (auto start = buffer_.get(), fin = buffer_.get() + cols_ * rows_; start < fin; ++start)
                *start *= num;
            return *this;
        }

        matrix_t& operator+=(const matrix_t& matrix) {
            if (cols_ != matrix.ncols() || rows_ != matrix.nrows())
                throw matrix_exceptions::MatrixesAreNotSameSize();

            for (auto start = buffer_.get(), fin = buffer_.get() + cols_ * rows_,
                      src = matrix.buffer_.get(); start < fin; ++start, ++src) {
                *start += *src;
            }
            return *this;
        }

        bool operator==(const matrix_t& matrix) const {
            if (cols_ != matrix.ncols() || rows_ != matrix.nrows()) return false;

            for (auto start = buffer_.get(), fin = buffer_.get() + cols_ * rows_,
                      src = matrix.buffer_.get(); start < fin; ++start, ++src) {
                if (!double_funcs::equal(*start, *src)) return false;
            }

            return true;
        }

        proxy_row_t operator[](size_t row) {
            if (row >= rows_) throw matrix_exceptions::MatrixOutOfRange();
            return proxy_row_t{data_[row], cols_};
        }

        const proxy_row_t operator[](size_t row) const {
            if (row >= rows_) throw matrix_exceptions::MatrixOutOfRange();
            return proxy_row_t{data_[row], cols_};
        }

        void dump(std::ostream& os) const {
        #ifndef NDEBUG
            for (size_t i = 0; i < rows_; ++i) {
                for (size_t j = 0; j < cols_; ++j) {
                    os << std::setw(9) << data_[i][j] << " ";
                }
                os << std::endl;
            }
            os << std::endl;
        #endif
        }
    };

    template <typename T>
    matrix_t<T> operator*(const matrix_t<T>& matrix, T num) {
        matrix_t<T> tmp{matrix}; tmp *= num;
        return tmp;
    }

    template <typename T>
    matrix_t<T> operator*(matrix_t<T>&& matrix, T num) {
        matrix *= num;
        return matrix;
    }

    template <typename T>
    matrix_t<T> operator*(T num, const matrix_t<T>& matrix) {
        matrix_t<T> tmp{matrix}; tmp *= num;
        return tmp;
    }

    template <typename T>
    matrix_t<T> operator*(T num, matrix_t<T>&& matrix) {
        matrix *= num;
        return matrix;
    }

    template <typename T>
    matrix_t<T> operator+(const matrix_t<T>& mtx1, const matrix_t<T>& mtx2) {
        matrix_t<T> tmp{mtx1}; tmp += mtx2;
        return tmp;
    }

    template <typename T>
    matrix_t<T> operator+(matrix_t<T>&& mtx1, const matrix_t<T>& mtx2) {
        mtx1 += mtx2; return mtx1;
    }

    template <typename T>
    matrix_t<T> operator+(const matrix_t<T>& mtx1, matrix_t<T>&& mtx2) {
        mtx2 += mtx1; return mtx2;
    }

    template <typename T>
    matrix_t<T> operator+(matrix_t<T>&& mtx1, matrix_t<T>&& mtx2) {
        mtx1 += mtx2; return mtx1;
    }
}
