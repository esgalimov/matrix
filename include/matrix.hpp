#pragma once

#include <iostream>
#include <utility>
#include <iomanip>
#include <type_traits>
#include "exceptions.hpp"
#include "double_funcs.hpp"


namespace matrix {
    template <typename T = int> class matrix_t;

    namespace detail {
        void add_submatrix_koeff_row(matrix_t<double>& matrix, size_t row, size_t add_row, double koeff = 1.0);
        size_t find_max_elem_submatrix_row(matrix_t<double>& matrix, size_t start);
    }


    template <typename T> class matrix_buffer_t {
    protected:
        T*  buffer_;
        T** data_;
        size_t cols_, rows_;

        matrix_buffer_t(const matrix_buffer_t&)            = delete;
        matrix_buffer_t& operator=(const matrix_buffer_t&) = delete;

        matrix_buffer_t(matrix_buffer_t&& rhs) noexcept :
            buffer_(rhs.buffer_), data_(rhs.data_),  cols_(rhs.cols_), rows_(rhs.rows_) {
                rhs.buffer_ = nullptr;
                rhs.data_ = nullptr;
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

            buffer_ = static_cast<T*> (::operator new(sizeof(T)  * cols * rows));
            data_   = static_cast<T**>(::operator new(sizeof(T*) * rows, std::nothrow));

            if (data_ == nullptr) {
                ::operator delete(buffer_);
                throw std::bad_alloc();
            }

            for (size_t i = 0; i < rows_; ++i)
                data_[i] = &buffer_[i * cols_];
        }

        ~matrix_buffer_t() {
            ::operator delete(buffer_);
            ::operator delete(data_);
        }
    };


    template <typename T> class matrix_t final : private matrix_buffer_t<T>  {
        static_assert(std::is_arithmetic<T>(), "Must be arithmetic type in matrix");

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
            for (size_t i = 0; i < rows; ++i)
                for (size_t j = 0; j < cols; ++j) data_[i][j] = val;
        }

        template<typename It>
        matrix_t(size_t cols, size_t rows, It start, It fin) : matrix_buffer_t<T>(cols, rows) {
            for (size_t i = 0; i < rows; ++i) {
                for (size_t j = 0; j < cols; ++j) {
                    if (start != fin) {
                        data_[i][j] = *start;
                        start++;
                    }
                    else data_[i][j] = T{};
                }
            }
        }

        matrix_t(const matrix_t& matrix) : matrix_buffer_t<T>(matrix.cols_, matrix.rows_) {
            for (size_t i = 0; i < rows_; ++i) {
                for (size_t j = 0; j < cols_; ++j) data_[i][j] = matrix.data_[i][j];
            }
        }

        template<typename U>
        matrix_t(const matrix_t<U>& matrix) : matrix_buffer_t<T>(matrix.ncols(), matrix.nrows()) {
            for (size_t i = 0; i < rows_; ++i)
                for (size_t j = 0; j < cols_; ++j)
                    data_[i][j] = static_cast<T>(matrix[i][j]);
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
            matrix_t tmp(*this);
            tmp *= -1;

            std::swap(*this, tmp);
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
                        elem += data_[i][k] * matrix[k][j];

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

        double gauss_algorithm() const {
            matrix_t<double> matrix{*this};

            double det = 1.0;

            for (size_t i = 0; i < cols_ - 1; ++i) {
                size_t max_elem_row = detail::find_max_elem_submatrix_row(matrix, i);

                if (double_funcs::equal(matrix[max_elem_row][i], 0)) return 0;
                else if (max_elem_row != i) {
                    matrix.swap_rows(max_elem_row, i);
                    det *= -1;
                }

                for (size_t row = i + 1; row < rows_; ++row)
                    detail::add_submatrix_koeff_row(matrix, row, i, -1 * matrix[row][i] / matrix[i][i]);
            }

            for (size_t i = 0; i < rows_; ++i)  det *= matrix[i][i];

            return det;
        }

    public:
        size_t ncols() const { return cols_; }
        size_t nrows() const { return rows_; }

        matrix_t& operator*=(T num) {
            for (size_t i = 0; i < rows_; ++i)
                for (size_t j = 0; j < cols_; ++j) data_[i][j] *= num;
            return *this;
        }

        matrix_t& operator+=(const matrix_t& matrix) {
            if (cols_ != matrix.ncols() || rows_ != matrix.nrows())
                throw matrix_exceptions::MatrixesAreNotSameSize();

            for (size_t i = 0; i < rows_; ++i)
                for (size_t j = 0; j < cols_; ++j) data_[i][j] += matrix.data_[i][j];
            return *this;
        }

        bool operator==(const matrix_t& matrix) const {
            if (cols_ != matrix.ncols() || rows_ != matrix.nrows()) return false;

            for (size_t i = 0; i < rows_; ++i)
                for (size_t j = 0; j < cols_; ++j)
                    if (!double_funcs::equal(data_[i][j], matrix.data_[i][j])) return false;

            return true;
        }

        proxy_row_t operator[](size_t row) const {
            if (row >= rows_) throw matrix_exceptions::MatrixOutOfRange();
            return proxy_row_t{data_[row], cols_};
        }

        void dump(std::ostream& os) const {
            for (size_t i = 0; i < rows_; ++i) {
                for (size_t j = 0; j < cols_; ++j) {
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
}
