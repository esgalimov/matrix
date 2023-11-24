#pragma once

#include <exception>

namespace matrix_exceptions {
    struct MatrixIsNotSquare : public std::runtime_error {
        MatrixIsNotSquare() : std::runtime_error("Matrix must be square to call this method") {}
    };

    struct MatrixesAreNotSameSize : public std::runtime_error {
        MatrixesAreNotSameSize() : std::runtime_error("Matrixes must be same size") {}
    };

    struct MatrixZeroColsOrRows : public std::runtime_error {
        MatrixZeroColsOrRows() : std::runtime_error("Matrix has cols_ = 0 or rows_ = 0") {}
    };

    struct MatrixOutOfRange : public std::out_of_range {
        MatrixOutOfRange() : std::out_of_range("Matrix cols_ or rows_ out of range") {}
    };
}
