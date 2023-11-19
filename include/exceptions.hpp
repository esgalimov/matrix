#pragma once

#include <exception>

namespace matrix_exceptions {
    struct MatrixIsNotSquare : public std::exception {
        const char* what() const noexcept override {
            return "Matrix must be square to call this method";
        }
    };

    struct MatrixesAreNotSameSize : public std::exception {
        const char* what() const noexcept override {
            return "Matrixes must be same size";
        }
    };
}
