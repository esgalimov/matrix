#include "matrix.hpp"

namespace matrix {
    namespace detail {

    void add_submatrix_koeff_row(matrix_t<double>& matrix, size_t row, size_t add_row, double koeff) {
        for (size_t i = add_row + 1, iend = matrix.ncols(); i < iend; ++i)
            matrix[row][i] += koeff * matrix[add_row][i];

        matrix[row][add_row] = 0;
    }

    size_t find_max_elem_submatrix_row(matrix_t<double>& matrix, size_t start) {
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
    }
}
