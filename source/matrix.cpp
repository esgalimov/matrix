#include <iostream>
#include <vector>
#include "matrix.hpp"

using namespace matrix;

int main() {
    int mtx_size = 0;

    if (!(std::cin >> mtx_size)) {
        std::cerr << "Bad matrix size" << std::endl;
        return 1;
    }

    std::vector<double> elem_vec;
    elem_vec.reserve(mtx_size * mtx_size);

    double curr = NAN;

    for (int i = 0, iend = mtx_size * mtx_size; i < iend; ++i) {
        if (!(std::cin >> curr)) {
            std::cerr << "Bad matrix elem" << std::endl;
            return 1;
        }
        elem_vec.push_back(curr);
    }

    matrix_t<double> mtx{mtx_size, mtx_size, elem_vec.begin(), elem_vec.end()};

    mtx.dump(std::cout);

    std::cout << mtx.determinant() << std::endl;

    return 0;
}
