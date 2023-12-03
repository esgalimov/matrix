#include <iostream>
#include <vector>
#include "matrix.hpp"

using namespace matrix;

int main() {
    size_t mtx_size = 0;

    if (!(std::cin >> mtx_size)) {
        std::cerr << "Bad matrix size" << std::endl;
        return 1;
    }

    std::vector<double> elem_vec;
    elem_vec.reserve(mtx_size * mtx_size);

    double curr = NAN;

    for (int i = 0, iend = mtx_size * mtx_size; i < iend; ++i) {
        if (!(std::cin >> curr)) {
            std::cerr << "Bad matrix elem on pos = " << i << std::endl;
            return 1;
        }
        elem_vec.push_back(curr);
    }

    matrix_t<double> mtx{mtx_size, mtx_size, elem_vec.begin(), elem_vec.end()};

    mtx.dump(std::cout);

    try {
        std::cout << mtx.determinant() << std::endl;
    }
    catch (std::bad_alloc&) {
        std::cerr << "Allocation failed" << std::endl;
    }
    catch (const std::runtime_error& mtx_exc) {
        std::cerr << mtx_exc.what() << std::endl;
    }

    return 0;
}
