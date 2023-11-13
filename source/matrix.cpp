#include <iostream>
#include <vector>
#include "matrix.hpp"

using namespace matrix;

int main() {
    //matrix_t mtx{10, 10, 1};
    //matrix_t mtx = matrix_t<>::eye(11);
    std::vector<int> vec = {1, 2, 3, 4, 5, 6, 56, 34, 23};
    matrix_t<int>* mtx = new matrix_t<int>{3, 4, vec.begin(), vec.end()};

    matrix_t mtxe;
    mtxe = *mtx + *mtx;
    mtxe = mtxe * 2;
    delete mtx;

    std::cout << mtxe[1][1] << std::endl;
    mtxe.dump(std::cout);

    return 0;
}
