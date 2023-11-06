#include <iostream>
#include "matrix.hpp"

using namespace matrix;

int main() {
    matrix_t mtx{10, 10, 1};
    //matrix_t mtx{10};
    mtx.dump(std::cout);
    return 0;
}
