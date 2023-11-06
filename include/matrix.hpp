#pragma once

#include <iostream>

namespace matrix {
    template <typename T = int> class matrix_t {
        T **data;
        T *buffer;
        int cols_, rows_;

        public:
            matrix_t(int cols = 1, int rows = 1, T val = T{}) : cols_(cols), rows_(rows) {
                data = new int*[rows];
                buffer = new int[rows * cols];

                for (int i = 0; i < rows; ++i) {
                    data[i] = &buffer[i * cols];

                    for (int j = 0; j < cols; ++j) data[i][j] = val;
                }
            }
            // add ctor from sequence
            template<typename It>
            matrix_t(int cols, int rows, It start, It fin) {
                data = new int*[rows];
                buffer = new int[rows * cols];

                for (int i = 0; i < rows; ++i) {
                    data[i] = &buffer[i * cols];

                    for (int j = 0; j < cols; ++j) {
                        if (start != fin) {
                            data[i][j] = *start;
                            start++;
                        }
                        else data[i][j] = 0;
                    }
                }
            }

            static matrix_t eye(int n = 1) {
                matrix_t mtx{n, n, 0};

                for (int i = 0; i < n; ++i) {
                    mtx.data[i][i] = 1;
                }
                return mtx;
            }

            void dump(std::ostream& os) const {
                for (int i = 0; i < rows_; ++i) {
                    for (int j = 0; j < cols_; ++j) {
                        os << data[i][j] << " ";
                    }
                    os << std::endl;
                }
            }
    };
}
