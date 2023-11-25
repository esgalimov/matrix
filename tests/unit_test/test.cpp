#include <gtest/gtest.h>
#include <vector>
#include "matrix.hpp"

using namespace matrix;


TEST(ctor_test, basic_ctor_test) {
    matrix_t<int> matrix{2, 3, 1};

    ASSERT_EQ(2, matrix.ncols());
    ASSERT_EQ(3, matrix.nrows());

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 2; ++j) {
            ASSERT_NO_THROW(matrix[i][j]);
            ASSERT_EQ(1, matrix[i][j]);
        }

    ASSERT_THROW(matrix[3][1], matrix_exceptions::MatrixOutOfRange);
    ASSERT_THROW(matrix[2][2], matrix_exceptions::MatrixOutOfRange);
}

TEST(ctor_test, zero_cols_rows_test) {
    ASSERT_THROW((matrix_t<int>{5, 0}), matrix_exceptions::MatrixZeroColsOrRows);
}

TEST(ctor_test, iterator_test) {
    std::vector<int> vec = {1, 2, 3, 4};

    matrix_t<int> matrix{3, 2, vec.begin(), vec.end()};

    ASSERT_EQ(3, matrix.ncols());
    ASSERT_EQ(2, matrix.nrows());

    ASSERT_EQ(1, matrix[0][0]);
    ASSERT_EQ(2, matrix[0][1]);
    ASSERT_EQ(3, matrix[0][2]);
    ASSERT_EQ(4, matrix[1][0]);
    ASSERT_EQ(0, matrix[1][1]);
    ASSERT_EQ(0, matrix[1][2]);
}

TEST(ctor_test, eye_test) {
    matrix_t<int> matrix = matrix_t<int>::eye(3);

    ASSERT_EQ(3, matrix.ncols());
    ASSERT_EQ(3, matrix.nrows());
    ASSERT_EQ(1, matrix[0][0]);
    ASSERT_EQ(1, matrix[1][1]);
    ASSERT_EQ(1, matrix[2][2]);
}

TEST(ctor_test, copy_ctor_test) {
    std::vector<int> vec = {1, 2, 3, 4};

    matrix_t<int>* matrix1 = new matrix_t<int>{2, 2, vec.begin(), vec.end()};
    matrix_t<int> matrix2{*matrix1};

    delete matrix1;

    ASSERT_EQ(2, matrix2.ncols());
    ASSERT_EQ(2, matrix2.nrows());

    ASSERT_EQ(1, matrix2[0][0]);
    ASSERT_EQ(2, matrix2[0][1]);
    ASSERT_EQ(3, matrix2[1][0]);
    ASSERT_EQ(4, matrix2[1][1]);
}

TEST(ctor_test, move_ctor_test) {
    std::vector<int> vec = {1, 2, 3, 4};

    matrix_t<int>* matrix1 = new matrix_t<int>{2, 2, vec.begin(), vec.end()};
    matrix_t<int> matrix2{std::move(*matrix1)};

    delete matrix1;

    ASSERT_EQ(2, matrix2.ncols());
    ASSERT_EQ(2, matrix2.nrows());

    ASSERT_EQ(1, matrix2[0][0]);
    ASSERT_EQ(2, matrix2[0][1]);
    ASSERT_EQ(3, matrix2[1][0]);
    ASSERT_EQ(4, matrix2[1][1]);
}

TEST(assig_test, copy_assig_test) {
    std::vector<int> vec = {1, 2, 3, 4};

    matrix_t<int>* matrix1 = new matrix_t<int>{2, 2, vec.begin(), vec.end()};
    matrix_t<int> matrix2;
    matrix2 = *matrix1;

    delete matrix1;

    ASSERT_EQ(2, matrix2.ncols());
    ASSERT_EQ(2, matrix2.nrows());

    ASSERT_EQ(1, matrix2[0][0]);
    ASSERT_EQ(2, matrix2[0][1]);
    ASSERT_EQ(3, matrix2[1][0]);
    ASSERT_EQ(4, matrix2[1][1]);
}

TEST(assig_test, move_assig_test) {
    std::vector<int> vec = {1, 2, 3, 4};

    matrix_t<int>* matrix1 = new matrix_t<int>{2, 2, vec.begin(), vec.end()};
    matrix_t<int> matrix2;
    matrix2 = std::move(*matrix1);

    delete matrix1;

    ASSERT_EQ(2, matrix2.ncols());
    ASSERT_EQ(2, matrix2.nrows());

    ASSERT_EQ(1, matrix2[0][0]);
    ASSERT_EQ(2, matrix2[0][1]);
    ASSERT_EQ(3, matrix2[1][0]);
    ASSERT_EQ(4, matrix2[1][1]);
}

TEST(operator_test, add_assig_test) {
    std::vector<int> vec = {1, 2, 3, 4};
    matrix_t<int> matrix1{2, 2, vec.begin(), vec.end()};
    matrix_t<int> matrix2{2, 2, 1};

    matrix1 += matrix2;

    ASSERT_EQ(2, matrix1[0][0]);
    ASSERT_EQ(3, matrix1[0][1]);
    ASSERT_EQ(4, matrix1[1][0]);
    ASSERT_EQ(5, matrix1[1][1]);
}

TEST(operator_test, add_test) {
    std::vector<int> vec = {1, 2, 3, 4};
    matrix_t<int> matrix1{2, 2, vec.begin(), vec.end()},
                  matrix2{2, 2, 1};
    matrix_t<int> matrix3 = matrix1 + matrix2;

    ASSERT_EQ(2, matrix3[0][0]);
    ASSERT_EQ(3, matrix3[0][1]);
    ASSERT_EQ(4, matrix3[1][0]);
    ASSERT_EQ(5, matrix3[1][1]);
}

TEST(operator_test, add_not_same_size_test) {
    std::vector<int> vec = {1, 2, 3, 4};
    matrix_t<int> matrix1{2, 2, vec.begin(), vec.end()},
                  matrix2{2, 3, 1};

    ASSERT_THROW(matrix1 + matrix2, matrix_exceptions::MatrixesAreNotSameSize);
}

TEST(operator_test, num_mult_test) {
    std::vector<int> vec = {1, 2, 3, 4, 5, 6};
    matrix_t<int> matrix1{2, 3, vec.begin(), vec.end()}, matrix2, matrix3;
    matrix2 = matrix1 * 5;
    matrix3 = 4 * matrix1;

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 2; ++j) {
            ASSERT_EQ(matrix1[i][j] * 5, matrix2[i][j]);
            ASSERT_EQ(matrix1[i][j] * 4, matrix3[i][j]);
    }
}

TEST(operator_test, multiply_test) {
    std::vector<int> vec1 = {1, 2, 3, 4, 5, 6},
                     vec2 = {1, -1, -1, 1};
    matrix_t<int> matrix1{2, 3, vec1.begin(), vec1.end()},
                  matrix2{2, 2, vec2.begin(), vec2.end()};

    ASSERT_NO_THROW(matrix1.multiply(matrix2));

    ASSERT_EQ(2, matrix1.ncols());
    ASSERT_EQ(3, matrix1.nrows());

    for (int i = 0; i < 3; ++i) ASSERT_EQ(-1, matrix1[i][0]);
    for (int i = 0; i < 3; ++i) ASSERT_EQ(1,  matrix1[i][1]);

    matrix2 = matrix_t<int>{3, 3};

    ASSERT_THROW(matrix1.multiply(matrix2), matrix_exceptions::MatrixesCannotBeMultiplied);
}

TEST(operator_test, equal_test) {
    std::vector<int> vec = {1, 2, 3, 4, 5, 6};
    matrix_t<int> matrix1{2, 3, vec.begin(), vec.end()}, matrix2 = matrix1;

    ASSERT_TRUE(matrix1 == matrix2);

    matrix1[1][1] = 10;
    ASSERT_FALSE(matrix1 == matrix2);

    matrix2 = matrix_t<int>::eye(2);
    ASSERT_FALSE(matrix1 == matrix2);
}

TEST(matrix_methods_test, tranpose_test) {
    std::vector<int> vec = {1, 2, 3, 4, 5, 6};
    matrix_t<int> matrix1{2, 3, vec.begin(), vec.end()}, matrix2 = matrix1;

    matrix2.transpose();

    ASSERT_EQ(3, matrix2.ncols());
    ASSERT_EQ(2, matrix2.nrows());

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 2; ++j)
            ASSERT_TRUE(matrix1[i][j] == matrix2[j][i]);
}

TEST(matrix_methods_test, trace_test) {
    std::vector<int> vec = {1, 2, 3, 4, 5, 6};
    matrix_t<int> matrix1{2, 3, vec.begin(), vec.end()};

    ASSERT_THROW(matrix1.trace(), matrix_exceptions::MatrixIsNotSquare);

    matrix1 = matrix_t<int>{2, 2, vec.begin(), vec.end()};

    ASSERT_EQ(5, matrix1.trace());
}

TEST(matrix_methods_test, negate_test) {
    std::vector<int> vec = {1, 2, 3, 4, 5, 6};
    matrix_t<int> matrix1{2, 3, vec.begin(), vec.end()}, matrix2 = matrix1;

    matrix2.negate();

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 2; ++j)
            ASSERT_EQ(-matrix1[i][j], matrix2[i][j]);
}

TEST(matrix_methods_test, det_test) {
    std::vector<int> vec = {1, 2, 3, 4, 5, 6};
    matrix_t<int> matrix1{2, 3, vec.begin(), vec.end()};

    ASSERT_THROW(matrix1.determinant(), matrix_exceptions::MatrixIsNotSquare);

    matrix1 = matrix_t<int>{2, 2, vec.begin(), vec.end()};
    ASSERT_EQ(-2, matrix1.determinant());

    matrix1.transpose();
    ASSERT_EQ(-2, matrix1.determinant());
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
