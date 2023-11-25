import random
import numpy as np

MATRIX_SIZE = 100

mtx = np.eye(MATRIX_SIZE, dtype=np.int32)

for i in range(6):
    pos = int(round(np.random.uniform(low=0, high=MATRIX_SIZE-1, size=None)))
    mtx[pos][pos] = 2

for i in range(3):
    pos = int(round(np.random.uniform(low=0, high=MATRIX_SIZE-1, size=None)))
    mtx[pos][pos] = 3



pos = int(round(np.random.uniform(low=0, high=MATRIX_SIZE-1, size=None)))
mtx[pos][pos] = -2

det = 1
for i in range(MATRIX_SIZE):
    det *= mtx[i][i]


def add_cols(m, n, koeff):
    global mtx
    for i in range(MATRIX_SIZE):
        mtx[i][m] += koeff * mtx[i][n]

def add_rows(m, n, koeff):
    global mtx
    mtx[m] += koeff * mtx[n]

def random_props():
    row1 = round(np.random.uniform(low=0, high=MATRIX_SIZE-1, size=None))
    row2 = round(np.random.uniform(low=0, high=MATRIX_SIZE-1, size=None))

    while row2 == row1:
        row2 = round(np.random.uniform(low=0, high=MATRIX_SIZE-1, size=None))

    koeff = round(np.random.uniform(low=-2, high=2, size=None))

    return row1, row2, koeff

for i in range(750):
    add_rows(*random_props())
    add_cols(*random_props())


with open("13.dat", "w") as mtx_file:
    print(MATRIX_SIZE, file=mtx_file)
    for i in range(MATRIX_SIZE):
        for j in range(MATRIX_SIZE):
            print(mtx[i][j], file=mtx_file, end="  ")
        print(file=mtx_file)

print(det)
#print(np.linalg.det(mtx))
