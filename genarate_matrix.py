import sys
import os
import random
import numpy as np

base_path = sys.path[0]

n = 5 #int(input("Digite a dimensão da matriz: "))


file = open(os.path.join(base_path, 'output.txt'), 'w')
file.write(f'{n}\n')

matrix = []
for i in range(0, n):
    line = []
    for j in range(0, n):
        line.append(0)
    matrix.append(line)
det = -1
while det < 1:
    for i in range(0, n):
        for j in range(0, n):
            matrix[i][j] == 0

    for i in range(0, n):
        for j in range(0, n):
            number_str = float(random.random())
            number_str = round(number_str, 2)
            if matrix[i][j] == 0:
                matrix[i][j] = number_str
                matrix[j][i] = number_str
    for i in range(0, n):
        for j in range(0, n):
            number_str = float(random.randint(5, 15))
            matrix[i][i] = number_str

    matriz = np.array(matrix)
    det = np.linalg.det(matriz)

for i in range(0, n):
    for j in range(0, n):
        file.write(f'{matrix[i][j]} ')
    file.write(f'\n')

file.close()