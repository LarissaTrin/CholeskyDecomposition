import sys
import os
import random
import numpy as np

n = 10000
m = 0
mat = np.zeros((n,n))

base_path = sys.path[0]

for i in range (0,n):
  if (i > 0): m=1
  for j in range(0,n):
    if(i==j):
      mat[i][j] = mat[i-m][i-m] + random.randint(5.0, 25.0) * 15
    else:
      mat[i][j] = random.randint(5.0, 25.0)
      mat[j][i] = random.randint(5.0, 25.0)

file = open(os.path.join(base_path, 'output_10000.txt'), 'w')
file.write(f'{n}\n')

for i in range(0, n):
    for j in range(0, n):
        file.write(f'{mat[i][j]} ')
    file.write(f'\n')

file.close()