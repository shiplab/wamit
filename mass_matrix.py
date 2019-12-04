import numpy as np
import sys

rho = float(sys.argv[1])
Vol = float(sys.argv[2])
L = float(sys.argv[3])
B = float(sys.argv[4])

# rho = 1.025
# Vol = 10000
# L = 100
# B = 20

print('rho = ' + str(rho))
print('Vol = ' + str(Vol))
print('L = ' + str(L))
print('B = ' + str(B))
print('')

M = rho*Vol
MM = np.diag([M, M, M, M*(.35*B)**2, M*(.25*L)**2, M*(.25*L)**2] )
MM.shape

np.savetxt(sys.stdout.buffer, MM, fmt="%.5E")
