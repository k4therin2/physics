"""metric pieces (g, ginv, gh, fields) shared by scripts; identical to einstein.py's construction."""
from tpoly import *
Psi = field_series('Psi'); Phi = field_series('Phi'); delta = field_series('delta')
om = [field_series('om%d' % i, orders=(2, 3)) for i in (1, 2, 3)]
v = [field_series('v%d' % i) for i in (1, 2, 3)]
hnames = {(0, 0): 'h11', (0, 1): 'h12', (0, 2): 'h13', (1, 1): 'h22', (1, 2): 'h23', (2, 2): 'h33'}
def hcomp(i, j):
    return field_series(hnames[(min(i, j), max(i, j))], orders=(2, 3))
gh = [[None] * 4 for _ in range(4)]
gh[0][0] = -(const(1) + Psi * 2)
for i in range(3):
    gh[0][i + 1] = gh[i + 1][0] = om[i]
    for j in range(3):
        gh[i + 1][j + 1] = ((const(1) - Phi * 2) if i == j else const(0)) + hcomp(i, j)
g = [[gh[a][b] * eta_pow(2) for b in range(4)] for a in range(4)]
eta0 = [[const(-1 if a == 0 else 1) if a == b else const(0) for b in range(4)] for a in range(4)]
H = [[gh[a][b] - eta0[a][b] for b in range(4)] for a in range(4)]
def matmul(A, B):
    return [[sum((A[a][c] * B[c][b] for c in range(4)), Poly({})) for b in range(4)] for a in range(4)]
M = matmul(H, eta0)
term = [[const(1) if a == b else const(0) for b in range(4)] for a in range(4)]
S = [[const(1) if a == b else const(0) for b in range(4)] for a in range(4)]
negM = [[-M[a][b] for b in range(4)] for a in range(4)]
for n in range(1, MAXORD + 1):
    term = matmul(term, negM)
    S = [[S[a][b] + term[a][b] for b in range(4)] for a in range(4)]
ghinv = matmul(eta0, S)
ginv = [[ghinv[a][b] * eta_pow(-2) for b in range(4)] for a in range(4)]
