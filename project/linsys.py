import numpy as np

def solve_upper(A, b):
    n = A.shape[0]
    res = np.zeros(n)

    for i in range(n - 1, -1, -1):
        res[i] = b[i]
        for j in range(i + 1, n):
            res[i] -= A[i][j] * res[j]
        res[i] /= A[i][i]

    return res


def solve_lower(M_, b_):
    M = M_[::-1][:, ::-1]
    b = b_[::-1]
    res_upper = solve_upper(M, b);
    return res_upper[::-1]


def lu(A, b):
    n = A.shape[0]
    L = np.zeros((n, n))
    U = np.zeros((n, n))
    
    for k in range(n):
        for i in range(k, n):
            L[i][k] = A[i][k]
            for m in range(k):
                L[i][k] -= L[i][m] * U[m][k]
        
        U[k][k] = 1
        if (L[k][k] == 0):
            print("LU error: singularity")
            return None
        
        for j in range(k + 1, n):
            U[k][j] = A[k][j]
            for m in range(k):
                U[k][j] -= L[k][m] * U[m][j]
            U[k][j] /= L[k][k]

    return solve_upper(U, solve_lower(L, b))
