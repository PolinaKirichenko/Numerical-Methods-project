import numpy as np

def integrate_simpson(tab_fun):
    n = tab_fun.shape[1]
    if n < 3:
        print("Error: Integral Simpson: too low grid size")
    if n % 2 == 0:
        print("Error: Integral Simpson: even number of intervals needed")
    if not np.all(np.diff(tab_fun[0], 2) < 1e-15):
        print("Error: Integral Simpson: grid must be uniform")
    N = n - 1

    h = tab_fun[0][1] - tab_fun[0][0]
    idx_4 = np.arange(1, N, 2)
    idx_2 = np.arange(2, N - 1, 2)
    res = h / 3 * (tab_fun[1][0] + 4 * tab_fun[1][idx_4].sum() + 2 * tab_fun[1][idx_2].sum() + tab_fun[1][N])

    return res

def integrate_trapez(tab_fun):
    n = tab_fun.shape[1]
    if n == 1:
        return 0
    res = 0
    h = tab_fun[0][1] - tab_fun[0][0]
    for i in range(n - 1):
        res += h / 2 * (tab_fun[1][i] + tab_fun[1][i + 1])
    return res