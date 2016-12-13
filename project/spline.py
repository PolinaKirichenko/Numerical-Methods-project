import numpy as np 

def tridiagonal(A, C, B, F):
    k = A.size
    
    alpha = np.zeros(k)
    beta = np.zeros(k)

    for i in range(1, k):
        denom = (A[i - 1] * alpha[i - 1] + C[i - 1])
        alpha[i] = -B[i - 1] / denom
        beta[i] = (F[i - 1] - A[i - 1] * beta[i - 1]) / denom

    res = np.zeros(k)
    res[k - 1] = (F[k - 1] - A[k - 1] * beta[k - 1]) / (C[k - 1] + A[k - 1] * alpha[k - 1])

    for i in range(k - 2, -1, -1):
        res[i] = alpha[i + 1] * res[i + 1] + beta[i + 1]

    return res

class CubicSpline:
    def __init__(self, x, y): 
        print("CubicSpline edge", x[0], x[-1])
        self.x = x
        self.y = y
        n = x.size
        self.splines = np.zeros((n - 1, 4))

        for i in range(n - 1):
            self.splines[i][0] = self.y[i + 1]

        A = np.zeros((n - 2))
        C = np.zeros((n - 2))
        B = np.zeros((n - 2))
        F = np.zeros((n - 2))

        # A[0] и B[n-2] не будут использованы #
        for i in range(1, n - 1):
            h_i = x[i] - x[i - 1]
            h_i1 = x[i + 1] - x[i];
            A[i - 1] = h_i;
            C[i - 1] = 2 * (h_i + h_i1)
            B[i - 1] = h_i1
            F[i - 1] = 6 * ((self.y[i + 1] - self.y[i]) / h_i1 - (self.y[i] - self.y[i - 1]) / h_i);

        self.splines[:, 2] = np.append(tridiagonal(A, C, B, F), 0)

        for i in range(n - 2, -1, -1):
            h_i = self.x[i + 1] - self.x[i]
            self.splines[i][1] = (self.y[i + 1] - self.y[i]) / h_i + \
                                  h_i * (2 * self.splines[i][2] + self.splines[i - 1][2]) / 6
            self.splines[i][3] = (self.splines[i][2] - self.splines[i - 1][2]) / h_i

        
    def at(self, x):
        idx = np.searchsorted(self.x, x)
        if (idx == 0 or idx == self.x.size) and (x != self.x[0] and x != self.x[-1]):
            print("Spline error: point", x, "out of interval", self.x[0], self.x[-1])
            return None
        z = x - self.x[idx]
        
        if idx == 0:
            idx += 1
        a = self.splines[idx-1][0]
        b = self.splines[idx-1][1]
        c = self.splines[idx-1][2]
        d = self.splines[idx-1][3]
        
        return a + (b + (c / 2 + d * z / 6) * z) * z


    def grad(self, x):
        idx = np.searchsorted(self.x, x)
        if (idx == 0 or idx == self.x.size) and (x != self.x[0] and x != self.x[-1]):
            print("Spline grad error: point", x, "out of interval", self.x[0], self.x[-1])
            return None
        z = x - self.x[idx]
        
        if idx == 0:
            idx += 1
        b = self.splines[idx-1][1]
        c = self.splines[idx-1][2]
        d = self.splines[idx-1][3]
        
        return b + 2 * c * z + 3 * d * z**2
