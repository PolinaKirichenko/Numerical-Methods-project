import numpy as np 

def runge_kutta(xn, n, y0, z0, f, g):
    """
    System of differential equations in the form:
    dy/dx = f(x, y, z)
    dz/dx = g(x, y, z)
    """
    h = xn / (n-1)
    x = [0]; y = [y0]; z = [z0]
    for i in range(1, n):
        k = np.zeros(4); l = np.zeros(4)
        k[0] = h * f(x[-1], y[-1], z[-1])
        l[0] = h * g(x[-1], y[-1], z[-1])
        
        k[1] = h * f(x[-1] + h / 2, y[-1] + k[0] / 2, z[-1] + l[0] / 2)
        l[1] = h * g(x[-1] + h / 2, y[-1] + k[0] / 2, z[-1] + l[0] / 2)
        
        k[2] = h * f(x[-1] + h / 2, y[-1] + k[1] / 2, z[-1] + l[1] / 2)
        l[2] = h * g(x[-1] + h / 2, y[-1] + k[1] / 2, z[-1] + l[1] / 2)
        
        k[3] = h * f(x[-1] + h, y[-1] + k[2], z[-1] + l[2])
        l[3] = h * g(x[-1] + h, y[-1] + k[2], z[-1] + l[2])
        
        y.append(y[-1] + np.sum(k * np.array([1, 2, 2, 1])) / 6)
        z.append(z[-1] + np.sum(l * np.array([1, 2, 2, 1])) / 6)
        x.append(i * h)
    return np.array(y), np.array(z)
