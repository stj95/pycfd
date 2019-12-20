import numpy as np
import matplotlib.pyplot as plt
import time, sys


if __name__ == "__main__":


    """
    define grid
    """
    nx = 81
    dx = 2 / (nx - 1)
    nt = 25
    dt = 0.025
    c = 1

    """
    Initial velocity profile
    """
    u = np.ones(nx)
    u[int(.5 / dx): int(1 / dx + 1)] = 2

    """
    Finite difference
    """

    for n in range(nt):
        un = u.copy()
        for i in range(1, nx):
            u[i] = un[i] - c * dt / dx * (un[i] - un[i-1])

        plt.plot(np.linspace(0, 2, nx), u)

    plt.show()

    # change in shape caused by error in finite difference