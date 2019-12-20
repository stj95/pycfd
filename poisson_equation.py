import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D


def plot2D(x, y, p):

    fig = plt.figure(figsize=(11, 7), dpi=100)
    ax = fig.gca(projection='3d')
    X, Y = np.meshgrid(x, y)
    surface = ax.plot_surface(X, Y, p, rstride=1, cstride=1, cmap=cm.viridis, linewidth=0, antialiased=False)
    ax.set_xlim(0, 2)
    ax.set_ylim(0, 1)
    ax.view_init(30, 225)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')

    return surface

def poisson2d(p, b, nx, ny, dx, dy, nt):

    n = 1

    while n < nt:

        pn = p.copy()

        p[1:-1, 1:-1] = ((dy**2 * (pn[1:-1, 2:] + pn[1:-1, :-2]) +
                          dx**2 * (pn[2:, 1:-1] + pn[:-2, 1:-1]) -
                          dx**2 * dy**2 * b[1:-1, 1:-1]) /
                         (2*(dx**2 + dy**2)))


        # set boundary conditions
        p[:, 0] = 0
        p[ny-1, :] = 0
        p[:, 0] = 0
        p[:, nx-1] = 0

        n += 1

    return p


if __name__ == "__main__":

    # declare variables
    nx = 31
    ny = 31
    nt = 100
    dx = 2 / (nx - 1)
    dy = 1 / (ny - 1)

    # IC's:
    p = np.zeros((ny, nx))
    b = np.zeros((ny, nx))
    b[int(ny / 4), int(nx / 4)] = 100
    b[int(3 * ny / 4), int(3 * nx / 4)] = -100

    # set dimensions:
    x = np.linspace(0, 2, nx)
    y = np.linspace(0, 1, ny)

    p_conv = poisson2d(p, b, nx, ny, dx, dy, 10)
    plot2D(x, y, p_conv)

    plt.show()