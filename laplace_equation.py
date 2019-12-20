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

def laplace2d(p, bc, dx, dy, target):

    norm = 1

    while norm > target:

        pn = p.copy()

        p[1:-1, 1:-1] = ((dy**2 * (pn[1:-1, 2:] + pn[1:-1, :-2]) +
                          dx**2 * (pn[2:, 1:-1] + pn[:-2, 1:-1])) /
                         (2*(dx**2 + dy**2)))


        # set boundary conditions
        p[:, 0] = 0
        p[:, -1] = bc
        p[0, :] = p[1, :]
        p[-1, :] = p[-2, :]

        norm = np.sum(np.abs(p[:]) - np.abs(pn[:])) / np.sum(np.abs(pn[:]))

    return p


if __name__ == "__main__":

    # declare variables
    nx = 31
    ny = 31
    dx = 2 / (nx - 1)
    dy = 2 / (ny - 1)

    # IC's:
    p = np.zeros((ny, nx))

    # set dimensions:
    x = np.linspace(0, 2, nx)
    y = np.linspace(0, 1, ny)

    p_conv = laplace2d(p, y, dx, dy, 1e-4)
    plot2D(x, y, p_conv)

    plt.show()