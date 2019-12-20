from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np



def non_linear_convection(u, v, dx, dy, dt, nt, c):
    """
    discritizes the non-linear convection equation allowing the change of some initial profile u0 with time

    :param u0: Initial profile (also defines the dimensions of the domain)
    :param dx: distance between x values (defines the resolution of the domain)
    :param ny: distance between y values (defines the resolution of the domain)
    :param dt: amount of time between each step (defines how far ahead with the wave propagation you want to look)
    :param nt: the number of time steps
    :param c:  wave speed
    :return:
    """

    for n in range(nt + 1):

        un = u.copy()
        vn = v.copy()

        u[1:, 1:] = un[1:, 1:] - (un[1:, 1:] * (un[1:, 1:] - un[1:, :-1]) * c * dt / dx) \
                               - (vn[1:, 1:] * (un[1:, 1:] - un[:-1, 1:]) * c * dt / dy)

        v[1:, 1:] = vn[1:, 1:] - (un[1:, 1:] * (vn[1:, 1:] - vn[1:, :-1]) * c * dt / dx) \
                               - (vn[1:, 1:] * (vn[1:, 1:] - vn[:-1, 1:]) * c * dt / dy)

        # set the boundary conditions
        u[0, :] = 1
        u[-1, :] = 1
        u[:, 0] = 1
        u[:, -1] = 1

        v[0, :] = 1
        v[-1, :] = 1
        v[:, 0] = 1
        v[:, -1] = 1

    return u, v

if __name__ == "__main__":

    # defining the initial conditions
    nx = 101
    ny = 101
    nt = 80
    dx = 2 / (nx - 1)
    dy = 2 / (ny - 1)
    c = 1
    sigma = .2
    dt = sigma * dx

    # define the grid points
    x = np.linspace(0, 2, nx)
    y = np.linspace(0, 2, ny)

    # define the initial wave profile u0, v0
    u0 = np.ones((nx, ny))
    v0 = np.ones((nx, ny))
    u0[int(.5 / dy): int(1 / dy + 1), int(.5 / dx): int(1 / dx + 1)] = 2
    v0[int(.5 / dy): int(1 / dy + 1), int(.5 / dx): int(1 / dx + 1)] = 2

    # find the solution after nt time steps
    u, v = non_linear_convection(u0, v0, dx, dy, dt, nt, c)

    # initialise the figure and axis
    fig = plt.figure(figsize=(10, 3.5), dpi=100)
    ax1 = fig.add_subplot(1, 2, 1, projection='3d')
    ax2 = fig.add_subplot(1, 2, 2, projection='3d')

    # define the grid and plot the u, v surfaces
    X, Y = np.meshgrid(x, y)
    ax1.plot_surface(X, Y, u, cmap='viridis', rstride=2, cstride=2)
    ax2.plot_surface(X, Y, v, cmap='viridis', rstride=2, cstride=2)

    # set labels
    ax1.set_xlabel('$x$')
    ax1.set_ylabel('$y$')
    ax2.set_xlabel('$x$')
    ax2.set_ylabel('$y$')

    plt.show()


