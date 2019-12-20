from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.animation import FuncAnimation
import numpy as np



def diffusion(u, dx, dy, dt, nt, nu):
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


        u[1: -1, 1: -1] = (un[1: -1, 1: -1] + (nu * dt / dx ** 2 * (un[1: -1, 2:] + un[1: -1, 0: -2] - 2 * un[1: -1, 1: -1]))
                                            + (nu * dt / dy ** 2 * (un[2:, 1: -1] + un[0: -2, 1: -1] - 2 * un[1: -1, 1: -1])))


        # set the boundary conditions
        u[0, :] = 1
        u[-1, :] = 1
        u[:, 0] = 1
        u[:, -1] = 1

    return u

if __name__ == "__main__":

    # defining the initial conditions
    nx = 31
    ny = 31
    nt = 17
    dx = 2 / (nx - 1)
    dy = 2 / (ny - 1)
    nu = .05
    sigma = .25
    dt = sigma * dx

    # define the grid points
    x = np.linspace(0, 2, nx)
    y = np.linspace(0, 2, ny)

    # define the initial wave profile u0, v0
    u0 = np.ones((nx, ny))
    u0[int(.5 / dy): int(1 / dy + 1), int(.5 / dx): int(1 / dx + 1)] = 2

    # initialise the figure and axis
    fig = plt.figure(figsize=(11, 7), dpi=100)
    ax = fig.gca(projection='3d')


    # define the grid and plot the u, v surfaces
    X, Y = np.meshgrid(x, y)

    def animate(nt, ax, X, Y, u0, dx, dy, dt, nu):

        ax.clear()
        ax.set_zlim(1, 2.5)
        ax.set_xlabel('$x$')
        ax.set_ylabel('$y$')

        u = diffusion(u0, dx, dy, dt, nt, nu)
        surface = ax.plot_surface(X, Y, u[:], cmap='viridis')

        return surface

    anim = FuncAnimation(fig, animate, fargs=(ax, X, Y, u0, dx, dy, dt, nu), interval=200, blit=False)
    anim.save('Diffusion-2d.gif', writer='imagemagick')

