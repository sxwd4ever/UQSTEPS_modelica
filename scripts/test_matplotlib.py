import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def show_xy_lim():

    x = np.arange(0, 10, 0.005)
    y = np.exp(-x/2.) * np.sin(2*np.pi*x)

    x_lim = [(0, 10), (10, 20)]
    y_lim = [(-1, 1), (-1, 2)]

    fig, ax = plt.subplots()
    ax.plot(x, y)
    ax.set_xlim(x_lim[0])
    ax.set_ylim(y_lim[0])

    points =[(0, 0), (5, 0), (0, -1), (10, 1)]

    print('see the value of "display" coordinates for points in the "data" coordinates:')

    for p in points:
        print("{0}_data = {1}_disp ".format(p, ax.transData.transform(p)))


    print('now see the change of "display" coordinates if we change "the data" coordinates:')

    p = points[0]

    print('before set, ylim = {0}'.format(y_lim[0]))

    print("{0}_data = {1}_disp ".format(p, ax.transData.transform(p)))

    ax.set_ylim(y_lim[1])

    print('after set, ylim = {0} (reset afterwards)'.format(y_lim[1]))

    print("{0}_data = {1}_disp ".format(p, ax.transData.transform(p)))

    ax.set_ylim(y_lim[0])

    print('before set, xlim = {0}'.format(x_lim[0]))

    print("{0}_data = {1}_disp ".format(p, ax.transData.transform(p)))

    ax.set_xlim(x_lim[1])

    print('after set, xlim = {0} (reset afterwards)'.format(x_lim[1]))

    print("{0}_data = {1}_disp ".format(p, ax.transData.transform(p)))

    ax.set_xlim(x_lim[0])

def physical_cor():
    fig, ax = plt.subplots(figsize=(5, 4))
    x, y = 10*np.random.rand(2, 1000)
    ax.plot(x, y*10., 'go', alpha=0.2)  # plot some data in data coordinates
    # add a circle in fixed-coordinates
    circ = mpatches.Circle((2.5, 2), 1.0, transform=fig.dpi_scale_trans,
                        facecolor='blue', alpha=0.75)
    ax.add_patch(circ)
    plt.show()

def main():
    physical_cor()

if __name__== '__main__' :
    main()
