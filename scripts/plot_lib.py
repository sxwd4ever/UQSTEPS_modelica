'''
The plot lib from hal's steps
'''

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np

def openplot(figno):
    plt.close(1)
    f1=plt.figure()
    ax = f1.add_subplot(111)
    return (f1, ax)

def plotanno(ax,xlabel="", ylabel="", xlim=[], ylim=[], grid="on", \
                legendloc="", title=""):
    if ylim!=[]:
        ax.set_ylim(ylim)
    if xlim!=[]:
        ax.set_xlim(xlim)
    if xlabel!="":
        ax.set_xlabel(xlabel)
    if ylabel!="":
        ax.set_ylabel(ylabel)
    if grid!="":
        if grid=="on" or grid=="off":
            ax.grid(grid)
        else:
            ax.grid(True, grid)
    if legendloc!="":
        ax.legend(loc=legendloc)
    if title!="":
        ax.set_title(title)

def imgplot(X, Y, XXYY, cs, imgfile="", ax=[], xticks=[], yticks=[]):

    if ax==[]:
        img=mpimg.imread(imgfile)
        f1=plt.figure()
        ax = f1.add_subplot(111)
        ax.imshow(img)
    else:
        f1=plt.gcf()
    (y1,y2)=ax.get_ylim()
    (x1,x2)=ax.get_xlim()
    xxyy=[(x1, x2),(y1,y2)]
    (x,y)=ptransform((X,Y), XXYY, xxyy)
    ax.plot(x,y,cs)
#    ax.yaxis.set_major_locator(plt.NullLocator())
#    ax.xaxis.set_major_formatter(plt.NullFormatter())
    if xticks!=[]:
        n=len(xticks)
        xtickpos=np.zeros(n)
        dx=(x2-x1)/(n-1)
        for k in range(0,n):
            xtickpos[k]=k*dx
#        ax.set_xticks([x1, x2])
        ax.set_xticks(xtickpos)
        ax.set_xticklabels(xticks)
#        print("(x2-x1)/(n-1)=(%f-%f)/%d=%f XTICKPOS="%(x2,x1,n-1,dx), end="")
        # print(xtickpos)
    if yticks!=[]:
        n=len(yticks)
        ytickpos=np.zeros(n)
        dy=abs(y2-y1)/(n-1)
        for k in range(0,n):
            ytickpos[k]=y1-k*dy
        ax.set_yticks(ytickpos)
        ax.set_yticklabels(yticks)
        # print("(y2-y1)/(n-1)=(%f-%f)/%d=%f yTICKPOS="%(y2,y1,n-1,dy), end="")
        # print(ytickpos)
    return (f1, ax)

def ptransform(XY, XXYY, xxyy):
    X1=XXYY[0][0]
    X2=XXYY[0][1]
    Y1=XXYY[1][0]
    Y2=XXYY[1][1]
    x1=xxyy[0][0]
    x2=xxyy[0][1]
    y1=xxyy[1][0]
    y2=xxyy[1][1]
    X=XY[0]
    Y=XY[1]
    x=(X-X1)/(X2-X1)*(x2-x1)+x1
    # y=y2-(Y2-Y)/(Y2-Y1)*(y2-y1)
    y = y2 + (y1 - y2) * (Y - Y1) / (Y2 - Y1)
    return (x,y)

def cal_rel_axes(x, x_range):
    x_min = min(x_range)
    x_max = max(x_range)
    return (x - x_min) / ( x_max - x_min)

def plot_xy_series(x, y, range_x, range_y, img_file = "", color_style = 'r*', axes_cur = []):
    ''' 
    plot a series of (x,y) on an figure
    '''  
    if axes_cur==[]:
        img = mpimg.imread(img_file)
        # flip the image upside down and with set orign ='lower' when call imshow
        # to make the coordinates works as 'normal' way
        img = np.flipud(img)
        f1 = plt.figure()
        axes_cur = f1.add_subplot(111)
        axes_cur.imshow(img, origin='lower')
        
        # (y_img_1, y_img_2) = axes_cur.get_ylim()
        # (x_img_1, x_img_2) = axes_cur.get_xlim()
        #plt.show()
    else:
        f1 = plt.gcf()

    # calculate x, y in axes coordinates
    x_axes = cal_rel_axes(x, range_x)
    y_axes = cal_rel_axes(y, range_y)

    # combine it into array of turple
    xy_axes = [(x_axes[i], y_axes[i]) for i in range(0, len(x_axes))]
    # transfrom axes coordinates to display coordinates
    xy_disp = axes_cur.transAxes.transform(xy_axes)
    # transfrom from display coordinates to data coordinates
    xy_data = axes_cur.transData.inverted().transform(xy_disp)
    # unpack into seperate x, y series 
    x_data = [xy_data[i][0] for i in range(0, len(xy_data))]
    y_data = [xy_data[i][1] for i in range(0, len(xy_data))]
    
    # draw data
    axes_cur.plot(x_data, y_data, color_style)

    return (f1, axes_cur)
