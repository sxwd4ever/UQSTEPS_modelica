'''
The plot lib from hal's steps
'''

import matplotlib.image as mpimg
def imgplot(X, Y, XXYY, cs, imgfile="", ax=[], xticks=[], yticks=[]):
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
        y=y2-(Y2-Y)/(Y2-Y1)*(y2-y1)
        return (x,y)
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