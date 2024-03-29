'''
The plot lib from hal's steps
'''

import matplotlib.image as mpimg
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np

from enum import IntEnum
from matplotlib.figure import Figure
from utils import mkdir_filepath

class AxisType(IntEnum):
    Primary = 0,
    Secondary = 1
class DataSeries(object):
    ''' 
    2D data Series used to contain a series of (x,y) data
    '''
    def __init__(self, name, x:dict, y:dict, cs = 'r-s', ax_type=AxisType.Primary):
        self.x = x['value']
        self.y = y['value']
        self.name = name
        self.range_x = x['axis_range']
        self.range_y = y['axis_range']
        self.cs = cs
        self.ax_type = ax_type
        self.ticks = [None] * 2
        self.ticks_label = [None] * 2 
        self.label = [x['label'], y['label']]

    def set_tick(self, tick, tick_label, orient = 0):
        '''
        orient: 0: x-axis; 1: y-axis
        '''
        self.ticks[orient] = tick
        self.ticks_label[orient] = tick_label

        return self

class PlotManager(object):

    def __init__(self, title):        
        self.series_set = dict()
        self.title = title

    def add(self, series: DataSeries):
        name = series.name

        if name in self.series_set:
            raise ValueError(' existing data series with key=' + name)
        
        self.series_set[name] = series               

    def draw(self, img_file: str = "", dest_file: str = ""):

        fig, ax = None, None
        
        for key, series in self.series_set.items():

            if(fig == None):
                (fig, ax) = self.plot_xy_series(series=series, img_file=img_file)
            else:
                (fig, ax) = self.plot_xy_series(series=series, ax=ax)

        fig.suptitle(self.title)
        # save the figure 

        mkdir_filepath(dest_file)
        fig.savefig(dest_file)

        plt.close() # close the window once it is done. 

    def draw_ticks(self, ax: mp.pyplot.axes, series: DataSeries):
        '''
        xy = draw x ticks or y ticks, 0 x_axis; 1 y_axis
        '''
        i = 0
        for i in range(0, len(series.ticks)):
            tick = series.ticks[i]
            ticks_label = series.ticks_label[i]
            label = series.label[i]

            if i == 0:
                if series.ax_type == AxisType.Primary:
                    ax.set_xticks(tick)
                    ax.set_xticklabels(ticks_label)
                    ax.set_xlabel(label)
                else:
                    # ignore this top xtick for now
                    # secax = ax.secondary_xaxis('top')
                    # secax.set_xticks(tick)
                    # secax.set_xticklabels(label)
                    pass
            else:
                if series.ax_type == AxisType.Primary:
                    ax.set_yticks(tick)
                    ax.set_yticklabels(ticks_label)
                    ax.set_ylabel(label)
                else:
                    secax = ax.secondary_yaxis('right')
                    secax.set_yticks(tick)
                    secax.set_yticklabels(ticks_label)            
                    secax.set_ylabel(label)

    def plot_xy_series(self, series: DataSeries, ax: mp.pyplot.axes = [], img_file = []):
        ''' 
        plot a series of (x,y) on an figure
        ''' 
        x = series.x
        y = series.y
        range_x = series.range_x
        range_y = series.range_y
        color_style = series.cs
    
        if ax==[]:
            img = mpimg.imread(img_file)
            # flip the image upside down and set orign ='lower' when call imshow
            # to make the coordinates works in 'normal' sense
            img = np.flipud(img)
            f1 = plt.figure()
            ax = f1.add_subplot(111)
            ax.imshow(img, origin='lower')
            
        else:
            f1 = plt.gcf()

        # calculate x, y in axes coordinates
        x_axes = self.cal_coor_axes(x, range_x)
        y_axes = self.cal_coor_axes(y, range_y)
        
        # transfrom from axes coordinates to data coordinates
        x_data, y_data = self.coor_trans_axes_2_data(ax, x_axes, y_axes) 
        
        # draw data
        ax.plot(x_data, y_data, color_style)

        num_seg = 10

        tick, label = self.cal_tick_data(ax=ax, border= range_x, num_seg=10)
        series.set_tick(tick, label)

        tick, label = self.cal_tick_data(ax=ax, border=range_y, num_seg = 10, orient=1)
        series.set_tick(tick, label, orient=1)

        self.draw_ticks(ax, series)
            
        return (f1, ax)

    def cal_tick_data(self, ax, border, num_seg, orient = 0):
        
        (l, r) = (min(border), max(border))
        
        # tick = [0, 0.5, 1]
        # num_seg = len(tick)

        num_tick = num_seg + 1

        tick = (np.linspace(l, r, num_tick) - l) / (r - l)
        label = [''] * len(tick)
        label[0] = str(l)
        label[-1] = str(r)
        m = num_tick // 2
        label[m] = str(tick[m] * (r - l) + l)

        # label = list(map(str, [l, (l + r) / 2, r]))

        if orient == 0:
            (tick_data, t) = self.coor_trans_axes_2_data(ax, tick, np.zeros(num_tick))
        else:
            (t, tick_data) = self.coor_trans_axes_2_data(ax, np.zeros(num_tick), tick)        

        return tick_data, label

    def cal_coor_axes(self, x, x_r = []):
        '''
            calculate axes coordinates for x in range x_r or in its own range
        '''
        l, r = (min(x), max(x))
        if x_r != []:
            (l, r) = (min(x_r), max(x_r))

        return (np.array(x) - l) / (r - l)

    def coor_trans_axes_2_data(self, ax, x_axes, y_axes):
        '''
            transfrom points axes coordinates [(0, 0), (1, 1)]-> data coorinates [(0,0), (max_x, max_y)]
        '''
        # combine it into array of turple
        xy_axes = [(x_axes[i], y_axes[i]) for i in range(0, len(x_axes))]

        # transfrom axes coordinates to display coordinates
        xy_disp = ax.transAxes.transform(xy_axes)

        # transfrom from display coordinates to data coordinates
        xy_data = ax.transData.inverted().transform(xy_disp)

        x_data = [xy_data[i][0] for i in range(0, len(xy_data))]
        y_data = [xy_data[i][1] for i in range(0, len(xy_data))]
        
        return x_data, y_data