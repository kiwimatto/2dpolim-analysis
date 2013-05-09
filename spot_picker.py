import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np


class Annotate(object):
    def __init__(self,ax=None):
        if ax==None:
            self.ax = plt.gca()
        else:
            self.ax = ax
        self.rect = Rectangle((0,0), 0, 0, facecolor='yellow', edgecolor='yellow', alpha=.5, \
                                  zorder=10, clip_on=False)
        self.x0 = None
        self.y0 = None
        self.x1 = None
        self.y1 = None
        self.is_pressed = False
        self.ax.add_patch(self.rect)
        self.ax.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.ax.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.ax.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)
    def on_press(self, event):
        self.is_pressed = True
        print 'press'
        self.x0 = event.xdata
        self.y0 = event.ydata    
        self.x1 = event.xdata
        self.y1 = event.ydata
        if not np.any( [c==None for c in [self.x0,self.y0,self.x1,self.y1]] ):
            self.rect.set_width(self.x1 - self.x0)
            self.rect.set_height(self.y1 - self.y0)
            self.rect.set_xy((self.x0, self.y0))
            self.rect.set_linestyle('dashed')
        self.ax.figure.canvas.draw()
    def on_motion(self,event):
        if self.is_pressed is True:
            self.x1 = event.xdata
            self.y1 = event.ydata
            if not np.any( [c==None for c in [self.x0,self.y0,self.x1,self.y1]] ):
                self.rect.set_width(self.x1 - self.x0)
                self.rect.set_height(self.y1 - self.y0)
                self.rect.set_xy((self.x0, self.y0))
                self.rect.set_linestyle('dashed')
                self.ax.figure.canvas.draw()
    def on_release(self, event):
        self.is_pressed = False
        print 'release'
        self.x1 = event.xdata
        self.y1 = event.ydata
        if not np.any( [c==None for c in [self.x0,self.y0,self.x1,self.y1]] ):
            self.rect.set_width(self.x1 - self.x0)
            self.rect.set_height(self.y1 - self.y0)
            self.rect.set_xy((self.x0, self.y0))
            self.rect.set_linestyle('solid')
            self.ax.figure.canvas.draw()
        print self.x0,self.y0,self.x1,self.y1


def pick_a_spot( frame ):
    plt.imshow( frame )
    a = Annotate()
    plt.show()

    coords = np.round( np.array( [a.x0, a.y0, a.x1, a.y1] ) )
    print coords 

    return coords


if __name__=="__main__":
    pick_a_spot( np.random.random((512,512)) )
