
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import get_test_data
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
from .C_Diagram import *

class Diagram3D (Diagram):

    def __init__ (self):
        pass


    def set_figure_size (self,size=(8,8)):
        self.figure = plt.figure ( figsize=size,facecolor='white')
        self.axes   = self.figure.add_subplot (1, 1, 1, projection ='3d')
        plt.subplots_adjust(left=0.17, bottom=0.1, right=0.95, top=0.85, hspace = 0.05)

    def add_plot (self,x_, y_, z_, label_='', linestyle_="-",color_="black", \
                  marker_=None,markersize_=5):
        self.axes.plot3D (x_, y_, z_, label=label_, linestyle=linestyle_, \
                         color = color_, \
                         marker=marker_, markersize=markersize_)

    def set_z_label (self,label_,fontsize_=18):
        self.axes.set_zlabel (label_,fontsize=fontsize_)