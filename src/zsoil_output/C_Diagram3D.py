"""Build three-dimensional result diagrams with Matplotlib.

This module is part of ``zsoil-output``, a Python SDK developed for use with ZSoil FEM to
read model data and output-result files.

The :class:`Diagram3D` wrapper extends :class:`Diagram` with three-dimensional curves and
view settings for spatial result visualization.
"""


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import get_test_data
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
from .C_Diagram import *


class Diagram3D (Diagram):

    def __init__(self):
        """Initialize an empty 3D diagram."""
        pass

    def set_figure_size(self, size=(8, 8)):
        """Set figure size for subsequent calculations.

        :param size: Size used by this operation.
        :type size: int
        """
        self.figure = plt.figure(figsize=size, facecolor='white')
        self.axes = self.figure.add_subplot(1, 1, 1, projection='3d')
        plt.subplots_adjust(left=0.17, bottom=0.1, right=0.95, top=0.85, hspace=0.05)

    def add_plot(self, x_, y_, z_, label_='', linestyle_="-", color_="black", marker_=None, markersize_=5):
        """Add plot to the current diagram.

        :param x_: X used by this operation.
        :type x_: Union[float, array_like]
        :param y_: Y used by this operation.
        :type y_: Union[float, array_like]
        :param z_: Z used by this operation.
        :type z_: Union[float, array_like]
        :param label_: Label used to identify the item.
        :type label_: str
        :param linestyle_: Linestyle used by this operation.
        :type linestyle_: str
        :param color_: Color used by this operation.
        :type color_: str
        :param marker_: Marker used by this operation.
        :type marker_: Optional[str]
        :param markersize_: Markersize used by this operation.
        :type markersize_: float
        """
        self.axes.plot3D(x_, y_, z_, label=label_, linestyle=linestyle_,
                         color=color_,
                         marker=marker_, markersize=markersize_)

    def set_z_label(self, label_, fontsize_=18):
        """Set z label for subsequent calculations.

        :param label_: Label used to identify the item.
        :type label_: str
        :param fontsize_: Fontsize used by this operation.
        :type fontsize_: float
        """
        self.axes.set_zlabel(label_, fontsize=fontsize_)
