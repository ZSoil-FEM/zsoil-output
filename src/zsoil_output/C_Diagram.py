"""Build and format two-dimensional result diagrams with Matplotlib.

This module is part of ``zsoil-output``, a Python SDK developed for use with ZSoil FEM to
read model data and output-result files.

The :class:`Diagram` wrapper manages curves, axes, legends, annotations, styles, and image
export for plots produced from ZSoil results.
"""


import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib.lines as mlines


class Diagram ():

    def __init__(self):
        """Initialize the diagram object for the supplied project data."""
        self.figure = None
        self.axes = None
        self.cs = None
        self.x_list = []
        self.y_list = []

    def get_markers(self, size):
        """Return a repeating sequence of plot markers of the requested length.

        :param size: Size used by this operation.
        :type size: int
        :return: A repeating sequence of plot markers of the requested length.
        :rtype: list[str]
        """
        markers = ['o', 's', 'v', '^', 'p', 'P', 'h', 'D']
        if size <= len(markers):
            markersout = markers[:size]
        else:
            markersout = markers[:]
            while len(markersout) < size:
                markersout.append(markersout[-1])
        return markersout

    def get_colors(self, size, with_black=True):
        """Return a repeating sequence of plot colors of the requested length.

        :param size: Size used by this operation.
        :type size: int
        :param with_black: With black used by this operation.
        :type with_black: bool
        :return: A repeating sequence of plot colors of the requested length.
        :rtype: list[str]
        """
        colors = ['black', 'red', 'blue', 'fuchsia', 'sienna', 'olivedrab', 'darkgoldenrod', 'darkgreen']
        if not with_black:
            colors[0] = 'forestgreen'
        if size <= len(colors):
            colorsout = colors[:size]
        else:
            colorsout = colors[:]
            while len(colorsout) < size:
                colorsout.append(colorsout[-1])
        return colorsout

    def set_figure_size(self, size=(8, 8)):
        """Set figure size for subsequent calculations.

        :param size: Size used by this operation.
        :type size: int
        """
        self.figure = plt.figure(figsize=size, facecolor='white')
        self.axes = self.figure.add_subplot(1, 1, 1)
        plt.subplots_adjust(left=0.17, bottom=0.15, right=0.95, top=0.85, hspace=0.05)

    def add_plot(self, x_, y_, label_="", linestyle_="-", color_="black", lw_=2,
                 marker_=None, markersize_=5, markevery_=1, markerfacecolor_='black', markeredgecolor_='black'):
        """Add plot to the current diagram.

        :param x_: X used by this operation.
        :type x_: Union[float, array_like]
        :param y_: Y used by this operation.
        :type y_: Union[float, array_like]
        :param label_: Label used to identify the item.
        :type label_: str
        :param linestyle_: Linestyle used by this operation.
        :type linestyle_: str
        :param color_: Color used by this operation.
        :type color_: str
        :param lw_: Lw used by this operation.
        :type lw_: float
        :param marker_: Marker used by this operation.
        :type marker_: Optional[str]
        :param markersize_: Markersize used by this operation.
        :type markersize_: float
        :param markevery_: Markevery used by this operation.
        :type markevery_: int
        :param markerfacecolor_: Markerfacecolor used by this operation.
        :type markerfacecolor_: str
        :param markeredgecolor_: Markeredgecolor used by this operation.
        :type markeredgecolor_: str
        """
        self.axes.plot(x_, y_, label=label_, linestyle=linestyle_,
                       color=color_, lw=lw_,
                       marker=marker_, markersize=markersize_, markevery=markevery_,
                       markerfacecolor=markerfacecolor_, markeredgecolor=markeredgecolor_)
        self.x_list.append(x_)
        self.y_list.append(y_)

    def add_grid(self, lw_=0.4, ls_="-", which_="both"):
        """Add grid to the current diagram.

        :param lw_: Lw used by this operation.
        :type lw_: float
        :param ls_: Ls used by this operation.
        :type ls_: str
        :param which_: Which used by this operation.
        :type which_: str
        """
        self.axes.grid(True, which=which_, ls=ls_, lw=lw_)

    def add_legend(self, loc_="best", fontsize_=14):
        """Add legend to the current diagram.

        :param loc_: Loc used by this operation.
        :type loc_: str
        :param fontsize_: Fontsize used by this operation.
        :type fontsize_: float
        """
        self.axes.legend(loc=loc_, fontsize=fontsize_)

    def set_x_label(self, label_, fontsize_=18):
        """Set x label for subsequent calculations.

        :param label_: Label used to identify the item.
        :type label_: str
        :param fontsize_: Fontsize used by this operation.
        :type fontsize_: float
        """
        self.axes.set_xlabel(label_, fontsize=fontsize_)

    def set_y_label(self, label_, fontsize_=18):
        """Set y label for subsequent calculations.

        :param label_: Label used to identify the item.
        :type label_: str
        :param fontsize_: Fontsize used by this operation.
        :type fontsize_: float
        """
        self.axes.set_ylabel(label_, fontsize=fontsize_)

    def save_fig(self, out_png_):
        """Save fig to the requested output file.

        :param out_png_: Path of the image file to create.
        :type out_png_: str
        """
        self.figure.savefig(out_png_)

    def invert_x_axis(self):
        """Invert x axis on the current diagram."""
        self.axes.invert_xaxis()

    def invert_y_axis(self):
        """Invert y axis on the current diagram."""
        self.axes.invert_yaxis()

    def set_x_min(self, xmin_):
        """Set x min for subsequent calculations.

        :param xmin_: Xmin used by this operation.
        :type xmin_: Union[float, array_like]
        """
        self.axes.set_xlim(xmin=xmin_)

    def set_x_max(self, xmax_):
        """Set x max for subsequent calculations.

        :param xmax_: Xmax used by this operation.
        :type xmax_: Union[float, array_like]
        """
        self.axes.set_xlim(xmax=xmax_)

    def set_y_min(self, ymin_):
        """Set y min for subsequent calculations.

        :param ymin_: Ymin used by this operation.
        :type ymin_: Union[float, array_like]
        """
        self.axes.set_ylim(ymin=ymin_)

    def set_y_max(self, ymax_):
        """Set y max for subsequent calculations.

        :param ymax_: Ymax used by this operation.
        :type ymax_: Union[float, array_like]
        """
        self.axes.set_ylim(ymax=ymax_)

    def set_x_axes_tick_params(self, **kwargs):
        """Set x axes tick params for subsequent calculations.

        .. seealso::
           `Matplotlib Axes.tick_params <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.tick_params.html>`_

        :param **kwargs: Additional keyword options passed to the underlying operation.
        :type **kwargs: dict
        """
        self.axes.tick_params(**kwargs)
        #self.axes.tick_params (labelsize=25)

    def set_x_log_scale(self):
        """Set x log scale for subsequent calculations."""
        plt.xscale('log')

    def set_y_log_scale(self):
        """Set y log scale for subsequent calculations."""
        plt.yscale('log')

    def fill_between_curves(self, index1, index2, color):
        """Fill between curves on the current diagram.

        :param index1: Index1 used by this operation.
        :type index1: int
        :param index2: Index2 used by this operation.
        :type index2: int
        :param color: Color used by this operation.
        :type color: str
        """
        self.axes.fill_between(self.x_list[index1], self.y_list[index1], self.y_list[index2],
                               where=self.y_list[index1] >= self.y_list[index2],
                               facecolor=color, interpolate=True)
        self.axes.fill_between(self.x_list[index1], self.y_list[index2], self.y_list[index1],
                               where=self.y_list[index2] >= self.y_list[index1],
                               facecolor=color, interpolate=True)

    def add_text(self, x, y, text_, fontsize_=16, horizontalalignment_='center', verticalalignment_='center'):
        """Add text to the current diagram.

        :param x: X used by this operation.
        :type x: Union[float, array_like]
        :param y: Y used by this operation.
        :type y: Union[float, array_like]
        :param text_: Text used by this operation.
        :type text_: str
        :param fontsize_: Fontsize used by this operation.
        :type fontsize_: float
        :param horizontalalignment_: Horizontalalignment used by this operation.
        :type horizontalalignment_: str
        :param verticalalignment_: Verticalalignment used by this operation.
        :type verticalalignment_: str
        """
        plt.text(x, y, text_, fontsize=fontsize_, horizontalalignment=horizontalalignment_,
                 verticalalignment=verticalalignment_)

    def add_arrow(self, x_, y_, dx_, dy_, head_with_, head_length_, color_='red', double_side_=False):
        """Add arrow to the current diagram.

        :param x_: X used by this operation.
        :type x_: Union[float, array_like]
        :param y_: Y used by this operation.
        :type y_: Union[float, array_like]
        :param dx_: Dx used by this operation.
        :type dx_: float
        :param dy_: Dy used by this operation.
        :type dy_: float
        :param head_with_: Head with used by this operation.
        :type head_with_: float
        :param head_length_: Head length used by this operation.
        :type head_length_: float
        :param color_: Color used by this operation.
        :type color_: str
        :param double_side_: Double side used by this operation.
        :type double_side_: bool
        """
        self.axes.arrow(x_, y_, dx_, dy_, head_width=head_with_, head_length=head_length_, color=color_)
        if double_side_:
            self.axes.arrow(x_+dx_, y_+dy_, -dx_, -dy_, head_width=head_with_, head_length=head_length_, color=color_)

    def set_title(self, title_):
        """Set title for subsequent calculations.

        :param title_: Title shown on the generated figure.
        :type title_: str
        """
        self.axes.set_title(title_)

    def add_contour(self, xi_, yi_, zi_, map_levels_=10, cmap_=pylab.cm.jet):
        """Add contour to the current diagram.

        :param xi_: Xi used by this operation.
        :type xi_: array_like
        :param yi_: Yi used by this operation.
        :type yi_: array_like
        :param zi_: Zi used by this operation.
        :type zi_: array_like
        :param map_levels_: Map levels used by this operation.
        :type map_levels_: int or sequence of float
        :param cmap_: Cmap used by this operation.
        :type cmap_: matplotlib.colors.Colormap
        """
        self.cs = self.axes.contourf(xi_, yi_, zi_, map_levels_, cmap=cmap_)
        self.axes.contour(xi_, yi_, zi_, map_levels_, colors='k')

    def add_contour_bar(self):
        """Add contour bar to the current diagram."""
        if self.cs != None:
            self.figure.colorbar(self.cs)


def main():
    """Run this module's example or command-line workflow."""
    dg = Diagram()
    dg.set_figure_size((6, 6))
    dg.set_title("M_XX-(A-W) (-)")


# if __name__ == '__main__':
#    main()
