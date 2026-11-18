
import matplotlib.pyplot as plt
import matplotlib.pylab  as pylab
import matplotlib.lines as mlines

class Diagram ():

    def __init__ (self):
        self.figure = None
        self.axes   = None
        self.cs     = None
        self.x_list = []
        self.y_list = []

    def get_markers (self, size):
        markers = ['o','s','v','^','p','P','h','D']
        if size <= len(markers):
            markersout = markers [:size]
        else:
            markersout = markers [:]
            while len(markersout) < size:
                markersout.append ( markersout [-1])
        return markersout

    def get_colors (self, size, with_black=True):
        colors = ['black','red','blue','fuchsia','sienna','olivedrab','darkgoldenrod','darkgreen']
        if not with_black:
            colors [0] = 'forestgreen'
        if size <= len(colors):
            colorsout = colors [:size]
        else:
            colorsout = colors [:]
            while len(colorsout) < size:
                colorsout.append ( colorsout [-1])
        return colorsout


    def set_figure_size (self,size=(8,8)):
        self.figure = plt.figure ( figsize=size,facecolor='white')
        self.axes   = self.figure.add_subplot (1, 1, 1)
        plt.subplots_adjust(left=0.17, bottom=0.15, right=0.95, top=0.85, hspace = 0.05)


    def add_plot (self,x_,y_,label_="",linestyle_="-",color_="black", lw_=2,  \
                  marker_=None,markersize_=5,markevery_=1,markerfacecolor_='black', markeredgecolor_='black'):
        self.axes.plot (x_, y_, label=label_, linestyle=linestyle_, \
                        color = color_, lw=lw_, \
                        marker=marker_, markersize=markersize_,markevery=markevery_,\
                        markerfacecolor=markerfacecolor_, markeredgecolor=markeredgecolor_)
        self.x_list.append (x_)
        self.y_list.append (y_)

    def add_grid (self, lw_=0.4,ls_="-",which_="both"):
        self.axes.grid (True,which=which_,ls=ls_,lw=lw_)

    def add_legend (self,loc_="best",fontsize_=14):
        self.axes.legend (loc=loc_,fontsize=fontsize_)

    def set_x_label (self,label_,fontsize_=18):
        self.axes.set_xlabel (label_,fontsize=fontsize_)

    def set_y_label (self,label_,fontsize_=18):
        self.axes.set_ylabel (label_,fontsize=fontsize_)

    def save_fig (self,out_png_):
        self.figure.savefig (out_png_)

    def invert_x_axis (self):
        self.axes.invert_xaxis()

    def invert_y_axis (self):
        self.axes.invert_yaxis()

    def set_x_min (self,xmin_):
        self.axes.set_xlim (xmin=xmin_)

    def set_x_max (self,xmax_):
        self.axes.set_xlim (xmax=xmax_)

    def set_y_min (self,ymin_):
        self.axes.set_ylim (ymin=ymin_)

    def set_y_max (self,ymax_):
        self.axes.set_ylim (ymax=ymax_)

    def set_x_axes_tick_params (self,**kwargs):
        #help: https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.axes.Axes.tick_params.html
        self.axes.tick_params (**kwargs)
        #self.axes.tick_params (labelsize=25)

    def set_x_log_scale (self):
        plt.xscale ('log')

    def set_y_log_scale (self):
        plt.yscale ('log')

    def fill_between_curves (self, index1, index2, color):
        self.axes.fill_between ( self.x_list [index1], self.y_list [index1], self.y_list [index2],\
                                 where=self.y_list [index1] >= self.y_list [index2],\
                                 facecolor=color, interpolate=True)
        self.axes.fill_between ( self.x_list [index1], self.y_list [index2], self.y_list [index1],\
                                 where=self.y_list [index2] >= self.y_list [index1],\
                                 facecolor=color, interpolate=True)

    def add_text (self,x,y,text_,fontsize_=16, horizontalalignment_='center',verticalalignment_='center'):
        plt.text (x,y,text_,fontsize=fontsize_, horizontalalignment=horizontalalignment_,\
                  verticalalignment=verticalalignment_)

    def add_arrow (self,x_,y_,dx_,dy_,head_with_,head_length_,color_='red',double_side_=False):
        self.axes.arrow (x_,y_,dx_,dy_,head_width=head_with_,head_length=head_length_,color=color_)
        if double_side_:
            self.axes.arrow (x_+dx_, y_+dy_, -dx_, -dy_, head_width=head_with_, head_length=head_length_, color=color_)

    def set_title (self,title_):
        self.axes.set_title (title_)

    def add_contour(self, xi_, yi_, zi_, map_levels_=10, cmap_=pylab.cm.jet):
        self.cs = self.axes.contourf(xi_, yi_, zi_, map_levels_, cmap=cmap_)
        self.axes.contour (xi_, yi_, zi_, map_levels_, colors='k')

    def add_contour_bar(self):
        if self.cs != None:
            self.figure.colorbar(self.cs)

def main():
    dg = Diagram ()
    dg.set_figure_size((6,6))
    dg.set_title ("M_XX-(A-W) (-)")


#if __name__ == '__main__':
#    main()