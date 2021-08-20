#! /usr/bin/env python
"""
Code for making a simple position plot of Vizier SED values.
"""

import sys
import math
import tkinter as Tk
import numpy
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
BGCOL = '#F8F8FF'


def make_position_plot(data_values, nset):
    """

    Parameters
    ----------
    data_values : structure
        A structure which includes the data values needed for the plot.
    nset : integer
        The index value for the values to plot.

    Returns
    -------
    root1 : Tkinter root variable
        The window root variable.
    myplot : plot object
        The object for the position plot.

    """
    root1 = Tk.Tk()
    root1.title('Positions Window')
    myplot = PositionPlot(root1)
    myplot.data = data_values[nset]
    myplot.plot_positions()
    return root1, myplot


class PositionPlot(Tk.Frame):
    """
    Class to produce a Tkinter window with a position plot.

    Parameters
    ----------

    Tk.Frame    A Tkinter root or top level variable to hold the plot

    Returns
    -------

    No values are returned by this routine.
    """
    def __init__(self, parent, **args):
        if sys.version_info[0] == 2:
            raise ValueError('Python version 2 is required for the code.')
        if parent is None:
            raise ValueError('A root or top level window is required in ' +
                             'the call to sed_plot_window')
        # initialize the window and make the plot area.
        Tk.Frame.__init__(self, parent, args)
        self.root = parent
        self.data = None
        self.make_window()
        self.plot_positions()

    def make_window(self):
        """
        Code to make the plot in a new window.
        """
        frame1 = Tk.Frame(self.root)
        frame1.pack(side=Tk.TOP)
        frame1.config(bg=BGCOL)
        figure_area = Figure(figsize=(5., 5.), dpi=100)
        figure_subplot = figure_area.add_subplot(1, 1, 1)
        position_label = Tk.Label(frame1, text='Position:')
        position_label.pack(side=Tk.TOP)
        position_label.config(bg=BGCOL)
        main_canvas = FigureCanvasTkAgg(figure_area, master=frame1)
        # Here are defined the events that the program responds to for
        # the plot area.
        main_canvas.mpl_connect("motion_notify_event", self.__position_text)
        main_canvas.mpl_connect("key_press_event", self.__key_commands)
        main_canvas.mpl_connect("button_press_event", self.__plot_marker_set)
        main_canvas.mpl_connect("button_release_event",
                                self.__plot_marker_release)
        main_canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH,
                                         expand=Tk.YES)
        main_canvas.draw()
        self.plot_vars = {}
        self.plot_vars['canvas'] = main_canvas
        self.plot_vars['figure'] = figure_area
        self.plot_vars['subplot'] = figure_subplot
        self.plot_vars['position_label'] = position_label
        self.plot_positions()
        button_frame = Tk.Frame(frame1)
        button_frame.pack(side=Tk.TOP)
        close_button = Tk.Button(button_frame, text='Close',
                                 command=self.root.destroy)
        close_button.pack(side=Tk.TOP)
        close_button.config(bg=BGCOL)

    def plot_positions(self):
        """
        This routine does the actual plotting of the position values.
        """
        try:
            data_values = self.data['distance']
            xvalues = numpy.squeeze(data_values[2, :])
        except (ValueError, TypeError):
            return
        xvalues = numpy.squeeze(data_values[2, :])
        yvalues = numpy.squeeze(data_values[3, :])
        subplot = self.plot_vars['subplot']
        canvas = self.plot_vars['canvas']
        subplot.plot(xvalues, yvalues, marker='o', markersize=2.0,
                     linestyle='none')
        xmin, xmax = subplot.get_xbound()
        ymin, ymax = subplot.get_ybound()
        newmax = max(abs(xmax), abs(ymax), abs(ymin), abs(ymax))
        newmin = -newmax
        subplot.set_xbound(newmin, newmax)
        subplot.set_ybound(newmin, newmax)
        subplot.invert_xaxis()
        subplot.tick_params(axis='x', direction='in')
        subplot.tick_params(axis='y', direction='in')
        subplot.tick_params(bottom=True, top=True)
        subplot.tick_params(left=True, right=True)
        subplot.set_xlabel('RA Offset (arc-seconds)')
        subplot.set_ylabel('Dec Offset (arc-seconds)')
        canvas.draw()

    def __key_commands(self, event):
        pass

    def __position_text(self, event):
        """

        Parameters
        ----------

        event :   A motplotlib event variable.

        Returns
        -------

        No values are returned.

        """
        if (event.xdata is None) or (event.ydata is None):
            return
        dmin, indmin = self.__match_point(event.xdata, event.ydata)
        if dmin is None:
            return
        str1 = 'Position: %.3f %.3f\n' % (event.xdata, event.ydata)
        str1 = str1 + 'Nearest point: %.3f %.3f\n' % (
            self.data['distance'][2, indmin], self.data['distance'][3, indmin])
        str1 = str1 + 'Point %d, distance %.3f label %s\n' % (
            indmin+1, dmin, self.data['filter_name'][indmin])
        radius = math.sqrt(self.data['distance'][2, indmin]**2 +
                           self.data['distance'][3, indmin]**2)
        str1 = str1 + 'Radius from centre: %.3f arc-seconds' % (radius)
        self.plot_vars['position_label'].config(text=str1)

    def __match_point(self, xdata, ydata):
        """

        Parameters
        ----------
        xdata:   a float value, the xposition for the matching

        ydata:   a float value, the y position for the matching

        Returns
        -------
        radius:   a float value, the radius of the nearest data point from
                  the cursor position

        minind:   an integer value, the index of the closest data point in
                  the data array

        """
        try:
            xvalues = numpy.squeeze(self.data['distance'][2, :])
            yvalues = numpy.squeeze(self.data['distance'][3, :])
        except TypeError:
            return None, None
        radius = numpy.sqrt((xvalues - xdata)**2 + (yvalues - ydata)**2)
        minind = numpy.argmin(radius)
        return radius[minind], minind

    def __plot_marker_set(self, event):
        pass

    def __plot_marker_release(self, event):
        pass
