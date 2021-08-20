#! /usr/bin/env python
#
"""
This is the main driver routien for the SED display and analysis.

The normal use would be to invoke this from the command line as in

sed_plot_interface.py

There are no parameters for the call.

This code requires the Python extinction package.  Installation of the
package is described at "https://extinction.readthedocs.io/en/latest/".

Other required packages:  tkinter, matplotlib, numpy, scipy, sys, os, math,
and astropy.  All these are common packages.

"""
import math
import sys
import os
import tkinter as Tk
import tkinter.ttk
import tkinter.filedialog
import tkinter.simpledialog
import tkinter.messagebox
from tkinter.colorchooser import askcolor
from tkinter.scrolledtext import ScrolledText
from tkinter.filedialog import askopenfilenames
from tkinter.simpledialog import askinteger
import numpy
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib
from matplotlib.figure import Figure
import read_vizier_sed_table
import sed_utilities
import tkinter_utilities
import extinction_code
import model_utilities
import position_plot
matplotlib.use('TkAgg')

# The following are "global" variables with line/marker information from
# matplotlib.  These are used, but not changed, in the code in more than one
# place hence I am using single variables here for the values.
MATPLOTLIB_SYMBOL_LIST = [
    None, "o", "v", "^", "<", ">", "1", "2", "3", "4", "8", "s", "p",
    "P", "*", "h", "H", "+", "x", "X", "D", "d", "|", "_", ".", "."]
MATPLOTLIB_SYMBOL_NAME_LIST = [
    'None', 'circle', 'triangle down', 'triangle up', 'triangle_left',
    'triangle right', 'tri_down', 'tri_up', 'tri_left', 'tri_right',
    'octagon', 'square', 'pentagon', 'plus (filled)', 'star',
    'hexagon1', 'hexagon2', 'plus', 'x', 'x (filled)', 'diamond',
    'thin_diamond', 'vline', 'hline', 'point', 'pixel']
MATPLOTLIB_LINE_LIST = ['-', '--', '-.', ':', None]
MATPLOTLIB_LINE_NAME_LIST = ['solid', 'dashed', 'dashdot', 'dotted', 'None']
# define a background colour for windows
BGCOL = '#F8F8FF'
# default colours
COLOUR_SET = ['blue', 'forestgreen', 'black', 'orange', 'red', 'purple',
              'cyan', 'lime', 'brown', 'violet', 'grey', 'gold']
WAVELNORM = 2.2
FLAMBDANORM = 5.e-15


def startup():
    """
    Startup.py is a wrapper for starting the plot tool.

    Parameters
    ----------
        None.

    Returns
    -------
        root :   The Tkinter window class variable for the plot window.

        plot_object :  The SED plot GUI object variable.

    """
    # Make a Tkinter window
    newroot = Tk.Tk()
    newroot.title("SED Fitting Tool")
    newroot.config(bg=BGCOL)
    # start the interface
    plot_object = PlotWindow(newroot)
    return newroot, plot_object


def strfmt(value):
    """
    Apply a format to a real value for writing out the data values.

    This routine is used to format an input real value either in
    exponential format or in floating point format depending on
    the magnitude of the input value.
    This works better for constant width columns than the Python g format.

    Parameters
    ----------
        value :   a real number value

    Returns
    -------
        outstring :  a format string segment

    """
    if (abs(value) > 1.e+07) or (abs(value) < 0.001):
        outstring = '%14.7e ' % (value)
        if value == 0.:
            outstring = '%14.6f ' % (value)
    else:
        outstring = '%14.6f ' % (value)
    outstring = outstring.lstrip(' ')
    return outstring


def get_data_values(xdata, ydata, mask, filternames):
    """
    Utility routine to apply a mask to (x,y) data.

    Parameters
    ----------
        xdata:  A numpy array of values, nominally float values

        ydata:  A numpy array of values, nominally float values

        mask:   A numpy boolean array for which points are to be returned

        filternames:  A numpy string array of the filter names

    Returns
    -------
        newxdata:   A numpy array of the xdata values where mask = True

        newydata:   A numpy array of the ydata values where mask = True

        newfilternames:  A numpy array of the filter names where mask = True

    """
    inds = numpy.where(mask)
    newxdata = numpy.copy(xdata[inds])
    newydata = numpy.copy(ydata[inds])
    try:
        newfilternames = numpy.copy(filternames[inds])
    except TypeError:
        newfilternames = []
        for loop in range(len(newxdata)):
            newfilternames.append('')
        newfilternames = numpy.asarray(newfilternames)
    inds = numpy.argsort(newxdata)
    newxdata = newxdata[inds]
    newydata = newydata[inds]
    newfilternames = newfilternames[inds]
    return newxdata, newydata, newfilternames


def unpack_vot(votdata):
    """
    Utility routine to unpack Vizier VOT photometry data.

    Parameters
    ----------
        votdata:  A table of Vizier photometry values from the read_vizier_vot
                  function

    Returns
    -------

        data_set:  A list containing a variety of values read from the VOT file

    """
    photometry_values = numpy.copy(votdata[0])
    error_mask = numpy.copy(votdata[1])
    filter_names = numpy.copy(votdata[2])
    references = numpy.copy(votdata[3])
    refpos = numpy.copy(votdata[4])
    data_set = {'wavelength': None, 'frequency': None, 'fnu': None,
                'flambda': None, 'lfl': None, 'l4fl': None,
                'dfnu': None, 'dflambda': None, 'dl4fl': None,
                'dlfl': None, 'mask': None, 'plot_mask': None,
                'filter_name': None, 'distance': None,
                'position': None, 'refpos': None,
                'references': None, 'plot': None, 'source': None,
                'colour_by_name': False}
    data_set['wavelength'] = numpy.squeeze(photometry_values[0, :])
    data_set['frequency'] = numpy.squeeze(photometry_values[1, :])
    data_set['fnu'] = numpy.squeeze(photometry_values[4, :])
    data_set['flambda'] = numpy.squeeze(photometry_values[6, :])
    data_set['lfl'] = numpy.squeeze(photometry_values[8, :])
    data_set['l4fl'] = numpy.squeeze(photometry_values[8, :]) *\
        (data_set['wavelength']**3)
    data_set['dfnu'] = numpy.squeeze(photometry_values[5, :])
    data_set['dflambda'] = numpy.squeeze(photometry_values[7, :])
    data_set['dlfl'] = numpy.squeeze(photometry_values[9, :])
    data_set['dl4fl'] = numpy.squeeze(photometry_values[9, :]) *\
        (data_set['wavelength']**3)
    data_set['filter_name'] = numpy.copy(filter_names)
    data_set['position'] = [numpy.squeeze(photometry_values[2, :]),
                            numpy.squeeze(photometry_values[3, :])]
    data_set['refpos'] = numpy.copy(refpos)
    data_set['references'] = numpy.copy(references)
    data_set['mask'] = numpy.copy(error_mask)
    data_set['plot_mask'] = numpy.copy(error_mask)
    for loop in range(len(error_mask)):
        data_set['plot_mask'][loop] = True
    data_set['distance'] = numpy.copy(photometry_values[10:14, :])
    return data_set


class PlotWindow(Tk.Frame):
    """
    This is the class for the plotting window.

    Parameters
    ----------

    Tk.Frame:     A Tkinter Frame, root or Toplevel variable to hold
                  the plot window (more usually one of the latter two)

    Returns
    -------

    The class variable is returned.
    """

    def __init__(self, parent, **args):
        """
        This routine sets a few variables for the interface and calls the
        main widget function.

        Parameters
        ----------
            parent:     A parameter giving the parent Tk root window name.

            **args      A possible list of additional arguments.  This is
                        currently not used.

        Returns
        -------
            No value is returned by this routine.

        """
        if sys.version_info[0] == 2:
            raise ValueError('Python version 2 is required for the code.')
        if parent is None:
            raise ValueError('A root or top level window is required in ' +
                             'the call to sed_plot_window')
        # initialize the window and make the plot area.
        Tk.Frame.__init__(self, parent, args)
        self.root = parent
        self.data_values = [{'wavelength': None, 'fnu': None,
                             'frequency': None,
                             'flambda': None, 'lfl': None, 'l4fl': None,
                             'dfnu': None, 'dflambda': None,
                             'dl4fl': None, 'dlfl': None, 'mask': None,
                             'filter_name': None, 'distance': None,
                             'position': None, 'refpos': None,
                             'references': None, 'plot': None,
                             'source': None, 'plot_mask': None,
                             'color_by_name': False}, ]
        self.stellar_models = [{'wavelength': None, 'fnu': None,
                                'flambda': None, 'lfl': None,
                                'l4fl': None, 'label': None}, ]
        self.tkinter_values = {'plot_frame': None, 'variables': [None, ],
                               'entry': None, 'button': None,
                               'radio': None, 'subplot': None,
                               'canvas': None, 'figure': None,
                               'set_parameters': None,
                               'modeltype': None}
        self.extinction_values = {'av': None, 'ebmv': None,
                                  'function': None, 'redden': True,
                                  'controls': None,
                                  'rl_wavelengths': None,
                                  'rl_extinction': None}
        self.flags = {'auto_scale': True, 'toggle_status': False,
                      'zoom_status': False, 'positions': [],
                      'model_scale': False, 'point_scale': False,
                      'template_spectrum': -1, 'use_ratio': False,
                      'mask_area': False}
        self.make_widget()

    def make_widget(self):
        """
        This routine makes the main plot window.

        The plot area and the control area are created here with Tkinter calls.
        Various object variables are defined in the process.

        Parameters
        ----------
            None

        Returns
        -------
            Nothing is returned from this routine.

        """
        menu_frame = Tk.Frame(self.root)
        menu_frame.pack(side=Tk.TOP, anchor=Tk.W)
        menu_frame.config(bg=BGCOL)
        self.__make_menus(menu_frame)
        control_frame = Tk.Frame(self.root)
        control_frame.pack(side=Tk.LEFT, fill=Tk.Y, expand=1)
        control_frame.config(bg=BGCOL)
        self.__make_controls(control_frame)
        tkinter_utilities.separator_line(self.root, [5, 750, 5],
                                         False, Tk.LEFT, BGCOL)
        plot_frame = Tk.Frame(self.root)
        plot_frame.pack(side=Tk.LEFT, fill=Tk.Y, expand=1)
        plot_frame.config(bg=BGCOL)
        self.__make_plot_area(plot_frame)
        self.tkinter_values['plot_frame'] = plot_frame

    def __make_menus(self, parent):
        """
        Create pull-down menus for tool functionality.

        Given a Tk Frame variable "parent" this routine makes a pull-down
        menu area within this frame.

        Parameters
        ----------
            parent     A Tk.Frame variable, that holds the menus

        Returns
        -------
            No values are returned by the routine.

        """
        menubutton1 = Tk.Menubutton(parent, text="Read/Save Data")
        menubutton1.pack(side=Tk.LEFT, fill=Tk.X, expand=2)
        menubutton1.config(bg=BGCOL)
        menu1 = Tk.Menu(menubutton1)
        menubutton1['menu'] = menu1
        menu1.add_command(label='Read Vizier SED file',
                          command=self.read_vizier_file)
        menu1.add_command(label='Read spectrum file',
                          command=self.read_spectrum_file)
        menu1.add_command(label='Save data set',
                          command=self.save_data_set)
        menu1.add_command(label='Read data set',
                          command=self.read_data_set)
        menubutton2 = Tk.Menubutton(parent, text="Data Operations")
        menubutton2.pack(side=Tk.LEFT, fill=Tk.X, expand=2)
        menubutton2.config(bg=BGCOL)
        menu2 = Tk.Menu(menubutton2)
        menubutton2['menu'] = menu2
        menu2.add_command(label='Delete/Undelete Data Points',
                          command=self.__toggle_data_status)
        menu2.add_command(label='Toggle Points in Region',
                          command=self.__toggle_mask_area)
        menu2.add_command(label='Unmask All Points',
                          command=self.unmask_all)
        menu2.add_command(label='Cursor Zoom',
                          command=self.__toggle_zoom_status)
        menu2.add_command(label='Auto-scale Plot',
                          command=self.__autoscale_plot)
        menu2.add_command(label='Freeze Range',
                          command=self.freeze_range)
        menu2.add_command(label='Average Vizier Data',
                          command=self.average_vizier_data)
        menu2.add_command(label='Toggle Vizier Filter Colours',
                          command=self.__toggle_colour_display)
        menu2.add_command(label='De-redden Data',
                          command=lambda: self.make_extinction_window(True))
        menu2.add_command(label='Calculate Flux',
                          command=self.calculate_flux)
        menu2.add_command(label='Show Position Offsets',
                          command=self.__plot_position_offsets)
        menubutton3 = Tk.Menubutton(parent, text="Stellar Models")
        menubutton3.pack(side=Tk.LEFT, fill=Tk.X, expand=2)
        menubutton3.config(bg=BGCOL)
        menu3 = Tk.Menu(menubutton3)
        menubutton3['menu'] = menu3
        menu3.add_command(label='Read Stelar Model',
                          command=self.stellar_model_window)
        menu3.add_command(label='Scale Model to Cursor',
                          command=self.__scale_model_to_cursor)
        menu3.add_command(label='Scale Model to Data Point',
                          command=self.__scale_model_to_point)
        menu3.add_command(label='Define Blackbody Spectrum',
                          command=self.__make_blackbody)
        menu3.add_command(label='Toggle Standard/Ratio Fit',
                          command=self.__toggle_ratio_flag)
        menu3.add_command(label='Fit Model to Data Points',
                          command=self.__fit_model_to_data)
        menu3.add_command(label='Redden Model',
                          command=lambda: self.make_extinction_window(False))
        menu3.add_command(label='Set Standard Spectrum',
                          command=self.set_standard_spectrum)
        menu3.add_command(label='Find Best Fit Model',
                          command=self.__get_model_names)

    def __make_blackbody(self):
        """
        Query the user for a temperature value, then make a blackbody spectrum
        using the standard spectral template wavelengths.  Add the spectrum to 
        the plot with the usual normalization.
        """
        temperature = tkinter.simpledialog.askfloat(
            'Query', 'Enter blackbody temperature:')
        if temperature <= 4.:
            str1 = 'The temperature %.3f is too low' % (temperature)
            tkinter.messagebox.showinfo('Error', str1)
            return
        try:
            path = os.environ['EXTINCTION_PATH']
            if path[-1] != '/':
                path = path + '/'
        except:
            path = './'
        try:
            wavelengths, spectrum = \
                model_utilities.read_calspec_model(path+'sirius_mod_004.fits')
        except:
            tkinter.messagebox.showinfo(
                'Error', 'Could not read in file sirius_mod_004.fits.')
        # The following are the standard radiation constants (CODATA 2018)
        # with wavelength units changed to microns
        c1 = 3.741771852e+08
        c2 = 14387.76877
        #
        factor = numpy.expm1(c2/(WAVELNORM*temperature))
        f0 = (c1/numpy.power(WAVELNORM, 5))/factor
        scale = FLAMBDANORM/f0
        exp1 = c2/(wavelengths*temperature)
        spectrum = scale*(c1/numpy.power(wavelengths, 5))/numpy.expm1(exp1)
        overflowinds = numpy.isnan(spectrum)
        spectrum[overflowinds] = 0.
        overflowinds=numpy.isinf(spectrum)
        spectrum[overflowinds] = 0.
        label = 'stellar model (blackbody %.3f K)' % (temperature)
        self.add_spectrum(wavelengths, spectrum, 0, 0, label)
        self.make_plot()

    def __toggle_colour_display(self):
        response = tkinter.simpledialog.askstring(
            'Input', 'Set to toggle colour status (or all):')
        if 'all' == response.lower():
            for loop in range(len(self.data_values)):
                if 'vizier_sed' in self.data_values[loop]['source']:
                    self.data_values[loop]['colour_by_name'] = not \
                        self.data_values[loop]['colour_by_name']
        else:
            try:
                loop = int(response) - 1
                if loop < len(self.data_values):
                    self.data_values[loop]['colour_by_name'] = not \
                        self.data_values[loop]['colour_by_name']
            except:
                pass
        self.make_plot()

    def __toggle_ratio_flag(self):
        self.flags['use_ratio'] = not self.flags['use_ratio']
        if self.flags['use_ratio']:
            str1 = 'The point by point signal ratio will be used' + \
                ' in model fitting.'
        else:
            str1 = 'The normal least squares calculation will be ' + \
                'used in model fitting.'
        tkinter.messagebox.showinfo('Information', str1)

    def __get_model_names(self):
        names = askopenfilenames()
        print(names)

    def __make_controls(self, parent):
        """
        Make the control area within the main window.

        This routine makes a control area within the main window, under
        frame "parent".  The overall root value is also passed here for
        closing the window.

        Parameters
        ----------
            parent :   A Tk.Frame variable for the holder of the controla

        Returns
        -------
            No values are returned by this routine.

        """
        holder = Tk.Frame(parent)
        holder.pack(side=Tk.TOP)
        holder.config(bg=BGCOL)
        label1 = Tk.Label(holder, text=' ')
        label1.pack(side=Tk.TOP, fill=Tk.X)
        label1.config(bg=BGCOL)
        button1 = Tk.Button(holder, text='(Re-)Plot', command=self.make_plot)
        button1.pack(side=Tk.TOP, fill=Tk.X)
        button1.config(bg=BGCOL)
        button2 = Tk.Button(holder, text='Auto-scale Plot',
                            command=self.__autoscale_plot)
        button2.pack(side=Tk.TOP, fill=Tk.X)
        button2.config(bg=BGCOL)
        button3 = Tk.Button(holder, text='Apply/Freeze Range',
                            command=self.freeze_range)
        button3.pack(side=Tk.TOP, fill=Tk.X)
        button3.config(bg=BGCOL)
        button4 = Tk.Button(
            holder, text='Set Properties', command=self.change_set_properties)
        button4.pack(side=Tk.TOP, fill=Tk.X)
        button4.config(bg=BGCOL)
        range_field = Tk.Frame(holder)
        range_field.pack(side=Tk.TOP)
        range_field.config(bg=BGCOL)
        label1 = Tk.Label(range_field, text='x min: ')
        label1.grid(column=0, row=0)
        label1.config(bg=BGCOL)
        xmin_field = Tk.Entry(range_field, width=15)
        xmin_field.grid(column=1, row=0)
        xmin_field.insert(0, ' ')
        label1 = Tk.Label(range_field, text='x max: ')
        label1.grid(column=0, row=1)
        label1.config(bg=BGCOL)
        xmax_field = Tk.Entry(range_field, width=15)
        xmax_field.grid(column=1, row=1)
        xmax_field.insert(0, ' ')
        label1 = Tk.Label(range_field, text='y min: ')
        label1.grid(column=0, row=2)
        label1.config(bg=BGCOL)
        ymin_field = Tk.Entry(range_field, width=15)
        ymin_field.grid(column=1, row=2)
        ymin_field.insert(0, ' ')
        label1 = Tk.Label(range_field, text='y max: ')
        label1.grid(column=0, row=3)
        label1.config(bg=BGCOL)
        ymax_field = Tk.Entry(range_field, width=15)
        ymax_field.grid(column=1, row=3)
        ymax_field.insert(0, ' ')
        option_area = Tk.Frame(holder)
        option_area.pack(side=Tk.TOP)
        option_area.config(bg=BGCOL)
        label1 = Tk.Label(option_area, text='y values type: ')
        label1.pack(side=Tk.TOP)
        ytypeflag = Tk.IntVar()
        label1.config(bg=BGCOL)
        opt1 = Tk.Radiobutton(option_area, text='lambda*F_lambda',
                              variable=ytypeflag, value=0,
                              command=self.make_plot)
        opt1.pack(side=Tk.TOP)
        opt1.config(bg=BGCOL)
        opt2 = Tk.Radiobutton(option_area, text='F_lambda',
                              variable=ytypeflag, value=1,
                              command=self.make_plot)
        opt2.pack(side=Tk.TOP)
        opt2.config(bg=BGCOL)
        opt3 = Tk.Radiobutton(option_area, text='F_nu',
                              variable=ytypeflag, value=2,
                              command=self.make_plot)
        opt3.pack(side=Tk.TOP)
        opt3.config(bg=BGCOL)
        opt4 = Tk.Radiobutton(option_area, text='lambda^4 * F_lambda',
                              variable=ytypeflag, value=3,
                              command=self.make_plot)
        opt4.pack(side=Tk.TOP)
        opt4.config(bg=BGCOL)
        ytypeflag.set(0)
        option_area = Tk.Frame(holder)
        option_area.pack(side=Tk.TOP)
        option_area.config(bg=BGCOL)
        label1 = Tk.Label(option_area, text='x values type: ')
        label1.pack(side=Tk.TOP)
        xtypeflag = Tk.IntVar()
        label1.config(bg=BGCOL)
        opt1 = Tk.Radiobutton(option_area, text='Wavelength',
                              variable=xtypeflag, value=0,
                              command=self.make_plot)
        opt1.pack(side=Tk.TOP)
        opt1.config(bg=BGCOL)
        opt2 = Tk.Radiobutton(option_area, text='Frequency',
                              variable=xtypeflag, value=1,
                              command=self.make_plot)
        opt2.pack(side=Tk.TOP)
        opt2.config(bg=BGCOL)
        xtypeflag.set(0)
        tkinter_utilities.separator_line(holder, [200, 25, 5], True,
                                         Tk.TOP, BGCOL)
        label1 = Tk.Label(holder, text='Apply Filter Extinction')
        label1.pack(side=Tk.TOP)
        label1.config(bg=BGCOL)
        h1 = Tk.Frame(holder)
        h1.pack(side=Tk.TOP)
        filter_extinction_flag = Tk.IntVar()
        opt1 = Tk.Radiobutton(h1, text='Yes',
                              variable=filter_extinction_flag, value=0)
        opt1.pack(side=Tk.LEFT)
        opt1.config(bg=BGCOL)
        opt2 = Tk.Radiobutton(h1, text='No',
                              variable=filter_extinction_flag, value=1)
        opt2.pack(side=Tk.TOP)
        opt2.config(bg=BGCOL)
        filter_extinction_flag.set(1)
        label1 = Tk.Label(holder, text='Use Filter Mean Flux Density')
        label1.pack(side=Tk.TOP)
        label1.config(bg=BGCOL)
        h1 = Tk.Frame(holder)
        h1.pack(side=Tk.TOP)
        filter_mean_flambda_flag = Tk.IntVar()
        opt1 = Tk.Radiobutton(h1, text='Yes',
                              variable=filter_mean_flambda_flag, value=0)
        opt1.pack(side=Tk.LEFT)
        opt1.config(bg=BGCOL)
        opt2 = Tk.Radiobutton(h1, text='No',
                              variable=filter_mean_flambda_flag, value=1)
        opt2.pack(side=Tk.TOP)
        opt2.config(bg=BGCOL)
        filter_mean_flambda_flag.set(1)
        tkinter_utilities.separator_line(holder, [200, 25, 5], True,
                                         Tk.TOP, BGCOL)
        button5 = Tk.Button(holder, text='Save as PNG',
                            command=self.save_as_png)
        button5.pack(side=Tk.TOP, fill=Tk.X)
        button5.config(bg=BGCOL)
        button6 = Tk.Button(holder, text='Save as Postscript',
                            command=self.save_as_postscript)
        button6.pack(side=Tk.TOP, fill=Tk.X)
        button6.config(bg=BGCOL)
        tkinter_utilities.separator_line(holder, [200, 25, 5], True,
                                         Tk.TOP, BGCOL)
        button7 = Tk.Button(holder, text='Clear Plot',
                            command=self.__clear_plot)
        button7.pack(side=Tk.TOP, fill=Tk.X)
        button7.config(bg=BGCOL)
        button8 = Tk.Button(holder, text='Close Window',
                            command=self.__close_main_window)
        button8.pack(side=Tk.TOP, fill=Tk.X)
        button8.config(bg=BGCOL)
        self.tkinter_values['entry'] = {
            'xmin': xmin_field, 'xmax': xmax_field,
            'ymin': ymin_field, 'ymax': ymax_field}
        if len(self.tkinter_values['variables']) == 1:
            self.tkinter_values['variables'].append(ytypeflag)
            self.tkinter_values['variables'].append(xtypeflag)
            self.tkinter_values['variables'].append(filter_extinction_flag)
            self.tkinter_values['variables'].append(filter_mean_flambda_flag)
        else:
            self.tkinter_values['variables'][1] = ytypeflag
            self.tkinter_values['variables'][2] = xtypeflag
            self.tkinter_values['variables'][3] = filter_extinction_flag
            self.tkinter_values['variables'][4] = filter_mean_flambda_flag

    def __close_main_window(self):
        response = tkinter.messagebox.askyesno(
            "Verify",
            "Do you want to quit the interface?")
        if not response:
            return
        self.root.destroy()

    def __clear_plot(self):
        response = tkinter.messagebox.askyesno(
            "Verify",
            "Do you want to abandon the current plot?")
        if not response:
            return
        self.data_values = [{'wavelength': None, 'fnu': None,
                             'frequency': None,
                             'flambda': None, 'lfl': None, 'l4fl': None,
                             'dfnu': None, 'dflambda': None, 'dl4fl': None,
                             'dlfl': None, 'mask': None,
                             'filter_name': None, 'distance': None,
                             'position': None, 'refpos': None,
                             'references': None, 'plot': None,
                             'source': None, 'plot_mask': None,
                             'colour_by_name': False}, ]
        self.stellar_models = [{'wavelength': None, 'fnu': None,
                                'flambda': None, 'lfl': None, 'l4fl': None}, ]
        self.extinction_values = {'av': None, 'ebmv': None,
                                  'function': None, 'redden': True,
                                  'controls': None,
                                  'rl_wavelengths': None,
                                  'rl_extinction': None}
        self.flags = {'auto_scale': True, 'toggle_status': False,
                      'zoom_status': False, 'positions': [],
                      'model_scale': False, 'point_scale': False,
                      'template_spectrum': -1, 'use_ratio': False,
                      'mask_area': False}
        self.make_plot()
        xmin, xmax, ymin, ymax = self.get_bounds(
            self.tkinter_values['subplot'])

    def average_vizier_data(self):
        for loop in range(len(self.data_values)):
            if self.data_values[loop]['source'] == 'vizier_sed':
                data_set = {
                    'wavelength': None, 'frequency': None,
                    'fnu': None, 'flambda': None, 'lfl': None, 'l4fl': None,
                    'dfnu': None, 'dflambda': None, 'dl4fl': None,
                    'dlfl': None, 'mask': None, 'plot_mask': None,
                    'filter_name': None, 'distance': None,
                    'position': None, 'refpos': None,
                    'references': None, 'plot': None, 'source': None,
                    'colour_by_name': False}
                xvalues = numpy.copy(self.data_values[loop]['wavelength'])
                yvalues = numpy.copy(self.data_values[loop]['flambda'])
                filter_names = numpy.copy(
                    self.data_values[loop]['filter_name'])
                inds = numpy.where(self.data_values[loop]['plot_mask'])
                xvalues = xvalues[inds]
                yvalues = yvalues[inds]
                filter_names = filter_names[inds]
                newxvalues, newyvalues, new_filter_names = \
                    sed_utilities.vizier_means(
                        xvalues, yvalues, filter_names)
                data_set['wavelength'] = newxvalues
                data_set['frequency'] = sed_utilities.trans_wavelength(
                    newxvalues, 0, 2)
                data_set['flambda'] = newyvalues
                data_set['fnu'] = sed_utilities.trans_flux_density(
                    newxvalues, newyvalues, 0, 5)
                data_set['lfl'] = sed_utilities.trans_flux_density(
                    newxvalues, newyvalues, 0, 2)
                data_set['l4fl'] = data_set['lfl'] * (
                    data_set['wavelength']**3)
                data_set['dflambda'] = newyvalues * 0.
                data_set['dfni'] = newyvalues * 0.
                data_set['dlfl'] = newyvalues * 0.
                data_set['d4lfl'] = newyvalues * 0.
                data_set['source'] = 'averaged_vizier_set'
                data_set['filter_name'] = new_filter_names
                mask = []
                for n1 in range(len(newxvalues)):
                    mask.append(True)
                data_set['mask'] = mask
                data_set['plot_mask'] = mask
                data_set['plot'] = {
                    'symbol': self.data_values[loop]['plot']['symbol'],
                    'size': self.data_values[loop]['plot']['size'],
                    'line': True,
                    'line_type': MATPLOTLIB_LINE_LIST[1],
                    'line_width': 1.0,
                    'colour': self.data_values[loop]['plot']['colour'],
                    'show': True}
                for newloop in range(len(COLOUR_SET)):
                    if data_set['plot']['colour'] == COLOUR_SET[newloop]:
                        n1 = newloop + 1
                if n1 == len(COLOUR_SET):
                    n1 = 0
                data_set['plot']['colour'] = COLOUR_SET[n1]
                self.data_values.append(data_set)
                self.make_plot()

    def calculate_flux(self):
        """
        Bring up a window for flux calculations.

        Parameters
        ----------

        None

        Returns
        -------

        Nothing
        """
        try:
            ndatasets = len(self.data_values)
            if (ndatasets == 0) or (self.data_values[0]['wavelength'] is None):
                return
        except ValueError:
            return
        window1 = Tk.Toplevel(self.root)
        window1.config(bg=BGCOL)
        window1.title('Flux Calculation Window')
        text_area = ScrolledText(window1, height=20, width=50,
                                 wrap=Tk.NONE)
        text_area.pack(side=Tk.TOP)
        text_area.config(font=('courier', 16, 'bold'))
        frame1 = Tk.Frame(window1)
        frame1.pack(side=Tk.TOP)
        frame1.config(bg=BGCOL)
        label1 = Tk.Label(frame1, text='Data set number:')
        label1.grid(column=0, row=0)
        set_number_menu = tkinter.ttk.Combobox(frame1, width=10)
        menu = []
        for loop in range(ndatasets):
            menu.append(str(loop+1))
        set_number_menu['values'] = menu
        set_number_menu.grid(column=1, row=0)
        set_number_menu.current(0)
        label1 = Tk.Label(frame1, text='Distance/Parallax:')
        label1.grid(column=0, row=1)
        parallax_entry = Tk.Entry(frame1, width=10)
        parallax_entry.grid(column=1, row=1)
        parallax_entry.insert(0, '1.0')
        label1 = Tk.Label(frame1, text='Distance/Parallax Uncertainty:')
        label1.grid(column=0, row=2)
        parallax_uncertainty_entry = Tk.Entry(frame1, width=10)
        parallax_uncertainty_entry.grid(column=1, row=2)
        parallax_uncertainty_entry.insert(0, '0.0')
        label1 = Tk.Label(frame1, text='Value:')
        label1.grid(column=0, row=3)
        parallax_variable = Tk.IntVar()
        button_frame = Tk.Frame(frame1)
        button_frame.grid(column=1, row=3)
        tkinter_utilities.put_yes_no(button_frame, parallax_variable,
                                     ['Distance (kpc)', 'Parallax (mas)'],
                                     False)
        apply_button = Tk.Button(
            frame1, text='Recalculate',
            command=lambda: self.calculate_fluxes(
                set_number_menu, parallax_variable,
                parallax_entry, parallax_uncertainty_entry,
                text_area))
        apply_button.grid(column=0, row=4)
        close_button = Tk.Button(
            frame1, text='Close', command=window1.destroy)
        close_button.grid(column=1, row=4)
        self.calculate_fluxes(set_number_menu, parallax_variable,
                              parallax_entry, parallax_uncertainty_entry,
                              text_area)

    def calculate_fluxes(self, set_number_menu, parallax_variable,
                         parallax_entry, parallax_uncertainty_entry,
                         text_area):
        """
        This routine calculates the flux from integration over the
        spectral energy distribution.

        Parameters
        ----------

        set_number menu:    A Tkinter menu variable

        parallax_variable:  A Tkinter int variable, gives the parallax
                            option (parallax in mas or distance in kpc)

        parallax_entry:     A Tkinter entry variable for the parallax or
                            distance

        parallax_uncertainty_enrtry:   A Tkinter entry variable for the
                                       parallax or distance uncertainty

        text_area:   A TKinter text box or scrolled text variable where
                     messages are posted
        """
        nset = set_number_menu.current()
        wl0 = numpy.copy(self.data_values[nset]['wavelength'])
        inds = numpy.where(self.data_values[nset]['plot_mask'])
        wl1 = wl0[inds]
        fl0 = numpy.copy(self.data_values[nset]['flambda'])
        fl1 = fl0[inds]
        inds = numpy.argsort(wl1)
        wl2 = wl1[inds]
        fl2 = fl1[inds]
        flux1 = sed_utilities.integrate_sed(wl2, fl2)
        outstring = '\nData set: %d\nEstimated flux: %13.6e W/m^2\n' % (
            nset+1, flux1)
        try:
            value = float(parallax_entry.get())
            value_uncertainty = float(parallax_uncertainty_entry.get())
            option = parallax_variable.get()
            if option == 0:
                parallax = value
                distance = 1./parallax
                outstring = outstring + \
                    'Assumed parallax: %.4f mas (distance %.4f kpc)\n' % (
                        parallax, distance)
                if value_uncertainty <= 0.:
                    dmin = distance
                    dmax = distance
                else:
                    dmax = 1./(parallax - value_uncertainty)
                    dmin = 1./(parallax + value_uncertainty)
            else:
                distance = value
                parallax = 1./distance
                outstring = outstring + \
                    'Assumed distance: %.4f kpc (parallax %.4f mas)\n' % (
                        distance, parallax)
                dmin = distance - value_uncertainty
                dmax = distance + value_uncertainty
            parsec = 3.0856781e+16
            scale = 4.*parsec*parsec*math.pi*1.e+06*distance*distance
            xlsun = 3.84134e+26
            luminosity = flux1*scale/xlsun
            outstring = outstring + 'Estimated luminosity: %.4f L_sun\n' % (
                luminosity)
            if dmin < dmax:
                scale1 = 4.*parsec*parsec*math.pi*1.e+06*dmax*dmax
                scale2 = 4.*parsec*parsec*math.pi*1.e+06*dmin*dmin
                lum_max = flux1*scale1/xlsun
                lum_min = flux1*scale2/xlsun
                outstring = outstring + \
                    'Luminosity range: %.4f to %.4f L_sun' % (
                        lum_min,
                        lum_max
                        ) + '\n  (distance uncertainties only)\n'
            xmbol = 4.74 - 2.5*math.log10(luminosity)
            ambol = xmbol + 5.0*math.log10(distance/0.010)
            outstring = outstring + 'Bolometric magnitude: %.3f\n' % (ambol)
            outstring = outstring + \
                'Absolute bolometric magnitude: %.3f\n' % (xmbol)
            if dmin < dmax:
                xmbol1 = 4.74 - 2.5*math.log10(lum_max)
                xmbol2 = 4.74 - 2.5*math.log10(lum_min)
                delmbol1 = xmbol - xmbol1
                delmbol2 = xmbol2 - xmbol
                outstring = outstring + \
                    'Bolometric magnitude range -%.3f +%.3f\n' % (
                        delmbol1, delmbol2)
        except:
            pass
        outstring = outstring + ' \n'
        tkinter_utilities.append_text(text_area, outstring)

    def __toggle_data_status(self):
        """
        Routine to toggle the data status flag.

        Points with the flag value set to True are plotted.
        """
        self.flags['toggle_status'] = True

    def __toggle_mask_area(self):
        """
        Routine to set the mask area flag.  If set, cursor key events will 
        define an area within which to mask all points.
        """
        self.flags['mask_area'] = True
        
    def unmask_all(self):
        """
        Routine to unmask all data points for plotting.

        This undoes any points that have been masked by mouse input.
        """
        for loop in range(len(self.data_values)):
            for n1 in range(len(self.data_values[loop]['plot_mask'])):
                self.data_values[loop]['plot_mask'][n1] = True
        self.make_plot()

    def __toggle_zoom_status(self):
        """
        Routine to set the zoom status flag.

        When set, cursor button press/release is used to zoom on the plot.
        """
        self.flags['zoom_status'] = True

    def save_as_png(self):
        """
        Save the current plot as a PNG file.
        """
        outfile = tkinter.filedialog.asksaveasfilename(
            filetypes=[('PNG', '.png')])
        if isinstance(outfile, type('string')):
            s1 = outfile.split('.')
            if 'png' not in s1[-1]:
                outfile = outfile+'.png'
        self.tkinter_values['figure'].savefig(outfile, format="PNG")

    def save_as_postscript(self):
        """
        Save the current plot as a postscript file.
        """
        outfile = tkinter.filedialog.asksaveasfilename(
            filetypes=[('PS', '.ps')])
        if isinstance(outfile, type('string')):
            s1 = outfile.split('.')
            if 'ps' not in s1[-1]:
                outfile = outfile+'.ps'
        self.tkinter_values['figure'].savefig(outfile, format="PS")

    def __make_plot_area(self, parent):
        """
        Set up the main figure area for the plot.

        This routine makes the figure area and the sub-plot, and then
        sets up some event call-backs.

        Parameters
        ----------
            parent :   A Tk.Frame variable that holds the plot

        Returns
        -------
            No values are returned by this routine.

        """
        # The size of the plot area is determined here (values in inches).
        # The x/y ratio is 10 to 7 as in xmgrace.  One could also use the
        # golden ratio 1.618 to 1.  Take whatever values seem to be the
        # most aesthetically pleasing.  The DPI value should probably not
        # be set below 100.
        figure_area = Figure(figsize=(7.142857, 5.), dpi=100)
        figure_subplot = figure_area.add_subplot(1, 1, 1)
        self.tkinter_values['subplot'] = figure_subplot
        # The following sets up a label above the plot, for the plot
        # position.
        position_label_text = Tk.StringVar()
        self.tkinter_values['variables'][0] = position_label_text
        position_label = Tk.Label(parent, textvariable=position_label_text)
        position_label.pack(side=Tk.TOP)
        position_label_text.set("Position:")
        position_label.config(bg=BGCOL)
        main_canvas = FigureCanvasTkAgg(figure_area, master=parent)
        # Here are defined the events that the program responds to for
        # the plot area.
        main_canvas.mpl_connect("motion_notify_event", self.__plot_position)
        main_canvas.mpl_connect("key_press_event", self.__key_commands)
        main_canvas.mpl_connect("button_press_event", self.__plot_marker_set)
        main_canvas.mpl_connect("button_release_event",
                                self.__plot_marker_release)
        main_canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH,
                                         expand=Tk.YES)
        main_canvas.draw()
        self.tkinter_values['canvas'] = main_canvas
        self.tkinter_values['figure'] = figure_area
        self.make_plot()

    def get_bounds(self, subplot):
        """
        Find the plot bounds, and optionally set the range fields.

        If the self.flags['auto_scale'] flag is set the plot range is set to
        whatever the matplotlib range values are.  Otherwise the four range
        fields are read and these values are used to set the plot range.

        Parameters
        ----------

        subplot:   A matplotlib subplot variable

        Returns
        -------

        No values are returned
        """
        xmin, xmax = subplot.get_xbound()
        ymin, ymax = subplot.get_ybound()
        if self.flags['auto_scale']:
            tkinter_utilities.set_entry(self.tkinter_values['entry'][
                'xmin'], strfmt(xmin))
            tkinter_utilities.set_entry(self.tkinter_values['entry'][
                'xmax'], strfmt(xmax))
            tkinter_utilities.set_entry(self.tkinter_values['entry'][
                'ymin'], strfmt(ymin))
            tkinter_utilities.set_entry(self.tkinter_values['entry'][
                'ymax'], strfmt(ymax))
        else:
            xmin = float(self.tkinter_values['entry']['xmin'].get())
            xmax = float(self.tkinter_values['entry']['xmax'].get())
            ymin = float(self.tkinter_values['entry']['ymin'].get())
            ymax = float(self.tkinter_values['entry']['ymax'].get())
            subplot.set_xbound(xmin, xmax)
            subplot.set_ybound(ymin, ymax)
        return xmin, xmax, ymin, ymax

    def make_plot(self):
        """
        This is the main routine to make the SED plot.

        It reads the data sets and the control values to set some aspects of
        the plot.  Other aspects are fixed (i.e. the plot is always a log-log
        plot on the assumption that a large wavelength range is covered).
        """
        subplot = self.tkinter_values['subplot']
        canvas = self.tkinter_values['canvas']
        xmin, xmax, ymin, ymax = self.get_bounds(subplot)
        self.tkinter_values['canvas'].draw()
        xlabels = [r'Wavelength ($\mu$m)', 'Frequency (GHz)']
        ylabels = [r'$\lambda$ F$_\lambda$ (W/m$^2$)',
                   r'F$_\lambda$ (W/m$^2$/$\mu$m)',
                   r'F$_\nu$ (Jy)',
                   r'$\lambda^4$ F$_\lambda$ ($\mu$m$^3$W/m$^2$)']
        xprop = ['wavelength', 'frequency']
        yprop = ['lfl', 'flambda', 'fnu', 'l4fl']
        eyprop = ['dlfl', 'dflambda', 'dfnu', 'dl4fl']
        subplot.clear()
        if self.data_values[0]['wavelength'] is not None:
            ytype = self.tkinter_values['variables'][1].get()
            xtype = self.tkinter_values['variables'][2].get()
            for loop in range(len(self.data_values)):
                x_data_values, y_data_values, filterset = get_data_values(
                    self.data_values[loop][xprop[xtype]],
                    self.data_values[loop][yprop[ytype]],
                    self.data_values[loop]['plot_mask'],
                    self.data_values[loop]['filter_name'])
                dx_data_values, dy_data_values, dummy = get_data_values(
                    self.data_values[loop][xprop[xtype]],
                    self.data_values[loop][eyprop[ytype]],
                    self.data_values[loop]['plot_mask'],
                    self.data_values[loop]['filter_name'])
                dx_data_values = dx_data_values * 0.
                if self.data_values[loop]['plot']['show']:
                    if self.data_values[loop]['plot']['line']:
                        if self.data_values[loop]['plot']['symbol'] is None:
                            subplot.loglog(
                                x_data_values, y_data_values,
                                color=self.data_values[loop]['plot'][
                                    'colour'],
                                linestyle=self.data_values[loop]['plot'][
                                    'line_type'],
                                linewidth=self.data_values[loop]['plot'][
                                    'line_width'])
                            subplot.errorbar(
                                x_data_values, y_data_values, dy_data_values,
                                dx_data_values, fmt='None', elinewidth=1.0,
                                ecolor=self.data_values[loop]['plot'][
                                    'colour'])
                        else:
                            subplot.loglog(
                                x_data_values, y_data_values,
                                color=self.data_values[loop]['plot'][
                                    'colour'],
                                marker=self.data_values[loop]['plot'][
                                    'symbol'],
                                linestyle=self.data_values[loop]['plot'][
                                    'line_type'],
                                linewidth=self.data_values[loop]['plot'][
                                    'line_width'],
                                markersize=self.data_values[loop]['plot'][
                                    'size'])
                            subplot.errorbar(
                                x_data_values, y_data_values, dy_data_values,
                                dx_data_values, fmt='None', elinewidth=1.0,
                                ecolor=self.data_values[loop]['plot'][
                                    'colour'])
                    else:
                        if self.data_values[loop]['colour_by_name']:
                            self.plot_by_name(
                                x_data_values, y_data_values, filterset,
                                dy_data_values, dx_data_values,
                                subplot, loop)
                        else:
                            subplot.loglog(
                                x_data_values, y_data_values,
                                color=self.data_values[loop]['plot']['colour'],
                                marker=self.data_values[loop]['plot'][
                                    'symbol'],
                                linestyle='None',
                                markersize=self.data_values[loop]['plot'][
                                    'size'])
                            subplot.errorbar(
                                x_data_values, y_data_values, dy_data_values,
                                dx_data_values, fmt='None', elinewidth=1.0,
                                ecolor=self.data_values[loop]['plot'][
                                    'colour'])
            xmin, xmax, ymin, ymax = self.get_bounds(subplot)
            subplot.set_xlabel(xlabels[xtype], size=16)
            subplot.set_ylabel(ylabels[ytype], size=16)
            if (xtype == 1) and (not subplot.xaxis_inverted()):
                subplot.invert_xaxis()
            if (xtype != 1) and subplot.xaxis_inverted():
                subplot.invert_xaxis()
            subplot.tick_params(axis='x', direction='in')
            subplot.tick_params(axis='y', direction='in')
            subplot.tick_params(axis='x', direction='in', which='minor')
            subplot.tick_params(axis='y', direction='in', which='minor')
            subplot.tick_params(left=True, right=True, which='both')
            subplot.tick_params(top=True, bottom=True, which='both')
            subplot.tick_params(which='both', labelsize=14)
        self.tkinter_values['figure'].tight_layout()
        canvas.draw()

    def plot_by_name(self, x_data_values, y_data_values, filterset,
                     dy_data_values, dx_data_values, subplot, set_number):
        """
        For Vizier SED data, plot the plots with different colours or
        symbols by filter name.

        Parameters
        ----------

        x_data_values:   A numpy float array of the x data values to plot

        y_data_values:   A numpy float array of the y data values to plot

        filterset:       A numpy string array of the filter names for the
                         points

        dy_data_values:  A numpy array of the y uncertaintes to plot

        dx_data_values:  A numpy array of (dummy) x uncertainties to plot

        subplot:         A Tkinter subplot variable where the plot is made

        set_number:      The number for the set (0 to len(self.data_values)-1)

        Returns
        -------

        No values are returned by this routine.
        """
        manycolours = [
            'blue', 'forestgreen', 'black', 'orange', 'red', 'pink',
            'purple', 'cyan', 'lime', 'brown', 'violet', 'grey', 'gold',
            'maroon', 'moccasin', 'goldenrod', 'lightcyan', 'azure', 'navy',
            'indigo', 'beige', 'steelblue', 'darkseagreen', 'thistle'
        ]
        manysymbols = [
            "o", "s", "v", "^", "<", ">", "1", "2", "3", "4",
            "8", "p", "P", "*", "h", "H", "+", "x", "X", "D", "d"]
        if self.data_values[set_number]['source'] == 'vizier_sed':
            subfilterset = numpy.copy(filterset)
            for loop in range(len(subfilterset)):
                values = subfilterset[loop].split(':')
                if len(values[0]) > 0:
                    subfilterset[loop] = values[0]
                else:
                    subfilterset[loop] = 'monohromatic'
            allfilters = numpy.unique(subfilterset)
            ncount = 0
            for fname in allfilters:
                inds = numpy.where(subfilterset == fname)
                ncol = ncount % len(manycolours)
                nsymbol = ncount // len(manycolours)
                nsize = self.data_values[set_number]['plot']['size']
                subplot.loglog(
                    x_data_values[inds], y_data_values[inds],
                    color=manycolours[ncol],
                    marker=manysymbols[nsymbol],
                    linestyle='None',
                    markersize=nsize)
                subplot.errorbar(
                    x_data_values[inds], y_data_values[inds],
                    dy_data_values[inds],
                    dx_data_values[inds], fmt='None', elinewidth=1.0,
                    ecolor=manycolours[ncol])
                ncount = ncount + 1
                if ncount == len(manycolours)*len(manysymbols):
                    ncount = 0
                    nsize = nsize + 5
        else:
            subplot.loglog(
                x_data_values, y_data_values,
                color=self.data_values[set_number]['plot'][
                    'colour'],
                marker=self.data_values[set_number]['plot'][
                    'symbol'],
                linestyle=self.data_values[set_number]['plot'][
                    'line_type'],
                linewidth=self.data_values[set_number]['plot'][
                    'line_width'],
                markersize=self.data_values[set_number]['plot'][
                    'size'])
            subplot.errorbar(
                x_data_values, y_data_values, dy_data_values,
                dx_data_values, fmt='None', elinewidth=1.0,
                ecolor=self.data_values[set_number]['plot']['colour'])

    def __autoscale_plot(self):
        """
        This routine replots the figure with the auto scaling option.
        """
        self.flags['auto_scale'] = True
        self.make_plot()

    def freeze_range(self):
        """
        This routine replots the figure with the Entry values range.
        """
        self.flags['auto_scale'] = False
        self.make_plot()

    def change_set_properties(self):
        """
        Bring up a window to control the data set plotting symbol properties.
        """
        try:
            ndatasets = len(self.data_values)
            if (ndatasets == 0) or (self.data_values[0]['wavelength'] is None):
                return
        except ValueError:
            return
        window1 = Tk.Toplevel(self.root)
        window1.config(bg=BGCOL)
        window1.title('Data Set Display Properties')
        frame1 = Tk.Frame(window1)
        frame1.pack(side=Tk.TOP)
        frame1.config(bg=BGCOL)
        label1 = Tk.Label(frame1, text='Data set number:')
        label1.grid(column=0, row=0)
        set_number_menu = tkinter.ttk.Combobox(frame1, width=10)
        menu = []
        for loop in range(ndatasets):
            menu.append(str(loop+1))
        set_number_menu['values'] = menu
        set_number_menu.grid(column=1, row=0)
        set_number_menu.current(0)
        label1 = Tk.Label(frame1, text='Symbol:')
        label1.grid(column=0, row=1)
        symbol_menu = tkinter.ttk.Combobox(frame1, width=10)
        symbol_menu['values'] = MATPLOTLIB_SYMBOL_NAME_LIST
        symbol_menu.grid(column=1, row=1)
        symbol_menu.current(0)
        label1 = Tk.Label(frame1, text='Symbol Size:')
        label1.grid(column=0, row=2)
        size_entry = Tk.Entry(frame1, width=10)
        size_entry.grid(column=1, row=2)
        size_entry.insert(0, '2.0')
        label1 = Tk.Label(frame1, text='Colour:')
        label1.grid(column=0, row=3)
        colour_menu = tkinter.ttk.Combobox(frame1, width=10)
        colour_menu['values'] = COLOUR_SET
        colour_menu.grid(column=1, row=3)
        colour_menu.current(0)
        label1 = Tk.Label(frame1, text='Line:')
        label1.grid(column=0, row=4)
        line_variable = Tk.IntVar()
        button_frame = Tk.Frame(frame1)
        button_frame.grid(column=1, row=4)
        tkinter_utilities.put_yes_no(button_frame, line_variable,
                                     ['yes', 'no'], False)
        label1 = Tk.Label(frame1, text='Line Type:')
        label1.grid(column=0, row=5)
        line_type_menu = tkinter.ttk.Combobox(frame1, width=10)
        line_type_menu['values'] = MATPLOTLIB_LINE_NAME_LIST
        line_type_menu.grid(column=1, row=5)
        line_type_menu.current(0)
        label1 = Tk.Label(frame1, text='Line Width:')
        label1.grid(column=0, row=6)
        line_width_entry = Tk.Entry(frame1, width=10)
        line_width_entry.grid(column=1, row=6)
        line_width_entry.insert(0, '1.0')
        label1 = Tk.Label(frame1, text='Show set:')
        label1.grid(column=0, row=7)
        plot_variable = Tk.IntVar()
        button_frame = Tk.Frame(frame1)
        button_frame.grid(column=1, row=7)
        if self.data_values[0]['plot']['show']:
            tkinter_utilities.put_yes_no(button_frame, plot_variable,
                                         ['yes', 'no'], True)
        else:
            tkinter_utilities.put_yes_no(button_frame, plot_variable,
                                         ['yes', 'no'], False)
        label1 = Tk.Label(frame1, text='Set label:')
        label1.grid(column=0, row=8)
        source_entry = Tk.Entry(frame1, width=20)
        source_entry.grid(column=1, row=8)
        source_entry.insert(0, self.data_values[0]['source'])
#
        self.tkinter_values['set_parameters'] = [
            set_number_menu, symbol_menu, size_entry, line_variable,
            line_type_menu, line_width_entry, colour_menu, plot_variable,
            source_entry]
        set_number_menu.bind(
            "<<ComboboxSelected>>", self.update_set_fields)
        self.update_set_fields(None)
        frame2 = Tk.Frame(window1)
        frame2.pack(side=Tk.TOP)
        frame2.config(bg=BGCOL)
        apply_button = Tk.Button(
            frame2, text='Apply',
            command=lambda: self.apply_plot_values(
                self.tkinter_values['set_parameters']))
        apply_button.pack(side=Tk.LEFT)
        apply_button.config(bg=BGCOL)
        exit_button = Tk.Button(frame2, text='Close', command=window1.destroy)
        exit_button.pack(side=Tk.LEFT)
        exit_button.config(bg=BGCOL)

    def update_set_fields(self, event):
        """
        Update the set parameters on the parameter control window.

        Parameters
        ----------

        event:  A tkinter event variable

        Returns
        -------

        No values are returned.
        """
        parameters = self.tkinter_values['set_parameters']
        nset = parameters[0].current()
        for loop in range(len(MATPLOTLIB_SYMBOL_LIST)):
            if self.data_values[nset]['plot']['symbol'] == \
                    MATPLOTLIB_SYMBOL_LIST[loop]:
                parameters[1].current(loop)
        tkinter_utilities.set_entry(parameters[2], str(
            self.data_values[nset]['plot']['size']))
        if self.data_values[nset]['plot']['line']:
            parameters[3].set(1)
        else:
            parameters[3].set(0)
        if self.data_values[nset]['plot']['line_type'] is None:
            parameters[4].current(4)
        else:
            for loop in range(len(MATPLOTLIB_LINE_NAME_LIST)-1):
                if self.data_values[nset]['plot']['line_type'] == \
                        MATPLOTLIB_LINE_LIST[loop]:
                    parameters[4].current(loop)
        tkinter_utilities.set_entry(
            parameters[5], str(self.data_values[nset]['plot']['line_width']))
        for loop in range(len(COLOUR_SET)):
            if self.data_values[nset]['plot']['colour'] == \
                    COLOUR_SET[loop]:
                parameters[6].current(loop)
        if self.data_values[nset]['plot']['show']:
            parameters[7].set(1)
        else:
            parameters[7].set(0)
        tkinter_utilities.set_entry(
            parameters[8], self.data_values[nset]['source'])

    def apply_plot_values(self, parameters):
        """
        Read the data set parameters and replot the figure.

        Parameters
        ----------
            parameters:   A list of 8 Tkinter objects in the parent window --
                          a mix of menu variables, radio button variables,
                          and entry field variables

        Returns
        -------
            Nothing
        """
        nset = parameters[0].current()
        nsymbol = parameters[1].current()
        symbol_size = float(parameters[2].get())
        lineflag = parameters[3].get()
        if lineflag == 1:
            lineflag = True
        else:
            lineflag = False
        nline = parameters[4].current()
        line_width = float(parameters[5].get())
        ncolour = parameters[6].current()
        plotflag = parameters[7].get()
        if plotflag == 1:
            plotflag = True
        else:
            plotflag = False
        if (symbol_size <= 0.) or (line_width <= 0.):
            tkinter.messagebox.showinfo(
                'Error', 'Bad input values.  Check the settings.')
            return
        self.data_values[nset]['plot']['symbol'] = \
            MATPLOTLIB_SYMBOL_LIST[nsymbol]
        self.data_values[nset]['plot']['size'] = symbol_size
        self.data_values[nset]['plot']['line'] = lineflag
        self.data_values[nset]['plot']['line_type'] = \
            MATPLOTLIB_LINE_LIST[nline]
        self.data_values[nset]['plot']['line_width'] = line_width
        self.data_values[nset]['plot']['colour'] = COLOUR_SET[ncolour]
        self.data_values[nset]['plot']['show'] = plotflag
        self.make_plot()

    def read_spectrum_file(self, filename=None, model=False):
        """
        Utility to bring up a window for reading in a spectrum.

        One can also read in photometry values if these are in flux density
        units in an input file, as long as there are two columns to read with
        the wavelengths/frequency values and the flux density values in units
        that the code can handle.

        Parameters
        ----------

        filename:    An optional string variable giving the file name.  This
                     string is put into the name field.

        model:   An optional boolean value for whether the spectrum read in
                 is taken to be a model (this affects only the selection of
                 data sets as opposed to model sets).  The default is False
                 in which case the spectrum is assumed to be an observed
                 spectrum.
        """
        window1 = Tk.Toplevel(self.root)
        window1.config(bg=BGCOL)
        frame1 = Tk.Frame(window1)
        frame1.pack(side=Tk.TOP)
        frame1.config(bg=BGCOL)
        name_field = Tk.Entry(frame1, width=40)
        name_field.grid(column=1, row=0)
        name_field.config(bg=BGCOL)
        if filename is None:
            name_field.insert(0, ' ')
        else:
            name_field.insert(0, filename)
        button1 = Tk.Button(
            frame1, text="Select File",
            command=lambda: self.select_file(name_field))
        button1.grid(column=0, row=0)
        button1.config(bg=BGCOL)
        label1 = Tk.Label(frame1, text='x column/label:')
        label1.grid(column=0, row=1)
        label1.config(bg=BGCOL)
        xcolumn_entry = Tk.Entry(frame1, width=20)
        xcolumn_entry.grid(column=1, row=1)
        xcolumn_entry.config(bg=BGCOL)
        xcolumn_entry.insert(0, '1')
        label1 = Tk.Label(frame1, text='y column/label:')
        label1.grid(column=0, row=2)
        label1.config(bg=BGCOL)
        ycolumn_entry = Tk.Entry(frame1, width=20)
        ycolumn_entry.grid(column=1, row=2)
        ycolumn_entry.config(bg=BGCOL)
        ycolumn_entry.insert(0, '2')
        label1 = Tk.Label(frame1, text='x units:')
        label1.grid(column=0, row=3)
        label1.config(bg=BGCOL)
        unit_menu1 = tkinter.ttk.Combobox(frame1, width=10)
        unit_menu1['values'] = ['microns', 'Angstroms', 'GHz', 'Hz']
        unit_menu1.grid(column=1, row=3)
        unit_menu1.current(0)
        label1 = Tk.Label(frame1, text='y units:')
        label1.grid(column=0, row=4)
        label1.config(bg=BGCOL)
        unit_menu2 = tkinter.ttk.Combobox(frame1, width=10)
        unit_menu2['values'] = ['W/m^2/micron', 'erg/s/cm^2/Angstrom',
                                'W/m^2', 'erg/s/cm^2', 'W/m^2/Hz', 'Jy']
        unit_menu2.current(0)
        unit_menu2.grid(column=1, row=4)
        frame2 = Tk.Frame(window1)
        frame2.pack(side=Tk.TOP)
        frame2.config(bg=BGCOL)
        parameters = [name_field, xcolumn_entry, ycolumn_entry,
                      unit_menu1, unit_menu2]
        apply_button = Tk.Button(
            frame2, text='Apply',
            command=lambda: self.apply_spectrum_options(parameters, model))
        apply_button.pack(side=Tk.LEFT)
        apply_button.config(bg=BGCOL)
        exit_button = Tk.Button(frame2, text='Close', command=window1.destroy)
        exit_button.pack(side=Tk.LEFT)
        exit_button.config(bg=BGCOL)

    def select_file(self, entry):
        """
        Ask for a file name and populate the entry field with the name.

        Parameters
        ----------
            entry:   The name of a Tkinter Entry variable which will be
                     populated with a file name

        Returns
        -------
            Nothing
        """
        name = tkinter.filedialog.askopenfilename()
        string1 = entry.get()
        entry.delete(0, len(string1))
        entry.insert(0, name)

    def apply_spectrum_options(self, parameters, model=False, filename=None):
        """
        Take the spectrum read-in window options and read in the data.

        This is also used for generic model spectrum read-in.

        Parameters
        ----------
            parameters:  The Tkinter variables for the five input fields in the
                         spectrum reading window.

            model:       An optional boolean flag for whether the spectrum is
                         assumed to be a model rather than data values.

        Returns
        -------
            Nothing

        If the spectrum is read in, the self.data_values list is updated with
        a new data set.  If it worked, the plot is regenerated.
        """
        try:
            filename = parameters[0].get()
            x_column_value = parameters[1].get()
            y_column_value = parameters[2].get()
            if '.fits' not in filename[-5:]:
                x_column_value = int(x_column_value)
                y_column_value = int(y_column_value)
                if x_column_value < 1:
                    raise ValueError
                if y_column_value < 1:
                    raise ValueError
            try:
                wavelengths, spectrum = model_utilities.read_generic_model(
                    filename, x_column_value, y_column_value)
                if (len(wavelengths) != len(spectrum)) or \
                   (numpy.min(wavelengths) <= 0.):
                    raise ValueError
            except:
                raise ValueError
            opt1 = parameters[3].current()
            opt2 = parameters[4].current()
            values = filename.split('/')
            if (not model) and (wavelengths is None):
                outstring = 'Error: file %s could not be read.\n' % (
                    values[-1])
                tkinter.messagebox.showinfo('Information', outstring)
                return
            if model:
                label = 'stellar model ' + values[-1]
            else:
                label = 'spectrum: ' + values[-1]
            count1 = len(self.data_values)
            if self.data_values[0]['wavelength'] is None:
                count1 = 0
            self.add_spectrum(wavelengths, spectrum, opt1, opt2, label)
            count2 = len(self.data_values)
            if self.data_values[0]['wavelength'] is None:
                count2 = 0
            if model and (count2 == count1+1):
                wavelengths = numpy.copy(
                    self.data_values[count1]['wavelength'])
                spectrum = numpy.copy(self.data_values[count1]['flambda'])
                wl0 = WAVELNORM
                fl0 = FLAMBDANORM
                fl1 = numpy.interp(wl0, wavelengths, spectrum)
                if fl1 > 0.:
                    scale = fl0/fl1
                    self.data_values[count1]['flambda'] = \
                        self.data_values[count1]['flambda']*scale
                    self.data_values[count1]['lfl'] = \
                        self.data_values[count1]['lfl']*scale
                    self.data_values[count1]['l4fl'] = \
                        self.data_values[count1]['l4fl']*scale
                    self.data_values[count1]['fnu'] = \
                        self.data_values[count1]['fnu']*scale
            self.make_plot()
        except ValueError:
            return

    def add_spectrum(self, xdata, ydata, opt1, opt2, label):
        """

        Parameters
        ----------
        xdata:  a numpy float aray of wavelength or frequency values

        ydata:  a numpy float array f flux density values

        opt1:   an integer flag for the xdata type

        opt2:   an integer flag for the ydata type 

        label:  a strong variable, the label for the data set 

        Returns
        -------
        None.

        """
        data_set = {'wavelength': None, 'frequency': None, 'fnu': None,
                    'flambda': None, 'lfl': None, 'l4fl': None,
                    'dfnu': None, 'dflambda': None, 'dl4fl': None,
                    'dlfl': None, 'mask': None, 'plot_mask': None,
                    'filter_name': None, 'distance': None,
                    'position': None, 'refpos': None,
                    'references': None, 'plot': None, 'source': None,
                    'colour_by_name': False}
        data_set['wavelength'] = sed_utilities.trans_wavelength(
            xdata, opt1, 0)
        data_set['frequency'] = sed_utilities.trans_wavelength(
            xdata, opt1, 2)
        data_set['flambda'] = sed_utilities.trans_flux_density(
            data_set['wavelength'], ydata, opt2, 0)
        data_set['fnu'] = sed_utilities.trans_flux_density(
            data_set['wavelength'], ydata, opt2, 5)
        data_set['lfl'] = sed_utilities.trans_flux_density(
            data_set['wavelength'], ydata, opt2, 2)
        data_set['l4fl'] = data_set['lfl']*(data_set['wavelength']**3)
        data_set['dwavelength'] = data_set['wavelength'] * 0.
        data_set['dfrequency'] = data_set['wavelength'] * 0.
        data_set['dflambda'] = data_set['wavelength'] * 0.
        data_set['dfnu'] = data_set['wavelength'] * 0.
        data_set['dlfl'] = data_set['wavelength'] * 0.
        data_set['dl4fl'] = data_set['wavelength'] * 0.
        mask = []
        for loop in range(len(data_set['wavelength'])):
            mask.append('True')
        data_set['plot_mask'] = numpy.asarray(mask)
        data_set['source'] = label
        self.add_data_set(data_set, False)

    def read_vizier_file(self):
        """
        A wrapper routine for reading in a Vizier photometry VOT file.
        """
        try:
            filename = tkinter.filedialog.askopenfilename()
            votdata = read_vizier_sed_table.read_vizier_vot(filename)
            if votdata[0] is not None:
                data_set = unpack_vot(votdata)
                data_set['source'] = 'vizier_sed'
                self.add_data_set(data_set, True)
                self.make_plot()
        except:
            pass

    def add_data_set(self, data_set, marker_flag):
        """
        A utility routine to add a data set to the object list.

        Parameters
        ----------
            data_set:   A dictionary variable with the data set properties

            marker_flag:   A boolean value for whether the data values should
                           be plotted as discrete symbols

        Returns
        -------
            Nothing
        """
        ndata = len(self.data_values)
        if (ndata == 1) and (
                self.data_values[0]['wavelength'] is None):
            data_set['plot'] = {'symbol': MATPLOTLIB_SYMBOL_LIST[1],
                                'size': 2.0, 'line': False,
                                'line_type': MATPLOTLIB_LINE_LIST[0],
                                'line_width': 1.0,
                                'colour': COLOUR_SET[0],
                                'show': True}
            if not marker_flag:
                data_set['plot']['symbol'] = None
                data_set['plot']['line'] = True
            self.data_values[0] = data_set
        else:
            ncolour = ndata % len(COLOUR_SET)
            nsymbol = ndata % (len(MATPLOTLIB_SYMBOL_LIST) - 1)
            nsymbol = nsymbol + 1
            data_set['plot'] = {
                'symbol': MATPLOTLIB_SYMBOL_LIST[nsymbol],
                'size': 2.0, 'line': False,
                'line_type': 'solid',
                'line_width': 1.0,
                'colour': COLOUR_SET[ncolour],
                'show': True
            }
            if not marker_flag:
                data_set['plot']['symbol'] = None
                data_set['plot']['line'] = True
            self.data_values.append(data_set)

    def stellar_model_window(self):
        """
        This routine makes the window for reading in stellar model values.

        Parameters
        ----------
        None.

        Returns
        -------
        Nothing.

        """
        window1 = Tk.Toplevel(self.root)
        window1.config(bg=BGCOL)
        window1.title('Stellar Model Read-In')
        text_area = ScrolledText(window1, height=6, width=50, wrap=Tk.NONE)
        text_area.pack(side=Tk.TOP)
        text_area.config(font=('courier', 16, 'bold'))
        frame1 = Tk.Frame(window1)
        frame1.pack(side=Tk.TOP)
        frame1.config(bg=BGCOL)
        filename_entry = Tk.Entry(frame1, width=40)
        filename_entry.grid(row=0, column=1)
        button1 = Tk.Button(frame1, text='Select File',
                            command=lambda: self.get_model_name(
                                filename_entry, path_entry))
        button1.grid(row=0, column=0)
        lab1 = Tk.Label(frame1,text=' file path:')
        lab1.grid(row=1, column=0)
        path_entry = Tk.Entry(frame1, width=40)
        path_entry.grid(row=1, column=1)
        modeltype = Tk.IntVar()
        labels = ['BOSZ model', 'Allard Phoenix model',
                  'Husser Phoenix Grid', 'Calspec', 'Generic']
        modeltype = Tk.IntVar()
        for loop in range(len(labels)):
            b1 = Tk.Radiobutton(frame1, text=labels[loop], variable=modeltype,
                                value=loop)
            b1.config(bg=BGCOL)
            b1.grid(row=loop+2, column=0, columnspan=2)
        modeltype.set(0)
        self.tkinter_values['modeltype'] = modeltype
        resolution = Tk.IntVar()
        frame1 = Tk.Frame(window1)
        frame1.pack(side=Tk.TOP)
        sl1 = Tk.Scale(frame1, orient=Tk.HORIZONTAL, length=301,
                       from_=90, to=3010., resolution=10,
                       variable=resolution, label="Smoothing Resolution")
        sl1.config(bg=BGCOL)
        resolution.set(3010)
        sl1.pack(side=Tk.TOP)
        frame2 = Tk.Frame(window1)
        frame2.pack(side=Tk.TOP)
        frame2.config(bg=BGCOL)
        apply_button = Tk.Button(
            frame2, text='Apply',
            command=lambda: self.read_stellar_model(
                filename_entry, path_entry,
                self.tkinter_values['modeltype'],
                text_area, resolution)
            )
        apply_button.pack(side=Tk.LEFT)
        apply_button.config(bg=BGCOL)
        exit_button = Tk.Button(frame2, text='Close', command=window1.destroy)
        exit_button.pack(side=Tk.LEFT)
        exit_button.config(bg=BGCOL)

    def get_model_name(self, filename_entry, path_entry):
        """
        Routine to query for an input file name for a stellar model and put
        the name into an Entry field.

        Parameters
        ----------
        filename_entry:    A Tkinter Entry variable, will hold the input
                           file name.

        path_entry:        A Tkinter Entry variable, with hold the path 
                           to the file name

        Returns
        -------
        None.

        """
        try:
            filename = tkinter.filedialog.askopenfilename()
            values = filename.split('/')
            tkinter_utilities.set_entry(filename_entry, values[-1])
            pathname = filename.replace(values[-1], '')
            tkinter_utilities.set_entry(path_entry, pathname)
        except:
            pass

    def read_stellar_model(self, file_entry, path_entry, model_type,
                           text_area, resolution):
        """
        Read in a spectrum from a stellar model file.

        Parameters
        ----------

        file_entry:   A Tkinter entry field variable from which the file name
                      is read.

        path_entry:   A Tkinter entry field variable from which the file name
                      directory path is read.

        model_type:   An integer flag for the type of stellar model

        text_area:   A Tkinter text area variable to which messages are
                     written.

        resolution:  A TKinter variable giving the smoothing factor, only
                     used for the Phoenix models.

        Returns
        -------

        No values are returned.  A new data set is added if the model is
        read in successfully.
        """
        try:
            name = file_entry.get()
            path = path_entry.get()
            if path[-1] != '/':
                path = path+'/'
            modeltype = model_type.get()
            if not os.path.isfile(path+name):
                outstring = 'Error: file %s was not found.\n' % (name)
                tkinter_utilities.append_text(text_area, outstring)
                return
            if modeltype == 0:
                wavelengths, spectrum = model_utilities.read_bosz_model(name)
            elif modeltype == 1:
                r1 = resolution.get()
                wavelengths, spectrum = model_utilities.read_phoenix_model(
                    name, r1)
            elif modeltype == 2:
                r1 = resolution.get()
                wavelengths, spectrum = \
                    model_utilities.read_phoenix_grid_model(name, r1)
            elif modeltype == 3:
                wavelengths, spectrum = \
                    model_utilities.read_calspec_model(name)
            else:
                self.read_spectrum_file(name, True)
            if modeltype < 4:
                values = name.split('/')
                if wavelengths is None:
                    outstring = 'Error: file %s could not be read.\n' % (
                        values[-1])
                    tkinter_utilities.append_text(text_area, outstring)
                    return
                outstring = 'Have read in %d points from file %s\n' % (
                    len(wavelengths), values[-1])
                outstring = outstring + \
                    'Wavelength range: %f to %f microns\n' % (
                        wavelengths[0], wavelengths[-1])
                outstring = outstring + \
                    'Flux density range: %.6e to %.6e W/m^2/micron\n' % (
                        numpy.min(spectrum), numpy.max(spectrum))
                label = 'stellar model ' + values[-1]
                self.add_spectrum(wavelengths, spectrum, 0, 0, label)
                tkinter_utilities.append_text(text_area, outstring)
            self.make_plot()
        except:
            outstring = 'Error: file %s could not be read.\n' % (name)
            tkinter_utilities.append_text(text_area, outstring)

    def __scale_model_to_cursor(self):
        self.flags['model_scale'] = True

    def __scale_model_to_point(self):
        self.flags['point_scale'] = True

    def __fit_model_to_data(self):
        """
        Driver routien for matching a stellar model to a data set with the 
        best least-squares scaling.
        """
        label_names = [
            "2MASS:H", "2MASS:J", "2MASS:Ks", "Gaia:G", "Gaia:Gbp",
            "Gaia:Grp", "Bessell:B", "Bessell:I", "Bessell:R",
            "Bessell:U", "Bessell:V", "Cousins:I", "Cousins:R",
            "Johnson:B", "Johnson:I", "Johnson:R", "Johnson:U",
            "Johnson:V", "Stromgren:b", "Stromgren:u",
            "Stromgren:v", "Stromgren:y", "SDSS:g", "SDSS:g'",
            "SDSS:i", "SDSS:i'", "SDSS:r", "SDSS:r", "SDSS:u",
            "SDSS:u'", "SDSS:z", "SDSS:z'", "UKIDSS:H",
            "UKIDSS:J", "UKIDSS:K", "UKIDSS:Y", "UKIDSS:Z",
            "WISE:W1", "WISE:W2", "WISE:W3", "WISE:W4"]
        ndata = []
        nmodel = []
        for loop in range(len(self.data_values)):
            str1 = self.data_values[loop]['source']
            if 'stellar model' in str1[0:13]:
                nmodel.append(loop+1)
            if ('spectrum' in str1[0:8]) or ('vizier' in str1[0:6]):
                ndata.append(loop+1)
        if (len(ndata) == 0) or (len(nmodel) == 0):
            return
        for loop in ndata:
            frag = 'set %d' % (loop)
            nchar = len(frag)
            if frag in str1[0:nchar]:
                ndata.append(loop+1)
        for loop in nmodel:
            frag = 'model %d' % (loop)
            nchar = len(frag)
            if frag in str1[0:nchar]:
                nmodel.append(loop+1)
        if (len(ndata) > 1) or (len(nmodel) > 1):
            ndata, nmodel = self.__get_set_numbers(ndata, nmodel)
        datawl = []
        datafl = []
        dataerr = []
        labels = []
        if self.tkinter_values['variables'][4].get() == 1:
            for loop in range(len(ndata)):
                d1 = ndata[loop] - 1
                for ind in range(len(self.data_values[d1]['wavelength'])):
                    if self.data_values[d1]['plot_mask'][ind]:
                        datawl.append(self.data_values[d1]['wavelength'][ind])
                        datafl.append(self.data_values[d1]['flambda'][ind])
                        dataerr.append(self.data_values[d1]['dflambda'][ind])
                        labels.append(self.data_values[d1]['filter_name'][ind])
                datawl = numpy.asarray(datawl)
                datafl = numpy.asarray(datafl)
                dataerr = numpy.asarray(dataerr)
            weights = datawl*0. + 1.
            photon_mean_flag = False
        else:
            for loop in range(len(ndata)):
                d1 = ndata[loop] - 1
                for ind in range(len(self.data_values[d1]['wavelength'])):
                    if self.data_values[d1]['plot_mask'][ind]:
                        if self.data_values[d1]['filter_name'][ind] in \
                                label_names:
                            datawl.append(
                                self.data_values[d1]['wavelength'][ind])
                            datafl.append(
                                self.data_values[d1]['flambda'][ind])
                            dataerr.append(
                                self.data_values[d1]['dflambda'][ind])
                            labels.append(
                                self.data_values[d1]['filter_name'][ind])
                datawl = numpy.asarray(datawl)
                datafl = numpy.asarray(datafl)
                dataerr = numpy.asarray(dataerr)
            weights = datawl*0. + 1.
            photon_mean_flag = True
        for m1 in nmodel:
            scale, stats = model_utilities.match_data(
                datawl, datafl, dataerr, labels, weights,
                self.data_values[m1-1]['wavelength'],
                self.data_values[m1-1]['flambda'],
                True, photon_mean_flag, self.flags['use_ratio'])
            if scale > 0.:
                self.data_values[m1-1]['flambda'] = \
                    self.data_values[m1-1]['flambda']*scale
                self.data_values[m1-1]['lfl'] = \
                    self.data_values[m1-1]['lfl']*scale
                self.data_values[m1-1]['l4fl'] = \
                    self.data_values[m1-1]['l4fl']*scale
                self.data_values[m1-1]['fnu'] = \
                    self.data_values[m1-1]['fnu']*scale
                str1 = 'model %d with extinction' % (m1)
                for loop in range(len(self.data_values)):
                    if str1 in self.data_values:
                        self.data_values[loop]['flambda'] = \
                            self.data_values[loop]['flambda']*scale
                        self.data_values[loop]['lfl'] = \
                            self.data_values[loop]['lfl']*scale
                        self.data_values[loop]['l4fl'] = \
                            self.data_values[loop]['l4fl']*scale
                        self.data_values[loop]['fnu'] = \
                            self.data_values[loop]['fnu']*scale
                self.make_plot()

    def __get_set_numbers(self, ndata, nmodel):
        return ndata[0], nmodel[0]

    def __make_set_selection_window(self, ndata, nmodel):
        window1 = Tk.Toplevel(self.root)
        window1.config(bg=BGCOL)
        window1.title('Set Selection')
        text_area = ScrolledText(window1, height=6, width=50, wrap=Tk.NONE)
        text_area.pack(side=Tk.TOP)
        text_area.config(font=('courier', 16, 'bold'))
        frame1 = Tk.Frame(window1)
        frame1.pack(side=Tk.TOP)
        frame1.config(bg=BGCOL)
        label1 = Tk.Label(frame1, text='Data Sets to Fit:')
        label1.grid(row=0, column=0)
        label1.config(bg=BGCOL)
        field1 = tkinter.Entry(frame1, width=10)
        field1.grid(row=0, column=1)
        field1.config(bg=BGCOL)
        str1 = tkinter_utilities.list_sets_for_entry(ndata)
        field1.insert(0, str1)
        label2 = Tk.Label(frame1, text='Model to Scale:')
        label2.grid(row=1, column=0)
        label2.config(bg=BGCOL)
        field2 = tkinter.Entry(frame1, width=10)
        field2.grid(row=1, column=1)
        field2.config(bg=BGCOL)
        str1 = '%d' % (nmodel[0])
        field2.insert(0, str1)
        frame2 = Tk.Frame(window1)
        frame2.pack(side=Tk.TOP)
        frame2.config(bg=BGCOL)
        list_button = Tk.Button(frame2, text="List Sets",
                                command=lambda: self.list_sets(
                                    ndata, nmodel, text_area))
        list_button.pack(side=Tk.LEFT)
        list_button.config(bg=BGCOL)
        dlist = []
        modelout = -1
        apply_button = Tk.Button(
            frame2, text='Apply',
            command=lambda: self.read_set_fit_values(field1, field2,
                                                     dlist, modelout))
        apply_button.pack(side=Tk.LEFT)
        apply_button.config(bg=BGCOL)
        exit_button = Tk.Button(frame2, text='Close', command=window1.close)
        exit_button.pack(side=Tk.LEFT)
        exit_button.config(bg=BGCOL)

    def __plot_position(self, event):
        """
        Post the cursor position to the information line on the plot display.

        Routine to post the cursor position to the text area above the
        plot area.

        Parameters
        ----------
            event :   a motion-notify event from the image display window

        Returns
        -------
            No values are returned by this routine.

        """
        try:
            event.canvas.get_tk_widget().focus_set()
            if (event.xdata is not None) and (event.ydata is not None):
                ndata, npoints, dmin, xvalue, yvalue = self.__match_point(
                    event.xdata, event.ydata)
                outstring = "Position: x = %.6g y = %.6g" % (
                    event.xdata, event.ydata)
                if self.flags['toggle_status'] is True:
                    outstring = outstring + '  [delete/undelete active]'
                if ndata is not None:
                    str1 = "\nNearest point: x = %.6g " % (xvalue)
                    str1 = str1 + " y = %.6g [set %d, point %d]" % (
                        yvalue, ndata+1, npoints[0]+1)
                    if 'vizier_sed' in self.data_values[ndata]['source']:
                        str2 = ' (%s)]' % (
                            self.data_values[ndata]['filter_name'][npoints[0]])
                        str1 = str1.replace(']', str2)
                    outstring = outstring + str1
                    if not self.data_values[ndata]['plot_mask'][npoints[0]]:
                        outstring = outstring + ' (masked)'
                self.tkinter_values['variables'][0].set(outstring)
        except:
            pass

    def __key_commands(self, event):
        """
        A matplotlib callback for key press events.

        Parameters
        ----------

            event:  A matplotlib key press event variable

        Returns
        -------

            Nothing
        """
        if self.flags['model_scale'] or self.flags['point_scale']:
            self.match_model(event)
            self.flags['model_scale'] = False
            self.flags['point_scale'] = False
            return
        if (event.key == 'x') or (event.key == 'q'):
            self.flags['toggle_status'] = False
            self.make_plot()
        if event.key == 'd':
            if self.flags['toggle_status']:
                ndata, npoints, dmin, xpos, ypos = \
                    self.__match_point(event.xdata, event.ydata)
                if ndata is not None:
                    try:
                        for loop in range(len(npoints)):
                            self.data_values[ndata][
                                'plot_mask'][npoints[loop]] = False
                        self.make_plot()
                    except ValueError:
                        pass
        if event.key == 'u':
            if self.flags['toggle_status']:
                ndata, npoints, dmin, xpos, ypos = \
                    self.__match_point(event.xdata, event.ydata)
                if ndata is not None:
                    try:
                        for loop in range(len(npoints)):
                            self.data_values[ndata]['plot_mask'][
                                npoints[loop]] = True
                        self.make_plot()
                    except ValueError:
                        pass

    def match_model(self, event):
        """
        Code to scale a stellar model to a data set.

        Parameters
        ----------
        event:      A matplotlib Event variable, gives the position in the
                    plot for matching.

        Returns
        -------
        No values are returned.

        """
        set_numbers = []
        for loop in range(len(self.data_values)):
            if ('stellar' in self.data_values[loop]['source'][0:7]) or \
               ('model' in self.data_values[loop]['source'][0:6]):
                set_numbers.append(loop+1)
        if len(set_numbers) == 1:
            set_number = set_numbers[0] - 1
        else:
            set_number = None
            outstring = 'Select model to scale: ['
            for n1 in range(len(set_numbers)):
                if n1 < len(set_numbers) - 1:
                    outstring = outstring + '%d' % (set_numbers[n1]) + ', '
                else:
                    outstring = outstring + '%d' % (set_numbers[n1]) + ']'
            value = tkinter.simpledialog.askstring('Input', outstring)
            if value is not None:
                try:
                    i1 = int(value)
                    if i1 in set_numbers:
                        set_number = i1
                except:
                    pass
            if set_number is None:
                return
            self.flags['model_scale'] = False
            self.flags['point_scale'] = False
            set_number = set_number - 1
        xtype = self.tkinter_values['variables'][2].get()
        ytype = self.tkinter_values['variables'][1].get()
        scale = 1.0
        if xtype == 0:
            xvalues = self.data_values[set_number]['wavelength']
        else:
            xvalues = self.data_values[set_number]['frequency']
        if ytype == 1:
            yvalues = self.data_values[set_number]['flambda']
        elif ytype == 2:
            yvalues = self.data_values[set_number]['fnu']
        elif ytype == 3:
            yvalues = self.data_values[set_number]['l4fl']
        else:
            yvalues = self.data_values[set_number]['lfl']
        if self.flags['model_scale']:
            fl0 = numpy.interp(event.xdata, xvalues, yvalues)
            scale = event.ydata/fl0
        else:
            instring = self.tkinter_values['variables'][0].get()
            if 'Nearest' not in instring:
                xpoint = event.xdata
                ypoint = event.ydata
            else:
                values = instring.split()
                numbers = []
                for n1 in range(len(values)):
                    try:
                        a = float(values[n1])
                        numbers.append(a)
                    except:
                        pass
                xpoint = numbers[2]
                ypoint = numbers[3]
            fl0 = numpy.interp(xpoint, xvalues, yvalues)
            scale = ypoint/fl0
        self.data_values[set_number]['flambda'] = \
            self.data_values[set_number]['flambda']*scale
        self.data_values[set_number]['lfl'] = \
            self.data_values[set_number]['lfl']*scale
        self.data_values[set_number]['l4fl'] = \
            self.data_values[set_number]['l4fl']*scale
        self.data_values[set_number]['fnu'] = \
            self.data_values[set_number]['fnu']*scale
        # Check and see if there is a linked model with/without extinction,
        # and scale that model the same way.
        if 'with extinction' in self.data_values[set_number]['source']:
            values = self.data_values[set_number]['source'].split()
            other_number = int(values[1]) - 1
            self.data_values[other_number]['flambda'] = \
                self.data_values[other_number]['flambda']*scale
            self.data_values[other_number]['lfl'] = \
                self.data_values[other_number]['lfl']*scale
            self.data_values[other_number]['l4fl'] = \
                self.data_values[other_number]['l4fl']*scale
            self.data_values[other_number]['fnu'] = \
                self.data_values[other_number]['fnu']*scale
        else:
            str1 = 'model %d' % (set_number+1)
            for n1 in range(len(self.data_values)):
                if str1 in self.data_values[n1]['source']:
                    self.data_values[n1]['flambda'] = \
                        self.data_values[n1]['flambda']*scale
                    self.data_values[n1]['lfl'] = \
                        self.data_values[n1]['lfl']*scale
                    self.data_values[n1]['l4fl'] = \
                        self.data_values[n1]['l4fl']*scale
                    self.data_values[n1]['fnu'] = \
                        self.data_values[n1]['fnu']*scale
        self.make_plot()

    def __plot_marker_set(self, event):
        """
        A matplotlib callback for mouse button press events.

        Parameters
        ----------

            event:  A matplotlib mouse button press event variable

        Returns
        -------

            Nothing
        """
        if self.flags['model_scale'] or self.flags['point_scale']:
            self.match_model(event)
            return
        if self.flags['zoom_status']:
            self.flags['positions'].append(event.xdata)
            self.flags['positions'].append(event.ydata)
            return
        if self.flags['mask_area']:
            self.flags['positions'].append(event.xdata)
            self.flags['positions'].append(event.ydata)

    def __plot_marker_release(self, event):
        """
        A matplotlib callback for mouse button release events.

        Parameters
        ----------

            event:  A matplotlib mouse button release event variable

        Returns
        -------

            Nothing
        """
        if self.flags['zoom_status']:
            self.flags['positions'].append(event.xdata)
            self.flags['positions'].append(event.ydata)
            self.flags['zoom_status'] = False
            xmin, xmax, ymin, ymax = sed_utilities.sort_positions(
                self.flags['positions'])
            if xmin is not None:
                tkinter_utilities.set_entry(self.tkinter_values['entry'][
                    'xmin'], strfmt(xmin))
                tkinter_utilities.set_entry(self.tkinter_values['entry'][
                    'xmax'], strfmt(xmax))
                tkinter_utilities.set_entry(self.tkinter_values['entry'][
                    'ymin'], strfmt(ymin))
                tkinter_utilities.set_entry(self.tkinter_values['entry'][
                    'ymax'], strfmt(ymax))
                self.flags['auto_scale'] = False
                self.make_plot()
            self.flags['positions'] = []
            return
        if self.flags['mask_area']:
            self.flags['positions'].append(event.xdata)
            self.flags['positions'].append(event.ydata)
            self.flags['mask_area'] = False
            xmin, xmax, ymin, ymax = sed_utilities.sort_positions(
                self.flags['positions'])
            if xmin is not None:
                xprop = ['wavelength', 'frequency']
                yprop = ['lfl', 'flambda', 'fnu', 'l4fl']
                ytype = self.tkinter_values['variables'][1].get()
                xtype = self.tkinter_values['variables'][2].get()
                for loop in range(len(self.data_values)):
                   for n1 in range(len(self.data_values[loop]['wavelength'])):
                       if (self.data_values[loop][xprop[xtype]][n1] >= xmin) \
                          and \
                          (self.data_values[loop][xprop[xtype]][n1] <= xmax) \
                          and \
                          (self.data_values[loop][yprop[ytype]][n1] >= ymin) \
                          and \
                          (self.data_values[loop][yprop[ytype]][n1] <= ymax):
                           self.data_values[loop]['plot_mask'][n1] = not \
                               self.data_values[loop]['plot_mask'][n1]
                self.make_plot()
            self.flags['positions'] = []
            return
        if self.flags['toggle_status']:
            ndata, npoints, dmin, xpos, ypos = \
                self.__match_point(event.xdata, event.ydata)
            if ndata is not None:
                try:
                    for loop in range(len(npoints)):
                        self.data_values[ndata]['plot_mask'][npoints[loop]] = \
                            not \
                            self.data_values[ndata]['plot_mask'][npoints[loop]]
                    self.make_plot()
                except ValueError:
                    pass

    def __match_point(self, xposition, yposition):
        """
        Given a valid plot position find the closest point on the plot.

        Parameters
        ----------
            xposition:   an x plot position in data units (float)

            yposition:   a y plot position in data units (float)

        Returns
        -------
            nset:     an integer value for the data set number of the
                      closest point, or None if there is a problem

            npoints:  a list of integer values for the point numbers within
                      the data set for the closest point or points; there
                      may be more than one as the Vizier data tends to have
                      duplicate entries with the same values.  May be None
                      if there is an issue with the point matching

            dmin:     the minimum "distance" from the input position to the
                      nearest point (due to the log-log plotting this is not
                      the regular distance value...); a float value, or None
                      if there is an issue

            xvalue:   a float value, the x coordinate of the data point, or
                      None if there is a problem in the matching

            yvalue:   a float value, the y coordinate of the data point, or
                      None if there is a problem in the matching

        """
        if (xposition is None) or (yposition is None):
            return None, None, None, None, None
        if len(self.data_values) == 0:
            return None, None, None, None, None
        ytype = self.tkinter_values['variables'][1].get()
        xtype = self.tkinter_values['variables'][2].get()
        xprop = ['wavelength', 'frequency']
        yprop = ['lfl', 'flambda', 'fnu', 'l4fl']
        ndata = None
        for loop in range(len(self.data_values)):
            x_data_values = self.data_values[loop][xprop[xtype]]
            y_data_values = self.data_values[loop][yprop[ytype]]
            if x_data_values is not None:
                if loop == 0.:
                    dmin = math.sqrt((x_data_values[0]/xposition - 1.)**2 +
                                     (y_data_values[0]/yposition - 1.)**2)
                    dmin = 10.*dmin
                distance = numpy.sqrt((x_data_values/xposition-1.)**2 +
                                      (y_data_values/yposition-1.)**2)
                inds = numpy.argsort(distance)
                if distance[inds[0]] < dmin:
                    ndata = loop
                    npoints = []
                    dmin = distance[inds[0]]
                    for newloop in range(len(inds)):
                        if distance[inds[newloop]]/dmin < 1.0001:
                            npoints.append(inds[0])
                    xvalue = x_data_values[npoints[0]]
                    yvalue = y_data_values[npoints[0]]
        if ndata is None:
            return None, None, None, None, None
        return ndata, npoints, dmin, xvalue, yvalue

    def make_extinction_window(self, flagvalue=True):
        """
        This routine makes the extinction control window.  This is used to
        set extinction parameters.  One can then get the paraemters from the
        object for use in the main program.  Closing the window does not
        delete the extinction parameters from the object.

        The plot area and the control area are created here with Tkinter calls.
        Various object variables are defined in the process.

        Parameters
        ----------
            flagvalue:  An optional parameter for whether to redden or
                        deredden values.  The default is True for deredden.
                        Which value should be set depends on whether one is
                        dereddening data values or reddening a model.  The
                        routine can be called from a couple of places so the
                        initial setting depends on the source of the call.

        Returns
        -------
            Nothing is returned from this routine.

        """
        if (len(self.data_values) == 1) and (
                self.data_values[0]['wavelength'] is None):
            return
        if self.extinction_values['rl_wavelengths'] is None:
            wl1 = numpy.asarray([1., 2.])
            extinction, wl1, ext1 = extinction_code.get_extinction(wl1, 1.0,
                                                                   3.1)
            self.extinction_values['rl_wavelengths'] = numpy.copy(wl1)
            self.extinction_values['rl_extinction'] = numpy.copy(ext1)
        extinction_window = Tk.Toplevel(self.root)
        extinction_window.title('Extinction Parameters')
        main_frame = Tk.Frame(extinction_window)
        main_frame.pack(side=Tk.TOP)
        main_frame.config(bg=BGCOL)
        label1 = Tk.Label(main_frame, text='Extinction Function:')
        label1.grid(row=0, column=0)
        label1.config(bg=BGCOL)
        field1 = tkinter.ttk.Combobox(main_frame, width=40)
        field1.grid(row=0, column=1)
        field1['values'] = ["Rieke and Lebofsky (1985)",
                            "Cardelli, Clayton, and Mathis (1989)",
                            "O'Donnell (1994)", "Calzetti (2000)",
                            "Fitzpatrick (1999)",
                            "Fitzpatrick and Massa (2007)"]
        field1.current(0)
        label2 = Tk.Label(main_frame, text='A(V) Value:')
        label2.grid(row=1, column=0)
        label2.config(bg=BGCOL)
        field2 = Tk.Entry(main_frame, width=20)
        field2.grid(row=1, column=1)
        field2.insert(0, '0.0')
        ext_value = Tk.DoubleVar()
        self.field1 = [ext_value, field2]
        ext_slider = Tk.Scale(main_frame, orient=Tk.HORIZONTAL,
                              length=301, from_=0.0, to=20.0,
                              resolution=0.01,
                              variable=ext_value,
                              label="Exinction (mag.)",
                              command=self.apply_av)
        ext_slider.set(0.0)
        ext_slider.grid(row=2, column=0, columnspan=2)
        ext_slider.config(bg=BGCOL)
        label3 = Tk.Label(main_frame, text='R(V) Value:')
        label3.grid(row=3, column=0)
        label3.config(bg=BGCOL)
        field3 = Tk.Entry(main_frame, width=20)
        field3.grid(row=3, column=1)
        field3.insert(0, '3.1')
        r_value = Tk.DoubleVar()
        self.field2 = [r_value, field3]
        r_slider = Tk.Scale(main_frame, orient=Tk.HORIZONTAL,
                            length=301, from_=2.0, to=6.0,
                            resolution=0.01,
                            variable=r_value,
                            label="R(V) Parameter",
                            command=self.apply_rv)
        r_slider.set(3.1)
        r_slider.grid(row=4, column=0, columnspan=2)
        r_slider.config(bg=BGCOL)
        option_frame = Tk.Frame(main_frame)
        option_frame.grid(row=5, column=0, columnspan=2)
        redden_variable = Tk.IntVar()
        tkinter_utilities.put_yes_no(option_frame, redden_variable,
                                     ['De-redden Values', 'Redden Values'],
                                     flagvalue)
        label4 = Tk.Label(main_frame, text='Target Data Set:')
        label4.grid(row=6, column=0)
        label4.config(bg=BGCOL)
        field4 = Tk.Entry(main_frame, width=4)
        field4.grid(row=6, column=1)
        first_data = -1
        first_model = -1
        for loop in range(len(self.data_values)):
            if 'stellar' in self.data_values[loop]['source'][0:7]:
                if first_model < 0:
                    first_model = loop
            else:
                if first_data < 0:
                    first_data = loop
        if flagvalue:
            field4.insert(0, str(first_data+1))
        else:
            field4.insert(0, str(first_model+1))
        button_frame = Tk.Frame(main_frame)
        button_frame.grid(row=7, column=0, columnspan=2)
        close_button = Tk.Button(button_frame, text='Close Window',
                                 command=extinction_window.destroy)
        close_button.pack(side=Tk.TOP)
        close_button.config(bg=BGCOL)
        self.extinction_values['controls'] = {'function': field1,
                                              'av_field': ext_value,
                                              'ebmv_field': r_value,
                                              'flag_variable': redden_variable,
                                              'data_set_field': field4}

    def apply_av(self, event):
        """
        Code to apply the set extinction from a Tkinter slider.

        Parameters
        ----------
        event:   A Tkinter event variable, gives the slider value.

        Returns
        -------
        No value is returned.

        """
        variable = self.field1[0]
        field = self.field1[1]
        self.apply_reddening(variable, field)

    def apply_rv(self, event):
        """
        Code to apply the set extinction from a Tkinter slider.

        Parameters
        ----------
        event:   A Tkinter event variable, gives the slider value.

        Returns
        -------
        No value is returned.

        """
        variable = self.field2[0]
        field = self.field2[1]
        self.apply_reddening(variable, field)

    def apply_reddening(self, variable, field):
        """
        This routine reads the extinction fields and applies the values.

        Parameters
        ----------
        variable:    A float value, the extinction parameter.

        field:       A TKinter Entry field variable, displays the value.

        Returns
        -------
        No values are returned.

        """
        tkinter_utilities.set_field_with_variable(variable, field)
        try:
            set_number = int(
                self.extinction_values['controls']['data_set_field'].get())
            set_number = set_number - 1
            if set_number > len(self.data_values):
                return
            else:
                av = float(
                    self.extinction_values['controls']['av_field'].get())
                if av == 0.:
                    return
                ebmv = float(
                    self.extinction_values['controls']['ebmv_field'].get())
                function = self.extinction_values['controls'][
                    'function'].current()
                flag = self.extinction_values['controls'][
                    'flag_variable'].get()
                flag1 = (flag == 1)
                if function == 0:
                    extinction = extinction_code.get_extinction(
                        self.data_values[set_number]['wavelength'],
                        av, ebmv, function,
                        self.extinction_values['rl_wavelengths'],
                        self.extinction_values['rl_extinction'])
                else:
                    extinction = extinction_code.get_extinction(
                        self.data_values[set_number]['wavelength'],
                        av, ebmv, function)
                newflambda = extinction_code.apply_extinction(
                    self.data_values[set_number]['flambda'],
                    extinction, flag1)
                newfnu = extinction_code.apply_extinction(
                    self.data_values[set_number]['fnu'],
                    extinction, flag1)
                newlfl = extinction_code.apply_extinction(
                    self.data_values[set_number]['lfl'],
                    extinction, flag1)
                newl4fl = extinction_code.apply_extinction(
                    self.data_values[set_number]['l4fl'],
                    extinction, flag1)
                newdflambda = extinction_code.apply_extinction(
                    self.data_values[set_number]['dflambda'],
                    extinction, flag1)
                newdfnu = extinction_code.apply_extinction(
                    self.data_values[set_number]['dfnu'],
                    extinction, flag1)
                newdlfl = extinction_code.apply_extinction(
                    self.data_values[set_number]['dlfl'],
                    extinction, flag1)
                newdl4fl = extinction_code.apply_extinction(
                    self.data_values[set_number]['dl4fl'],
                    extinction, flag1)
                if self.tkinter_values['variables'][3].get() == 0:
                    if self.data_values[set_number]['filter_name'] is not None:
                        for loop in range(
                                len(self.data_values[set_number][
                                    'filter_name'])):
                            if self.flags['template_spectrum'] < 0:
                                value = extinction_code.filter_extinction(
                                    self.data_values[set_number][
                                        'filter_name'][loop], av, ebmv,
                                    function)
                            else:
                                wl1 = numpy.copy(
                                    self.data_values[
                                        self.flags['template_spectrum']][
                                            'wavelength'])
                                sp1 = numpy.copy(
                                    self.data_values[
                                        self.flags['template_spectrum']][
                                            'flambda'])
                                value = extinction_code.filter_extinction(
                                    self.data_values[set_number][
                                        'filter_name'][loop], av, ebmv,
                                    function, wl1, sp1)
                            if value > 1.:
                                newflambda[loop] = value * \
                                    self.data_values[set_number][
                                        'flambda'][loop]
                                newdflambda[loop] = value * \
                                    self.data_values[set_number][
                                        'dflambda'][loop]
                                newfnu[loop] = value * \
                                    self.data_values[set_number]['fnu'][loop]
                                newdfnu[loop] = value * \
                                    self.data_values[set_number]['dfnu'][loop]
                                newlfl[loop] = value * \
                                    self.data_values[set_number]['lfl'][loop]
                                newdlfl[loop] = value * \
                                    self.data_values[set_number]['dlfl'][loop]
                                newl4fl[loop] = value * \
                                    self.data_values[set_number]['l4fl'][loop]
                                newdl4fl[loop] = value * \
                                    self.data_values[set_number]['d4lfl'][loop]
                if 'stellar' in self.data_values[set_number]['source']:
                    label = 'model %d with extinction %.2f/%.2f/%d' % \
                            (set_number+1, av, ebmv, function+1)
                    shortlabel = 'model %d with extinction' % (set_number + 1)
                else:
                    label = 'set %d with extinction %.2f/%.2f/%d' % \
                            (set_number+1, av, ebmv, function+1)
                    shortlabel = 'set %d with extinction' % (set_number + 1)
                for loop in range(len(self.data_values)):
                    if shortlabel in self.data_values[loop]['source']:
                        self.data_values[loop]['flambda'] = numpy.copy(
                            newflambda)
                        self.data_values[loop]['fnu'] = numpy.copy(newfnu)
                        self.data_values[loop]['lfl'] = numpy.copy(newlfl)
                        self.data_values[loop]['l4fl'] = numpy.copy(newl4fl)
                        self.data_values[loop]['dflambda'] = numpy.copy(
                            newdflambda)
                        self.data_values[loop]['dfnu'] = numpy.copy(newdfnu)
                        self.data_values[loop]['dlfl'] = numpy.copy(newdlfl)
                        self.data_values[loop]['dl4fl'] = numpy.copy(newdl4fl)
                        self.data_values[loop]['source'] = label
                        self.make_plot()
                        return
                oldcolour = self.data_values[loop]['plot']['colour']
                for loop in range(len(COLOUR_SET)):
                    if oldcolour == COLOUR_SET[loop]:
                        n1 = loop + 1
                if n1 == len(COLOUR_SET):
                    n1 = 0
                newcolour = COLOUR_SET[n1]
                plot = {
                    'symbol': self.data_values[set_number]['plot']['symbol'],
                    'size': self.data_values[set_number]['plot']['size'],
                    'line': self.data_values[set_number]['plot']['line'],
                    'line_type': self.data_values[set_number][
                        'plot']['line_type'],
                    'line_width': self.data_values[set_number][
                        'plot']['line_width'],
                    'colour': newcolour,
                    'show': True}
                data_set = {
                    'wavelength': numpy.copy(
                        self.data_values[set_number]['wavelength']),
                    'frequency': numpy.copy(
                        self.data_values[set_number]['frequency']),
                    'fnu': newfnu, 'flambda': newflambda, 'lfl': newlfl,
                    'l4fl': newl4fl, 'dl4fl': newdl4fl,
                    'dfnu': newdfnu, 'dflambda': newdflambda, 'dlfl': newdlfl,
                    'mask': numpy.copy(self.data_values[set_number]['mask']),
                    'plot_mask': numpy.copy(
                        self.data_values[set_number]['plot_mask']),
                    'filter_name': self.data_values[set_number]['plot_mask'],
                    'distance': numpy.copy(self.data_values[set_number][
                        'distance']),
                    'position': self.data_values[set_number]['position'],
                    'refpos': self.data_values[set_number]['refpos'],
                    'references': self.data_values[set_number]['references'],
                    'plot': plot,
                    'source': label,
                    'colour_by_name': False}
                self.data_values.append(data_set)
                self.make_plot()
        except:
            pass

    def set_standard_spectrum(self):
        """
        This code allows the user to read in a stellar model as the standard
        spcetral template for extinction calculations.

        Parameters
        ----------
        None.

        Returns
        -------
        None.

        """
        set_numbers = []
        for loop in range(len(self.data_values)):
            if 'stellar' in self.data_values[loop]['source'][0:7]:
                set_numbers.append(loop+1)
        if len(set_numbers) == 0:
            return
        set_number = None
        outstring = \
            'Select standard spectrum shape (-1 for the default): [-1, '
        for n1 in range(len(set_numbers)):
            if n1 < len(set_numbers) - 1:
                outstring = outstring + '%d' % (set_numbers[n1]) + ', '
            else:
                outstring = outstring + '%d' % (set_numbers[n1]) + ']'
        value = tkinter.simpledialog.askstring('Input', outstring)
        if value is not None:
            try:
                i1 = int(value)
                if i1 in set_numbers:
                    set_number = i1
            except ValueError:
                pass
        if set_number is None:
            return
        if set_number > 0:
            self.flags['template_spectrum'] = set_number - 1

    def __plot_position_offsets(self):
        """
        """
        nsets = []
        for loop in range(len(self.data_values)):
            if 'vizier_sed' in self.data_values[loop]['source']:
                nsets.append(loop+1)
        if len(nsets) < 1:
            return
        root1, plot1 = position_plot.make_position_plot(self.data_values,
                                                        nsets[0]-1)
        root1.mainloop()

    def save_data_set(self):
        nsets = len(self.data_values)
        if nsets == 1:
            if self.data_values[0]['wavelength'] is None:
                tkinter.messagebox.showinfp(
                    'Information',
                    'No data sets are available to save.')
                return
            else:
                set_number = 1
        else:
            str1 = 'Enter'
            set_number = askinteger(
                'query','Enter the set number to save (1 or ):')
        if (set_number < 1) or (set_number > nsets):
            tkinter.messagebox.showinfp(
                'Information',
                'Bad data set number %d entered.' % (set_number))
            return
        outfile = tkinter.filedialog.asksaveasfilename()
        if outfile is not None:
            out1 = open(outfile, 'w')
            ndata = len(self.data_values[set_number-1]['wavelength'])
            print('# %s' % (
                self.data_values[set_number-1]['source']),file=out1)
            print('# %d' % (ndata), file=out1)
            print('# %s' % (
                self.data_values[set_number-1]['colour_by_name']),file=out1)
            if self.data_values[set_number-1]['refpos'] is None:
                print('# None None None', file=out1)
            else:
                print('# %f %f %f' % (
                    self.data_values[set_number-1]['refpos'][0],
                    self.data_values[set_number-1]['refpos'][1],
                    self.data_values[set_number-1]['refpos'][2]), file=out1)
            for loop in range(ndata):
                str1 = '%23.10f | ' % (
                    self.data_values[set_number-1]['wavelength'][loop])
                str1 = str1 + '%20.12e | ' % (
                    self.data_values[set_number-1]['frequency'][loop])
                str1 = str1 + '%23.15e | ' % (
                    self.data_values[set_number-1]['flambda'][loop])
                str1 = str1 + '%20.12e | ' % (
                    self.data_values[set_number-1]['dflambda'][loop])
                str1 = str1 + '%20.12e | ' % (
                    self.data_values[set_number-1]['fnu'][loop])
                str1 = str1 + '%20.12e | ' % (
                    self.data_values[set_number-1]['dfnu'][loop])
                str1 = str1 + '%20.12e | ' % (
                    self.data_values[set_number-1]['lfl'][loop])
                str1 = str1 + '%20.12e | ' % (
                    self.data_values[set_number-1]['dlfl'][loop])
                str1 = str1 + '%20.12e | ' % (
                    self.data_values[set_number-1]['l4fl'][loop])
                str1 = str1 + '%20.12e | ' % (
                    self.data_values[set_number-1]['dl4fl'][loop])
                try:
                    str1 = str1 + '%s | ' % (
                        self.data_values[set_number-1]['filter_name'][loop])
                except:
                    str1 = str1 + ' | '
                try:
                    str1 = str1 + '%10.5f | ' % (
                        self.data_values[set_number-1]['distance'][0, loop])
                except:
                    str1 = str1 + ' | '
                try:
                    str1 = str1 + '%10.5f | ' % (
                        self.data_values[set_number-1]['distance'][1, loop])
                except:
                    str1 = str1 + ' | '
                try:
                    str1 = str1 + '%10.5f | ' % (
                        self.data_values[set_number-1]['distance'][2, loop])
                except:
                    str1 = str1 + ' | '
                try:
                    str1 = str1 + '%10.5f | ' % (
                        self.data_values[set_number-1]['distance'][3, loop])
                except:
                    str1 = str1 + ' | '
                try:
                    str1 = str1 + '%20.12f | ' % (
                        self.data_values[set_number-1]['position'][loop][0])
                except:
                    str1 = str1 + ' | '
                try:
                    str1 = str1 + '%20.12f | ' % (
                        self.data_values[set_number-1]['position'][loop][1])
                except:
                    str1 = str1 + ' | '
                try:
                    str1 = str1 + '%s | ' % (
                        self.data_values[set_number-1]['references'][loop])
                except:
                    str1 = str1 + ' | '
                str1 = str1 + '%s | ' % (
                    self.data_values[set_number-1]['mask'][loop])
                str1 = str1 + '%s ' % (
                    self.data_values[set_number-1]['plot_mask'][loop])
                print(str1, file=out1)
            out1.close()
        else:
            return

    def read_data_set(self):
        data_set = {'wavelength': None, 'frequency': None, 'fnu': None,
                    'flambda': None, 'lfl': None, 'l4fl': None,
                    'dfnu': None, 'dflambda': None, 'dl4fl': None,
                    'dlfl': None, 'mask': None, 'plot_mask': None,
                    'filter_name': None, 'distance': None,
                    'position': None, 'refpos': None,
                    'references': None, 'plot': None, 'source': None,
                    'colour_by_name': False}
        filename = tkinter.filedialog.askopenfilename()
        try:
            infile = open(filename, 'r')
            lines = infile.readlines()
            nlines = len(lines)
            infile.close()
            source = lines[0].strip('\n')
            source = source.strip('# ')
            data_set['source'] = source
            frag = lines[1].strip('\n')
            frag = frag.strip('# ')
            ndata = int(frag)
            if nlines != ndata+4:
                raise ValueError
            if 'True' in lines[2]:
                data_set['colour_by_name'] = True
            else:
                data_set['colour_by_name'] = False
            data_set['refpos'] = numpy.zeros((3), dtype=numpy.float32)
            frag = lines[3].strip('\n')
            frag = frag.strip('# ')
            values = frag.split()
            if len(values) == 3:
                data_set['refpos'][0] = float(values[0])
                data_set['refpos'][1] = float(values[1])
                data_set['refpos'][2] = float(values[2])
            data_set['wavelength'] = numpy.zeros((ndata),dtype=numpy.float32)
            data_set['frequency'] = numpy.zeros((ndata),dtype=numpy.float32)
            data_set['flambda'] = numpy.zeros((ndata),dtype=numpy.float32)
            data_set['dflambda'] = numpy.zeros((ndata),dtype=numpy.float32)
            data_set['fnu'] = numpy.zeros((ndata),dtype=numpy.float32)
            data_set['dfnu'] = numpy.zeros((ndata),dtype=numpy.float32)
            data_set['lfl'] = numpy.zeros((ndata),dtype=numpy.float32)
            data_set['dlfl'] = numpy.zeros((ndata),dtype=numpy.float32)
            data_set['l4fl'] = numpy.zeros((ndata),dtype=numpy.float32)
            data_set['dl4fl'] = numpy.zeros((ndata),dtype=numpy.float32)
            data_set['distance'] = numpy.zeros((4, ndata),dtype=numpy.float32)
            data_set['position'] = []
            data_set['mask'] = []
            data_set['filter_name'] = []
            data_set['references'] = []
            data_set['plot_mask'] = []
            distance = numpy.zeros((4), dtype=numpy.float32)
            for loop in range(ndata):
                data_set['position'].append([0., 0.])
                data_set['mask'].append(False)
                data_set['plot_mask'].append(True)
                data_set['references'].append(' ')
                data_set['filter_name'].append(' ')
            for loop in range(ndata):
                line = lines[4+loop].strip('\n')
                values = line.split('|')
                data_set['wavelength'][loop] = float(values[0])
                data_set['frequency'][loop] = float(values[1])
                data_set['flambda'][loop] = float(values[2])
                data_set['dflambda'][loop] = float(values[3])
                data_set['fnu'][loop] = float(values[4])
                data_set['dfnu'][loop] = float(values[5])
                data_set['lfl'][loop] = float(values[6])
                data_set['dlfl'][loop] = float(values[7])
                data_set['l4fl'][loop] = float(values[8])
                data_set['dl4fl'][loop] = float(values[9])
                data_set['filter_name'][loop] = values[10].strip(' ')
                for ind in range(4):
                    try:
                        distance[ind] = float(values[11+ind])
                    except:
                        distance[ind] = 0.
                data_set['distance'][:,loop] = numpy.copy(distance)
                position = [0., 0.]
                try:
                    ra = float(values[15])
                    dec = float(values[16])
                    position = [ra, dec]
                except:
                    pass
                data_set['position'][loop] = position
                data_set['references'][loop] = values[17]
                if values[18].strip(' ') == 'True':
                    data_set['mask'][loop] = True
                else:
                    data_set['mask'][loop] = False
                if values[19].strip(' ') == 'True':
                    data_set['plot_mask'][loop] = True
                else:
                    data_set['plot_mask'][loop] = False
            self.add_data_set(data_set, True)
            self.make_plot()
        except:
            str1 = 'File %s could not be read.' % (filename)
            tkinter.messagebox.showinfo('Error', str1)
            return

        


if __name__ == "__main__":
    """
    This is the main program.

    Assuming that the program is called from the command line one defines
    the window and starts the plot widget.

    """
    root, gui_window = startup()
    root.mainloop()
