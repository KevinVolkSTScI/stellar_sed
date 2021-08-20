#! /usr/bin/env python
#
"""
This file contains a number of general utility routines for Tkinter.

put_yes_no     Make a two-section radio button

set_field_with_variable    Set a text field with a value from a variable

set_entry      Insert a string into an Entry field

separator_line   Create a line within a frame for separating items

append_text      Append a text string to Tkinter text box or Scrolled Text area

list_sets_for_entry    Make a string to interest into a label field of a
                       messagebox entry widget

"""
import tkinter as Tk
# define a background colour for windows
BGCOL = '#F8F8FF'


def put_yes_no(root, var, labels, flag):
    """
    Create a Tkineter yes/no radio button.

    This is a utility routine to make a yes/no radio button in Tkinter.
    The required variables are passed to the code and it produces the
    button pair.

    Parameters
    ----------
        root :  The Tk frame variable to hold the buttons

        var :   The Tk IntVar that is used to communicate with the
                buttons

        labels : A two element list with strings for the two states of
                 the buttons, first "yes" then "no".

        flag :  A boolean value that causes the code to set the yes field
                (the first of the two) if it is True.  If the value is
                False the no field (the second of the two) is set instead.

    Returns
    -------
        No values are returned by this routine.

    """
    yesfield = Tk.Radiobutton(root, text=labels[0], variable=var, value=1)
    yesfield.grid(row=0, column=0, sticky=Tk.W)
    yesfield.config(bg=BGCOL)
    nofield = Tk.Radiobutton(root, text=labels[1], variable=var, value=0)
    nofield.grid(row=0, column=1, sticky=Tk.W)
    nofield.config(bg=BGCOL)
    if flag:
        yesfield.select()
    else:
        nofield.select()


def set_field_with_variable(variable, field):
    """
    Utility routine: set an entry field with a Tkinter variable value

    Parameters
    ----------

    variable:   A TKinter variable (e.g. from Tk.IntVar or Tk.DoubleVar) to
                insert into the field.

    field:      A Tkinter Entry field variable (from Tk.Entry).  The current
                content is removed and the variable value is put in the field.

    Returns
    -------

    No values are returned.

    """
    value = variable.get()
    nchar = len(field.get())
    field.delete(0, last=nchar)
    field.insert(0, str(value))


def set_entry(entry, value):
    """
    Utility routine to set a Tkinter entry field to a string.

    Parameters
    ----------
        entry:  A Tkinter Entry variable

        value:  A string variable holding the new contents to put in the
                entry field

    Returns
    -------
        No values are returned from this routine.

    """
    try:
        old_string = entry.get()
        entry.delete(0, last=len(old_string))
        entry.insert(0, value)
    except:
        pass


def separator_line(parent, dims, flag, packvalue=Tk.TOP, bgcol='#d9d9d9'):
    """
    Create the Tkinter canvas object making a separator line.

    This is a utility routine to make a separator line within a Tkinter
    frame.

    Parameters
    ----------
        parent :  A Tkinter Frame variable, that will contain the line

        dims :    A three-element list with the line dimensions

            dims[0] :  An integer value for the line canvas width (pixels)

            dims[1] :  An integer value for the line canvas height (pixels)

            dims[2] :  An integer value for the line padding (pixels)

        flag :    A Boolean value, the line is horizontal if the
                  value is True, vertical otherwise

        packvalue :  An optional Tkinter pack direction (one of Tk.LEFT,
                     Tk.RIGHT, Tk.TOP, or Tk.BOTTOM) with default value
                     of Tk.TOP

        bgcol :  An optional colour value for the window background.

    Returns
    -------
        linecanvas:  the Tkinter Canvas variable for the line, in case
                     the user wishes to modify it later

    For a vertical line normally the height will be much larger than the
    width, as in

    sepline = separator_line(frame, [10, 300, 5], False)

    while for a horizontal line normally the width will be much larger
    than the height as in

    sepline = separator_line(frame, [500, 5, 5], True)

    """
    lincanvas = Tk.Canvas(parent, height=dims[1], width=dims[0])
    lincanvas.config(bg=bgcol)
    try:
        lincanvas.pack(side=packvalue, fill=Tk.BOTH, expand=Tk.YES)
    except:
        return None
    if flag:
        lincanvas.create_line(dims[2], dims[1]/2, dims[0] - dims[2],
                              dims[1]/2)
    else:
        lincanvas.create_line(dims[0]/2, dims[2], dims[0]/2,
                              dims[1] - dims[2])
    return lincanvas


def append_text(text_area, outstring):
    """
    Utility routine to add a string to the end of a text box.

    Parameters
    ----------

    text_area:  A tkinter text box or scrolled text variable

    outstring:  A string value to be appended to the text area

    Returns
    -------

    Nothing
    """
    text_area.insert(Tk.END, outstring)
    text_area.see(Tk.END)


def list_sets_for_entry(nsets):
    """

    Parameters
    ----------
    nsets:    List of integers, sets to select from

    Returns
    -------
    str1: sring variable, the list of sets (for a simple dialog)

    """
    str1 = '%d' % nsets[0]
    if len(nsets) > 1:
        for loop in range(1, len(nsets)):
            str1 = str1 + ', %d' % (nsets[loop])
    return str1
