# stellar_sed
Code to plot photometry and spectra for astronomical objects and compare with stellar models (or other types of models if the model is saved in an ascii file of the proper format).  The code brings up a tkinter/matplotlib interface within which one can read in photometry and spectral data for astronomical sources and plot them, along with options to read in several types of stellar model spectra for comparison.  The code is specifically intended to use the Vizier photometry information that can be found for astronomical sources at the SIMBAD Astronomical Database (http://simbad.u-strasbg.fr/simbad).  If one queries the VizieR SED service (http://vizier.u-strasbg.fr/vizier/sed/) and saves the output as a VOT table the code is able to read the information in and plot the values.  Indeed, the SED plots that this site produces are the model for the current code, with added functionality to also plot spectra and stellar models on the same graph.

The code requires the normal Python packages sys, os, math, matplotlib, and tkitner.  It also requires numpy, scipy, astropy, and extinction.  The last few packages can all be installed via "pip install numpy", "pip install scipy", "pip install astropy", and "pip install extinction".  The "extinction" package comes from https://github.com/kbarbary/extinction and the astropy package is found at https://www.astropy.org/.

The following files constitute the code:

extinction_code.py

model_utilities.py

position_plot.py

read_vizier_sed_table.py

sed_plot_interface.py

sed_utilities.py

tkinter_utilities.py

In addition the code looks for the file extinction.values which should be in the same directory as the code.  Further, a set of filter profiles are needed for the code to compare stellar models to observed photometry.  These are stored in the filter_subset directory.  If one runs the code from a directory other than where the code files are found, one needs to define the $EXTINCTION_PATH environment variable to point to the directory where the code and the filter_subset directory are found.

The code is run by starting the "sed_plot_interface.py" code from the command line.  Hence this file needs to be executable.  

The code requires python version 3.  It has been tested on a MacOS system and on a Linux/Unix system.  It has not been tested on a Microsoft Windows system.

Description of how to use the code is found in the documentation files sed_display_and_fitting_code.docx and sed_display_and_fitting_code.pdf.

