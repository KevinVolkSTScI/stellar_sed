#! /usr/bin/env python
#
"""
This file contains a number of utilities for dealing with spectral energy
distributions.

trans_wavelenth    transform wavelength/frequency units

trans_flux_density   transformflux density values to different units

sort_positions     sort posiion values

integrate_sed      integrate over the SED to get a total flux value

vizier_means       find the mean wavelength and flux density values for
                   unique Vizier filter names

"""
import numpy


def trans_wavelength(data_values, option, target):
    """
    A utility code to transform wavelength/frequency values.

    Parameters
    ----------

    data_values:   A numpy float array of wavelength or frequency values

    option:        An integer flag for the units of the input data values

    target:        An integer flag for the units of the output data values

    Returns
    -------

    new_data_values:  A numpy flat array of the same dimensions as the input
                      data_values array with the new wavelength or frequency
                      values

    The options and target flags can take the values

        0   for wavelength values in microns

        1   for wavelength values in Angstroms

        2   for frequency values in GHz

        3   for frequency values in Hz

    """
    c = 299792458.
    if option == target:
        return data_values
    if option == 0:
        if target == 1:
            return data_values * 10000.
        frequency = c / (data_values * 1.e-06)
        if target == 2:
            return frequency / 1.e+09
        return frequency
    if option == 1:
        if target == 0:
            return data_values / 10000.
        frequency = c / (data_values * 1.e-10)
        if target == 2:
            return frequency / 1.e+09
        return frequency
    if option == 2:
        if target == 3:
            return data_values * 1.e+09
        wavelengths = c / (data_values * 1.e+09)
        if target == 0:
            return wavelengths * 1.e-06
        return wavelengths * 1.e-10
    if option == 3:
        if target == 2:
            return data_values * 1.e-09
        wavelengths = c / data_values
        if target == 0:
            return wavelengths * 1.e-06
        return wavelengths * 1.e-10
    return data_values


def trans_flux_density(wavelengths, data_values, option, target):
    """
    A utility code to transform flux density values.

    Parameters
    ----------

    wavelengths:   A numpy float array of wavelength values in microns

    data_values:   A numpy float array of flux density values

    option:        An integer flag for the units of the data_values array

    target:        An integer flag for the output flux density units

    Returns
    -------

    new_data_values:   A numpy array of the same dimensions as data_values
                       giving the transformed flux density values

    The option/target values can be the following:

        0:   denotes wavelength flux density in Watt/meter^2/micron units

        1:   denotes wavelength flux density in
             erg/second/centimeterm^2/Angstrom units

        2:   denotes wavelength * wavelength flux density in
             Watt/meter^2 units

        3:   denotes wavelength * wavelength flux density in
             erg/second/centimeter^2 units

        4:   denotes frequency flux density in Watt/meter^2/Hz

        5:   denotes frequency flux density in Jansky (1.e-26 Watt/meter^2/Hz)

    """
    if option == target:
        return data_values
    c = 299792458.
    frequency = c / (wavelengths * 1.e-06)
    flambda = numpy.copy(data_values)
    if option == 1:
        flambda = data_values * 10.0
    elif option == 2:
        flambda = flambda / wavelengths
    elif option == 3:
        flambda = flambda / 1000. / wavelengths
    elif option == 4:
        flambda = data_values * frequency / wavelengths
    elif option == 5:
        flambda = data_values * 1.e-26 * frequency / wavelengths
    if target == 0:
        return flambda
    if target == 1:
        return flambda / 10.
    if target == 2:
        return flambda * wavelengths
    if target == 3:
        return flambda * wavelengths * 1000.
    if target == 4:
        return flambda * wavelengths / frequency
    if target == 5:
        return 1.e+26 * flambda * wavelengths / frequency
    return data_values


def sort_positions(positions):
    """
    Utility routine to take two corner values are return the x and y range.

    Parameters
    ----------

    positions:   A four element list or numpy array of float values; values
                 are assumed to be [x1. y1, x2, y2] for two corner positions
                 (x1, y1) and (x2, y2).

    Returns
    -------

    xmin:    A floating point value, the x minimum defining a rectangular area
             or None if there is an issue

    xmax:    A floating point value, the x maximum defining a rectangular area
             or None if there is an issue

    ymin:    A floating point value, the y minimum defining a rectangular area
             or None if there is an issue

    ymax:    A floating point value, the y maximum defining a rectangular area
             or None if there is an issue

    """
    if len(positions) != 4:
        return None, None, None, None
    for loop in range(4):
        if positions[loop] is None:
            return None, None, None, None
    xmin = min(positions[0], positions[2])
    xmax = max(positions[0], positions[2])
    ymin = min(positions[1], positions[3])
    ymax = max(positions[1], positions[3])
    return xmin, xmax, ymin, ymax


def integrate_sed(wavelength, flambda, wlmin=None, wlmax=None):
    """
    Calculate the flux in an SED by direct integration.

    A direct trapezoidal rule integration is carried out on the flambda values
    and the associated wavelength values.

    Parameters
    ----------

        wavelength:  A numpy float array of wavelength values, normally in
                     microns

        flambda:     A numpy float array of flux density values, normally
                     F_lambda in W/m^2/micron

        wlmin:       An optional float value for the minimum wavelength of
                     the calculation, or None to have no lower limit aside
                     from the data range

        wlmax:       An optional float value for the maximum wavelength of
                     the calculation, or None to have no upper limit aside
                     from the data range

    Returns
    -------

        flux1:   The float value, the estimated total flux, nominally in
                 W/m^2 if the input units are microns and W/m^2/micron; if
                 the wavelength range is bad or the two arrays do not match
                 in length a value of zero is returned

    """
    if len(wavelength) != len(flambda):
        return 0.
    if wlmin is None:
        xmin = 0.9 * numpy.min(wavelength)
    else:
        xmin = wlmin
    if wlmax is None:
        xmax = 1.1 * numpy.max(wavelength)
    else:
        xmax = wlmax
    if (xmin >= xmax) or (len(wavelength) < 2):
        return 0.
    inds = numpy.argsort(wavelength)
    newwavelength = numpy.copy(wavelength[inds])
    newflambda = numpy.copy(flambda[inds])
    if (xmin > numpy.min(wavelength)) or (xmax < numpy.max(wavelength)):
        fl1 = numpy.interp(xmin, wavelength, flambda)
        fl2 = numpy.interp(xmax, wavelength, flambda)
        newwavelength[newwavelength < xmin] = xmin
        newwavelength[newwavelength > xmax] = xmax
        newflambda[newwavelength < xmin] = fl1
        newflambda[newwavelength > xmax] = fl2
    flux = numpy.trapz(newflambda, newwavelength)
    return flux


def vizier_means(wavelengths, flambda, filter_names):
    """

    Parameters
    ----------
    wavelengths:  a numpy float array, of wavelength values

    flambda:      a numpy float array of wavelength flux density values

    filter_names:  a numpy string array, holding the filter names

    Returns
    -------
    wlout:     a numpy float array of output wavelength values
    flout:     a numpy float array of output wavelength flux density values
    outnames:  a numpy string array of output unique filter names

    """
    outnames = numpy.unique(filter_names)
    wlout = numpy.zeros((len(outnames)), dtype=numpy.float32)
    flout = numpy.zeros((len(outnames)), dtype=numpy.float32)
    for loop in range(len(outnames)):
        inds = numpy.where(filter_names == outnames[loop])
        subwl = wavelengths[inds]
        subfl = flambda[inds]
        wlout[loop] = subwl[0]
        flout1 = numpy.unique(subfl)
        flout[loop] = numpy.mean(flout1)
    if numpy.min(flout) <= 0.:
        inds = numpy.where(flout > 0.)
        wlout = wlout[inds]
        flout = flout[inds]
        outnames = outnames[inds]
    return wlout, flout, outnames
