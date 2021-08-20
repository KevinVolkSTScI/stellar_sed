#! /usr/bin/env python
#
"""
Various utility routine for reading in stellar models, or manipulating them.

Each routine from

read_bosz_model            read a BOSZ model (fits or ascii)
read_phoenix_model         read a France Allard phoenix model
read_phoenix_grid_model    read one of the Husser et al. (2013) phoenix grid
                           models
read_calspec_model         read one of the HST Calspec model fits files

takes as input a file name and returns the model spectrum values

wavelengths, flambda

where the wavelengths are in microns and the flambda values are in
Watts/meter squared/micron.  If there is a problem these routines
return None for both variables.

An additional routine

read_generic_model         read a generic ascii spectral model

takes as input a file name and the columns to read and returns the wavelength
and spectrum information.

Two other routines are included, which are used to smooth the Phoenix models
to lower wavelength resoution:

__convolve_spectrum   convolve an input spectrum with a Gaussian function that
                      has constant relative sigma value (this is done using
                      logarithmic sampling of the wavelengths)

__rebin_spectrum      rebin a convolved spectrum to a lower wavelength
                      resolution

Also two routines related to scaling models to data values are included:

match_data            given a set of input data values with uncertainties,
                      do a least-squares fit to scale a model to the data,
                      or do a mean ratio fit to scale a model to the data
                      as an alternative

__calculate_photon_mean_flux_density    calculate the photon wieghted mean flux
                                        density for a model spectrum, if the
                                        name of the filter is recognised so
                                        that the profile can be read in

"""
import math
import scipy.signal as signal
import numpy
import astropy.io.fits as fits
import os
from astropy.io.votable import parse_single_table
WAVELNORM = 2.2
FLAMBDANORM = 5.e-15


def read_bosz_model(filename, wlnorm=0.0, flnorm=0.):
    """
    Read in a BOSZ model from a file.  Scale to about K = 7.

    Given a file name, this routine attempts to read in a BOSZ model from
    the file.  The file can be a FITS table file or the ascii version.  The
    ascii files have suffix values either '.asc' or '.asc.bz2'.  The FITS
    files have suffix '.fits'.  If the file cannot be read, the return
    values are None.

    The input model wavelength range is 0.1 to 32 microns.  A simply power-law
    extrapolation is used to extend the spectrum to 1000 microns.  As the
    last flux point or two are sometimes lower than expected from an simple
    extrapolation of preceding few points, the tenth and ninth last points
    are used to define the power law index.

    Parameters
    ----------

    filename:   A string variable giving the file name to be read in.

    wlnorm:     An optional float value, the wavelength to normalize to in
                microns.  If no value is given or the wavelength is outside
                the model range a value of WAVELNORM microns is used.

    flnorm:     An optional float value, the wavelength flux density to
                normalize to in W/m^2/micron.  If the value is zero,
                negative, or larger than 1.e-05 a default value of FLAMBDANORM
                is used.

    Returns
    -------

    wavelengths:  A numpy float array of wavelength values in microns

    flambda:      A numpy array of F_lambda values in W/m^2/micron.  The
                  input model spectrum is normalized to the flnorm value at
                  wavelength wlnorm.
    """
    try:
        if '.fits' in filename:
            s1 = fits.open(filename)
            tab1 = s1[1].data
            s1.close()
            wavelengths = numpy.copy(tab1['WAVELENGTH']/10000.)
            flambda = numpy.copy(tab1['SPECIFICINTENSITY']*10.*3.14159265359)
            tab1 = 0.
        elif '.asc.bz2' in filename:
            os.system('bunzip2 '+filename)
            newfilename = filename.replace('.bz2', '')
            wavelengths = numpy.loadtxt(newfilename, usecols=(0,))/10000.
            flambda = numpy.loadtxt(newfilename, usecols=(1, )) * 10. * \
                3.14159265359
            os.system('bzip2 '+newfilename)
        elif '.asc' in filename[-4:]:
            wavelengths = numpy.loadtxt(newfilename, usecols=(0,))/10000.
            flambda = numpy.loadtxt(newfilename, usecols=(1,)) * 10. * \
                3.14159265359
        else:
            return None, None
        if (wlnorm < wavelengths[0]) or (wlnorm > wavelengths[-1]):
            wlnorm = WAVELNORM
        if (flnorm <= 0.) or (flnorm > 1.e-05):
            flnorm = FLAMBDANORM
        wl1 = wavelengths[-10]
        wl2 = wavelengths[-9]
        fl1 = flambda[-10]
        fl2 = flambda[-9]
        try:
            power = math.log10(fl2/fl1)/math.log10(wl2/wl1)
        except ValueError:
            power = -4.0
        if power > -3.:
            power = -4.0
        fl0 = numpy.interp(wlnorm, wavelengths, flambda)
        if fl0 <= 0.:
            return None, None
        scale = flnorm/fl0
        addwavelengths = numpy.arange(32.1, 100.15, 0.2)
        morewavelengths = numpy.arange(101., 1010.5, 10.)
        addwavelengths = numpy.append(addwavelengths, morewavelengths)
        addflambda = fl2*numpy.power(addwavelengths/wl2, power)
        wavelengths = numpy.append(wavelengths, addwavelengths)
        flambda = numpy.append(flambda, addflambda)
        flambda = flambda * scale
        return wavelengths, flambda
    except:
        return None, None


def read_phoenix_model(filename, resolution=1000.0):
    """
    Read in a France Allard  Phoenix model file (ascii).

    Try to read one of the Phoenix models tabuated in the Allard web pages,
    such as the NextGen or BTSettl models.  There are occational oddities
    in these files so success of reading is not certain, but it works for
    many of the files.

    The model is scaled down to be reasonable for a 2MASS detection.

    Parameters
    ----------

    filename:   A string value, the file name to read

    resolution:  An optional floating point value, the target resolution for
                 the model; if the value is between 100 and 3000 it is applied
                 to the Phoenix model.

    Returns
    -------

    wavelengths:   A numpy float array of the wavelengths in microns, or
                   None if there is a problem

    spectrum:      A numpy float array of the F_lambda values in W/m^2/micron,
                   or None if there is a problem

    """
    try:
        newfilename = filename
        if '.gz' in filename:
            os.system('gunzip '+filename)
            newfilename = filename.replace('.gz', '')
        if '.bz2' in filename:
            os.system('bunzip2 '+filename)
            newfilename = filename.replace('.bz2', '')
        if '.xz' in filename:
            os.system('unxz '+filename)
            newfilename = filename.replace('.xz', '')
        if ('.spec' in newfilename[-5:]):
            # This code is to read the older ".spec.gz" Allard models
            infile = open(newfilename, 'r')
            line = infile.readline()
            line = infile.readline()
            line = line.strip('\n')
            values = line.split()
            nwavel = int(values[0])
            nlines = int(float(nwavel)/4.)
            if 4*nlines < nwavel:
                nlines = nlines+1
            wavelengths = numpy.zeros((nwavel), dtype=numpy.float32)
            spectrum = numpy.zeros((nwavel), dtype=numpy.float32)
            n1 = 0
            for loop in range(nlines):
                line = infile.readline()
                line = line.strip('\n')
                values = line.split()
                try:
                    wavelengths[n1] = float(values[0])
                    n1 = n1+1
                    wavelengths[n1] = float(values[1])
                    n1 = n1+1
                    wavelengths[n1] = float(values[2])
                    n1 = n1+1
                    wavelengths[n1] = float(values[3])
                    n1 = n1+1
                except:
                    pass
            nlines = int(float(nwavel)/6.)
            if 6*nlines < nwavel:
                nlines = nlines+1
            n1 = 0
            wavelengths = wavelengths/10000.
            for loop in range(nlines):
                line = infile.readline()
                line = line.strip('\n')
                values = line.split()
                try:
                    spectrum[n1] = float(values[0])
                    n1 = n1+1
                    spectrum[n1] = float(values[1])
                    n1 = n1+1
                    spectrum[n1] = float(values[2])
                    n1 = n1+1
                    spectrum[n1] = float(values[3])
                    n1 = n1+1
                    spectrum[n1] = float(values[4])
                    n1 = n1+1
                    spectrum[n1] = float(values[5])
                    n1 = n1+1
                except:
                    pass
            infile.close()
        else:
            # this is for the ".7.gz" Allard grid models
            infile = open(newfilename, 'r')
            lines = infile.readlines()
            infile.close()
            nlines = len(lines)
            wavelengths = numpy.zeros((nlines), dtype=numpy.float32)
            spectrum = numpy.zeros((nlines), dtype=numpy.float32)
            for n1 in range(nlines):
                values = lines[n1].split()
                wavelengths[n1] = float(values[0])/10000.
                spectrum[n1] = 10.**(float(values[1].replace('D', 'E')))
        if '.bz2' in filename:
            os.system('bzip2 '+newfilename)
        elif '.xz' in filename:
            os.system('xz '+newfilename)
        elif '.gz' in filename:
            os.system('gzip '+newfilename)
        inds = numpy.where(spectrum > 0.)
        wavelengths = wavelengths[inds]
        spectrum = spectrum[inds]
        inds = numpy.argsort(wavelengths)
        wavelengths = wavelengths[inds]
        spectrum = spectrum[inds]
        delwl = wavelengths[1:] - wavelengths[0:-1]
        if wavelengths[-1] < 50.5:
            newwavelengths = numpy.arange(wavelengths[-1]+0.5, 300.1, 0.5)
            newspectrum = spectrum[-1]*numpy.power(
                (wavelengths[-1]/newwavelengths), 4)
            wavelengths = numpy.append(wavelengths, newwavelengths)
            spectrum = numpy.append(spectrum, newspectrum)
            newwavelengths = numpy.arange(301., 1000., 1.)
            newspectrum = spectrum[-1]*numpy.power(
                (wavelengths[-1]/newwavelengths), 4)
            wavelengths = numpy.append(wavelengths, newwavelengths)
            spectrum = numpy.append(spectrum, newspectrum)
        if delwl[-1] > 6.:
            inds = numpy.where(delwl > 6.)
            newwavelengths = numpy.zeros((0),dtype=numpy.float32)
            newspectrum = numpy.zeros((0),dtype=numpy.float32)
            for loop in range(inds[0][0], len(wavelengths)-1):
                wl0 = wavelengths[loop]
                fl0 = spectrum[loop]
                wl1 = wavelengths[loop+1]
                fl1 = spectrum[loop+1]
                power = math.log(fl1/fl0)/math.log(wl1/wl0)
                nwl = numpy.arange(wl0+0.5, wl1-0.05, 0.5)
                nfl = fl0*numpy.power((nwl/wl0), power)
                newwavelengths = numpy.append(newwavelengths, nwl)
                newspectrum = numpy.append(newspectrum, nfl)
            wavelengths = numpy.append(wavelengths, newwavelengths)
            spectrum = numpy.append(spectrum, newspectrum)
            inds = numpy.argsort(wavelengths)
            wavelengths = wavelengths[inds]
            spectrum = spectrum[inds]
        if (resolution >= 100.0) and (resolution <= 3000.0):
            spectrum = __convolve_spectrum(wavelengths, spectrum, resolution)
            wavelengths, spectrum = __rebin_spectrum(
                wavelengths, spectrum, resolution)
        wlnorm = WAVELNORM
        flnorm = FLAMBDANORM
        fl0 = numpy.interp(wlnorm, wavelengths, spectrum)
        if fl0 <= 0.:
            return None, None
        scale = flnorm/fl0
        spectrum = spectrum*scale
        return wavelengths, spectrum
    except:
        return None, None


def read_phoenix_grid_model(filename, resolution=1000.0):
    """
    Read in a Husser et al. (2013) Phoenix grid model file (fits table).

    The F_lambda values are read from the named file.  The wavelengths are
    in a separate fits file named 'WAVE_PHOENIX-ACES-AGSS-COND-2011.fits'
    that is assumed to be in the same directory as the model file.

    The model is scaled down to be reasonable for a 2MASS detection.  Also,
    additional long wavelength points are added assuming a Rayleigh-Jeans
    spectral shape.

    Parameters
    ----------

    filename:    A string value, the file name to read

    resolution:  An optional floating point value, the target resolution for
                 the model; if the value is between 100 and 3000 it is applied
                 to the Phoenix model.

    Returns
    -------

    wavelengths:   A numpy float array of the wavelengths in microns, or
                   None if there is a problem

    spectrum:      A numpy float array of the F_lambda values in W?m^2/micron,
                   or None if there is a problem

    """
    try:
        spectrum = fits.getdata(filename)
        values = filename.split('/')
        newfilename = filename.replace(values[-1],
                                       'WAVE_PHOENIX-ACES-AGSS-COND-2011.fits')
        wavelengths = fits.getdata(newfilename)
        # convert wavelengths to microns
        wavelengths = wavelengths/10000.
        fl0 = spectrum[-1]
        wl0 = wavelengths[-1]
        if (resolution >= 100.0) & (resolution <= 3000.0):
            spectrum = __convolve_spectrum(wavelengths, spectrum, resolution)
            wavelengths, spectrum = __rebin_spectrum(
                wavelengths, spectrum, resolution)
        # add points to longer wavelengths; assume a power law of -4.0
        addwavelengths = numpy.arange(5.6, 1006., 1.)
        addspectrum = fl0*numpy.power(addwavelengths/wl0, -4.0)
        wavelengths = numpy.append(wavelengths, addwavelengths)
        spectrum = numpy.append(spectrum, addspectrum)
        # renormalize the spectrum to a reasonable value
        wl0 = WAVELNORM
        fl0 = FLAMBDANORM
        flnorm = numpy.interp(wl0, wavelengths, spectrum)
        scale = fl0/flnorm
        if scale <= 0.:
            scale = 1.0
        spectrum = spectrum*scale
        return wavelengths, spectrum
    except:
        return None, None


def read_calspec_model(filename, resolution=10000):
    try:
        hdu1 = fits.open(filename)
        table = hdu1[1].data
        wavelengths = numpy.copy(table['wavelength'])/10000.0
        spectrum = numpy.copy(table['flux'])*10.
        if numpy.max(wavelengths) < 32.2:
            wl1 = wavelengths[-10]
            wl2 = wavelengths[-9]
            fl1 = spectrum[-10]
            fl2 = spectrum[-9]
            try:
                power = math.log10(fl2/fl1)/math.log10(wl2/wl1)
            except:
                power = -4.0
            if power < -3.:
                power = -4.0
            addwavelengths = numpy.arange(numpy.max(wavelengths)+0.1,
                                          100.15, 0.2)
            morewavelengths = numpy.arange(101., 1010.5, 10.)
            addwavelengths = numpy.append(addwavelengths, morewavelengths)
            addflambda = fl2*numpy.power(addwavelengths/wl2, power)
            wavelengths = numpy.append(wavelengths, addwavelengths)
            spectrum = numpy.append(spectrum, addflambda)
            if (resolution >= 100.0) & (resolution <= 3000.0):
                spectrum = __convolve_spectrum(wavelengths, spectrum,
                                               resolution)
                wavelengths, spectrum = __rebin_spectrum(
                    wavelengths, spectrum, resolution)
        return wavelengths, spectrum
    except:
        return None, None


def read_generic_model(filename, wlcol, spcol, comment=['#', '\\', '|']):
    """
    Read in a spectrum from a generic ascii or fits table file

    Parameters
    ----------

    filename:   A string variable giving the file name to be read.  If the
                name ends in '.fits' it is assumed to be a FITS table file,
                otherwise it is assumed to be in ascii format.

    wlcol:      Either a string variable for the FITS table column name for
                the wavelengths, or an integer column number for the
                wavelengths to be read from.

    spcol:      Either a string variable for the FITS table column name for
                the flux density values, or an integer column number for the
                flux density values to e read from.

    comment:    An optional list of characters which are taken to mark
                comment lines.  The default is for the characters '#',
                '\', and |' be taken to mark comment lines.

    Returns
    -------

    wavelengths:  A numpy float array of the wavelength (or possibly frequency)
                  values.  The units are not specified.  Any transformations
                  need to be done in the calling program.

    spectrum:     A numpy float array of the flux density values.  The units
                  are not specified.  Any transformations need to be done in
                  the calling program.
    """
    try:
        if '.fits' in filename:
            s1 = fits.open(filename)
            try:
                tab1 = s1[1].data
            except:
                tab1 = s1[0].data
            s1.close()
            wavelengths = tab1[wlcol]
            spectrum = tab1[spcol]
        else:
            wavelengths = numpy.loadtxt(filename, comments=comment,
                                        usecols=[wlcol-1, ])
            spectrum = numpy.loadtxt(filename, comments=comment,
                                     usecols=[spcol-1, ])
        inds = numpy.argsort(wavelengths)
        wavelengths = wavelengths[inds]
        spectrum = spectrum[inds]
        return wavelengths, spectrum
    except:
        return None, None


def __convolve_spectrum(wavelengths, spectrum, resolution):
    """
    Convolve an input spectrum to a fixed resolution.

    This file takes an assumed high resolution spectrum (most likely from a
    stellar model) and convolves it with a Gaussian kernel in logarithmic
    wavelength sampling so that the relative resolution is the same for all
    wavelengths.  The revised spectrum is returned, or if the resolution is
    out of range then the input spectrum is returned with no changes.

    This formulation assumes a large wavelength range for the input
    spectrum (several orders of magnitude) but with wide variations in
    the resolution.  This will not work well for a spectrum with a more
    restricted wavelength range, say a factor of 2 or so.  Hence it is not
    very useful for observed spectra.

    Parameters
    ----------

    wavelengths:   A numpy float array containing the original wavelength
                   wavelength sampling

    spectrum:      A numpy float array of the wavelength or frequency
                   flux density values for the wavelengths

    resolution:    A float or integer number for the effective resolution
                   that is to be applied, in the range from 10 to 10000.
                   This value is assumed to be the ratio of the mean and
                   sigma values of the Gaussian function used in the
                   convolution.

    Returns
    -------

    outspectrum:     A numpy float array with the new spectrum values.
                     The dimension is the same as for the original
                     spectrum.

    """
    if (resolution < 10.) or (resolution > 10000.):
        return spectrum
    wlmin = numpy.min(wavelengths)
    wlmax = numpy.max(wavelengths)
    delwl = 0.000005
    wl0 = 0.5
    a1 = math.log10(wlmin)
    a2 = math.log10(wlmax)
    dela = math.log10(1. + delwl/wl0)
    newwavelengths = numpy.power(10., numpy.arange(a1, a2, dela))
    newspectrum = numpy.interp(newwavelengths, wavelengths, spectrum)
    newspectrum[newspectrum < 0.] = 0.
    sigma = wl0/resolution
    b1 = math.log10(wl0 - 5.*sigma)
    b2 = math.log10(wl0 + 5.*sigma)
    gausswl = numpy.power(10., numpy.arange(b1, b2, dela))
    gauss = 1./numpy.exp((gausswl - wl0)*(gausswl - wl0)/(2.*sigma*sigma))
    sum1 = numpy.sum(gauss)
    gauss = gauss/sum1
    outspectrum = signal.fftconvolve(newspectrum, gauss, mode='same')
    ratio = outspectrum/newspectrum
    ratio[newspectrum == 0.] = 0.
    inds = numpy.where(newwavelengths >= 0.09)
    newwavelengths = newwavelengths[inds[0][0]:]
    outspectrum = numpy.copy(outspectrum[inds[0][0]:])
    outspectrum = numpy.interp(wavelengths, newwavelengths, outspectrum,
                               left=0., right=newspectrum[-1])
    return outspectrum


def __rebin_spectrum(wavelengths, spectrum, resolution):
    """
    Rebin an input spectrum at constant wavelength resolution.

    This file takes an assumed high resolution spectrum (most likely from a
    stellar model) convolved with a gaussian profile and reduces the
    wavelength sampling.

    Parameters
    ----------

    wavelengths:   A numpy float array containing the original wavelength
                   wavelength sampling

    spectrum:      A numpy float array of the wavelength or frequency
                   flux density values for the wavelengths

    resolution:    A float or integer number for the effective resolution
                   that is to be applied, in the range from 10 to 10000.

    Returns
    -------

    outwavelengths:  A numpy float array with the new wavelength values

    newspectrum:     A numpy float array with the new spectrum values

    """
    if resolution < 10. or resolution > 10000.:
        return wavelengths, spectrum
    wlmin = numpy.min(wavelengths)
    wlmax = numpy.max(wavelengths)
    a1 = math.log10(wlmin)
    a2 = math.log10(wlmax)
    dela = math.log10(1.+1./resolution)
    newwavelengths = numpy.power(10., numpy.arange(a1, a2+dela/2., dela))
    outwavelengths = (newwavelengths[1:] + newwavelengths[0:-1])/2.
    newspectrum = outwavelengths * 0.
    for loop in range(len(newwavelengths)-1):
        wl1 = newwavelengths[loop]
        wl2 = newwavelengths[loop+1]
        delwl = wl2 - wl1
        subwl = numpy.arange(wl1, wl2+delwl/1000., delwl/100.)
        subsp = numpy.interp(subwl, wavelengths, spectrum)
        newspectrum[loop] = numpy.trapz(subsp, subwl)/(wl2 - wl1)
    return outwavelengths, newspectrum


def match_data(wavelengths, values, errors, labels, weights,
               modelwl, modelvalues, errorflag=True,
               photon_mean_flag=False, ratio_flag=False):
    """
    Given an input set of data values, scale the model to match.

    This uses a straightforward least-squares fitting allowing for
    weighting.

    Parameters
    ----------

    wavelengths:   A numpy float array of data wavelengths (can be a mixture
                   of photometry and spectral values).

    values:        A numpy float array of flux density values at the
                   wavelengths.

    errors:        A numpy float array of uncertainties in the data values.
                   The inverse square of the uncertainties will be used in
                   the fitting if errorflag is True.  Error values of 0.0
                   are assigned a relative error value of 5% where needed in
                   the calculations.

    labels:        A list of strings containing the filter names (or None if
                   the input does not list the names).  Used if the photon
                   mean flux density option is selected.

    weights:       A numpy float array of relative weights to use in the
                   fitting, if errorflag is False.

    modelwl:       A numpy float array of the model wavelengths.

    modelvalues:   A numpy float array of the model flux density values.

    errorflag:     A boolean variable that determines whether the errors or
                   the weights are used in the fitting.

    photon_mean_flag:   A boolean variable that determines whether the
                        photon wieghted mean flux density values are
                        calculated for the fit.  This requires filter
                        profiles, so only a sub-set of filters can use this
                        option.

    ratio_flag:    A boolean variable that determines whether the ratio is 
                   used in the fitting or the arithmatic deviation is used 
                   in the fitting.  The default is to use the ratio, so 
                   the relative deviations are minimized.

    Returns
    -------

    scale:     The least-squares best fit scaling factor to match the model
               to the data points.

    stats:     A numpy float array giving goodness of fit information.

    """
    newvalues = numpy.interp(wavelengths, modelwl, modelvalues)
    if photon_mean_flag:
        for loop in range(len(labels)):
            mean_flux_density = __calculate_photon_mean_flux_density(
                modelwl, modelvalues, wavelengths[loop], labels[loop])
            if mean_flux_density > 0.:
                newvalues[loop] = mean_flux_density
    if errorflag:
        ratio = values/errors
        ratio[errors == 0.] = 20.
        newweights = ratio*ratio
    else:
        newweights = weights*1.
    newerrors = numpy.copy(errors)
    newerrors[errors == 0.] = 0.05*values[errors == 0.]
    if ratio_flag:
        scale = numpy.sum(values*newweights)/numpy.sum(newvalues*newweights)
    else:
        sigratio = values/newvalues
        sigratio[newvalues == 0.] = 1.
        scale = numpy.sum(sigratio*newweights)/numpy.sum(newweights)
    stats = numpy.zeros((4), dtype=numpy.float32)
    ratio = newvalues/values
    inds = numpy.where(values > 0.)
    ratio = ratio[inds]
    newvalues = newvalues[inds]
    # the following is the root mean square fractional deviation
    rms = numpy.sqrt(numpy.mean((1. - ratio)*(1. - ratio)))
    stats[0] = rms
    # calculate the RMS deviation in terms of quoted uncertainties
    ratio1 = ((values - newvalues)/newerrors)
    sigmarms = numpy.sqrt(numpy.mean(ratio1*ratio1))
    stats[1] = sigmarms
    # stats will have more values later
    return scale, stats


def __calculate_photon_mean_flux_density(wavelengths, spectrum, wl0, label):
    """
    Calculate the photon weighted mean flux density for a filter.

    Given as input a model spectrum, the filter wavelength, and the
    filter name, attempt to read in the filter profile from the SIP files
    and calculate the photon weighted mean flux density.  If there is an
    issue, the model is interpolated to the wavelength and this value is
    returned.

    Note that the filter profiles are assumed to be in a sub-directory named
    "filter_subset" under either the current directory or the directory that
    is defined by the $EXINCTION_PATH environment variable.  The latter takes
    precidence if the variable is defined.

    Parameters
    ----------

    wavelengths:    A numpy float array of the wavelength values in microns for
                    the model spectrum

    spectrum:   A numpy float arrat of the wavelength flux density values
                for the model spectrum (nominally in units of W/m^2/micron)

    wl0:        A float value, the filter wavelength in microns, used for
                interpolation if the filter profile cannot be read in

    label:      A string variable holding the Vizier filter name, or None

    Returns
    -------

    newvalue:   A floating point value, eEither the photon weighted mean
                flux density, or the interpolated flux density at wl0 if
                the photon weighted mean flux density calculation fails

    """
    newvalue = numpy.interp(wl0, wavelengths, spectrum)
    profile_names = ['filter_2MASS_2MASS_H.vot', 'filter_2MASS_2MASS_J.vot',
                     'filter_2MASS_2MASS_Ks.vot', 'filter_GAIA_GAIA2r_G.vot',
                     'filter_GAIA_GAIA2r_Gbp.vot',
                     'filter_GAIA_GAIA2r_Grp.vot',
                     'filter_Generic_Bessell_B.vot',
                     'filter_Generic_Bessell_I.vot',
                     'filter_Generic_Bessell_R.vot',
                     'filter_Generic_Bessell_U.vot',
                     'filter_Generic_Bessell_V.vot',
                     'filter_Generic_Cousins_I.vot',
                     'filter_Generic_Cousins_R.vot',
                     'filter_Generic_Johnson_B.vot',
                     'filter_Generic_Johnson_I.vot',
                     'filter_Generic_Johnson_R.vot',
                     'filter_Generic_Johnson_U.vot',
                     'filter_Generic_Johnson_V.vot',
                     'filter_Generic_Stromgren_b.vot',
                     'filter_Generic_Stromgren_u.vot',
                     'filter_Generic_Stromgren_v.vot',
                     'filter_Generic_Stromgren_y.vot',
                     'filter_SLOAN_SDSS_g.vot',
                     'filter_SLOAN_SDSS_gprime_filter.vot',
                     'filter_SLOAN_SDSS_i.vot',
                     'filter_SLOAN_SDSS_iprime_filter.vot',
                     'filter_SLOAN_SDSS_r.vot',
                     'filter_SLOAN_SDSS_rprime_filter.vot',
                     'filter_SLOAN_SDSS_u.vot',
                     'filter_SLOAN_SDSS_uprime_filter.vot',
                     'filter_SLOAN_SDSS_z.vot',
                     'filter_SLOAN_SDSS_zprime_filter.vot',
                     'filter_UKIRT_UKIDSS_H.vot', 'filter_UKIRT_UKIDSS_J.vot',
                     'filter_UKIRT_UKIDSS_K.vot', 'filter_UKIRT_UKIDSS_Y.vot',
                     'filter_UKIRT_UKIDSS_Z.vot', 'filter_WISE_WISE_W1.vot',
                     'filter_WISE_WISE_W2.vot', 'filter_WISE_WISE_W3.vot',
                     'filter_WISE_WISE_W4.vot']
    label_names = ["2MASS:H", "2MASS:J", "2MASS:Ks", "Gaia:G", "Gaia:Gbp",
                   "Gaia:Grp", "Bessell:B",  "Bessell:I",  "Bessell:R",
                   "Bessell:U", "Bessell:V", "Cousins:I", "Cousins:R",
                   "Johnson:B", "Johnson:I", "Johnson:R", "Johnson:U",
                   "Johnson:V", "Stromgren:b", "Stromgren:u", "Stromgren:v",
                   "Stromgren:y", "SDSS:g", "SDSS:g'", "SDSS:i", "SDSS:i'",
                   "SDSS:r", "SDSS:r", "SDSS:u", "SDSS:u'", "SDSS:z",
                   "SDSS:z'", "UKIDSS:H", "UKIDSS:J", "UKIDSS:K", "UKIDSS:Y",
                   "UKIDSS:Z", "WISE:W1", "WISE:W2", "WISE:W3", "WISE:W4"]
    try:
        path = os.environ['EXTINCTION_PATH']
        if path[-1] != '/':
            path = path + '/'
    except:
        path = './'
    fpath = path + 'filter_subset/'
    if label not in label_names:
        return newvalue
    for loop in range(len(label_names)):
        if label == label_names[loop]:
            tab1 = parse_single_table(fpath + profile_names[loop])
            filterwl = tab1.array['Wavelength'].data/10000.
            response = numpy.copy(tab1.array['Transmission'].data)
            midwl = numpy.mean(filterwl)
            i1 = numpy.where(wavelengths > midwl)
            delwl1 = wavelengths[i1[0][0]] - wavelengths[i1[0][0]-1]
            i2 = numpy.where(filterwl > midwl)
            delwl2 = filterwl[i2[0][0]] - filterwl[i2[0][0]-1]
            if delwl1 < delwl2:
                newwavelengths = numpy.copy(wavelengths)
                newspectrum = numpy.copy(spectrum)
                newresponse = numpy.interp(wavelengths, filterwl,
                                           response, left=0., right=0.)
            else:
                newresponse = response
                newspectrum = numpy.interp(filterwl, wavelengths, spectrum)
                newwavelengths = filterwl
            sum1 = numpy.trapz(newwavelengths*newresponse*newspectrum)
            sum2 = numpy.trapz(newwavelengths*newresponse)
            if (sum2 > 0.) and (sum1 > 0.):
                return sum1/sum2
            return newvalue
    return newvalue
