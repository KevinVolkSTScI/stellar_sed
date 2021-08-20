#! /usr/bin/env python
#
"""
Various utility codes to define/apply extinction corrections.

The current routines are:

get_extinction:   A routine to return extinctions values in magnitudes for
                  a set of wavelengths and the two extinction parameters
                  A_V and R_V.  Several possible extinction functions can
                  be used in this routine.  (In some cases the R_V parameter
                  is not used.)

apply_extinction:  A routine to apply extinction values from get_extinction
                   to a set of flux density values; can either redden or
                   deredden the values depending upon a flag.

filter_extinction:  A routine to calculate filter extinction for one of a
                    sub-set of widely used filters (e.g. Jonhnson UBV,
                    SDSS ugirz, 2MASS JHKs, WISE W1 to W4).  This assumes the
                    Sirius spectral shape and photon detection in the
                    calculations.  The actual extinction will vary to some
                    degree with the spectral type of the source, so this is an
                    approximation to the actual extinction correction.
"""
import os
import numpy
import extinction
from astropy.io.votable import parse_single_table
import astropy.io.fits as fits


def get_extinction(wavelengths, avpar, rvpar, flag=0,
                   rl_wavelengths=None, rl_extinction=None):
    """
    General routine to apply one of several extinction functions.

    Parameters
    ----------

    wavelengths:   A numpy float array of wavelength values in microns at
                   which the extinction is to be calculated, assumed to be
                   one-dimensional

    avpar:         A real number >= 0., the A_V parameter for the extiction
                   (magnitudes)

    rvpar:         A real number, the ratio of A_V to E(B-V) parameter for
                   the extinction (denoted R_V); this is not used for some
                   of the extinction laws

    flag:          An integer flag for the extinction law (optional)

                   0  -> Rieke and Lebofsky (1985)
                   1  -> extinction package, Cardelli, Clayton & Mathis (1989)
                         extinction function
                   2  -> extinction package, O’Donnell (1994) extinction
                         function
                   3  -> Calzetti et al. (2000) extinction function
                   4  -> extinction package, Fitzpatrick (1999) dust
                         extinction function
                   5  -> extinction package, Fitzpatrick & Massa (2007)
                         extinction model for R_V = 3.1.

    rl_wavelengths:   The Rieke and Lebofsky (1985) wavelength array.  The
                       array is read in if the value is None.

    rl_extinction:   The Rieke and Lebofsky (1985) normalized extinction
                     array.  The array is read in as if the value is None.

    Returns
    -------

    extinction_values:   A numpy float array of the same dimension as
                         wavelengths, giving the extinction values in
                         magnitudes

    There are two optional return values.

    rl_wavelengths:   If the rl_wavelengths value is None in the call
                      and flag = 0, the array is returned.  It is a numpy
                      float array.

    rl_extinction:   If the rl_extinction value is None in the call
                     and flag = 0, the array is returned.  It is a numpy
                     float array.

    If the extinction values are negative, they are set to zero.  It is not
    clear whether this can ever happen in the extinction package.  It will
    not happen for the Rieke and Lebofsky (1985) extinction curve.
    """
    if avpar <= 0.:
        if avpar < 0.:
            print('Error, the A_V value must be positive')
        return wavelengths * 0.
    if len(wavelengths.shape) > 1:
        print('Error, the wavelength array must be one-dimensional')
        return wavelengths * 0.
    return_arrays = False
    if (flag == 0) and (
            (rl_wavelengths is None) or (rl_extinction is None)):
        return_arrays = True
        try:
            path = os.environ['EXTINCTION_PATH']
            if path[-1] != '/':
                path = path + '/'
        except ValueError:
            path = './'
        try:
            rl_wavelengths = numpy.loadtxt(path+'extinction.values',
                                           usecols=[0, ])
            rl_extinction = numpy.loadtxt(path+'extinction.values',
                                          usecols=[1, ])
#            print('Have read in the Reike/Lebofky extinction arrays.')
        except:
            print('Error: could not read the Rieke/Lebofsky extinction values')
            print('Will use the Cardelli, Clayton, and Mathis (1989) values.')
            newwavelengths = wavelengths * 10000.
            extinction_values = extinction.ccm89(newwavelengths, avpar,
                                                 rvpar)
            return extinction_values
    newwavelengths = wavelengths * 10000.
    if flag == 0:
        values = rl_extinction * avpar
        extinction_values = numpy.interp(wavelengths, rl_wavelengths,
                                         values)
        if return_arrays:
            return extinction_values, rl_wavelengths, rl_extinction
        return extinction_values
    if (flag == 1) or ((flag == 0) and (rl_wavelengths[0] == 0)):
        extinction_values = extinction.ccm89(newwavelengths, avpar, rvpar)
    if flag == 2:
        extinction_values = extinction.odonnell94(newwavelengths, avpar, rvpar)
    if flag == 3:
        extinction_values = extinction.calzetti00(newwavelengths, avpar, rvpar)
    if flag == 4:
        extinction_values = extinction.fitzpatrick99(newwavelengths, avpar,
                                                     rvpar)
    if flag == 5:
        extinction_values = extinction.fm07(newwavelengths, avpar)
    extinction_values[extinction_values < 0.] = 0.
    return extinction_values


def apply_extinction(fluxes, extinction_values, deredden=True):
    """
    Apply extinction values to an array of flux density values.

    Parameters:

    fluxes:   A numpy float array of flux density values, assumed
              to be 1-dimensional

    extinction_values:  A numpy float array of extinction values in
                        magnitudes corresponding to the flux density
                        values

    deredden:   An optional boolean flag; if True the extinction is used
                to deredden the input fluxes (e.g. deredden observed
                photometry) or if False the extinction is used to redden
                the input fluxes (e.g. redden a model spectrum).  The
                default is to deredden the input flux densities

    Returns
    -------

    newfluxes:   A numpy float array of the corrected flux density values.
                 Since the extinction is a dimensionless number, the units
                 of the output array are the same as those of the input
                 fluxes array.

    """
    correction = numpy.power(10., 0.4*extinction_values)
    if deredden:
        newfluxes = fluxes*correction
    else:
        newfluxes = fluxes/correction
    return newfluxes


def filter_extinction(filter_name, avpar, ebmvpar, flag=0, wavelengths=None,
                      spectrum=None):
    """
    General routine to apply one of several extinction functions to a filter.

    Parameters
    ----------

    filter_name:   A string giving the Vizier VOT name of the filter

    avpar:         A real number >= 0., the A_V parameter for the extiction
                   (magnitudes)

    rvpar:         A real number, the ratio of A_V to E(B-V) parameter for
                   the extinction (denoted R_V); this is not used for some
                   of the extinction laws

    flag:          An optional integer flag for the extinction law

                   0  -> Rieke and Lebofsky (1985)
                   1  -> extinction package, Cardelli, Clayton & Mathis (1989)
                         extinction function
                   2  -> extinction package, O’Donnell (1994) extinction
                         function
                   3  -> Calzetti (2000) extinction function
                   4  -> extinction package, Fitzpatrick (1999) dust
                         extinction function
                   5  -> extinction pacakge, Fitzpatrick & Massa (2007)
                         extinction model for R_V = 3.1.

    wavelengths:   An optional numpy float array of wavelengths for the
                   standard spectral shape; if it is None then the Sirius
                   spectral template is used

    spectrum:      An optional numpy float array of wavelength flux density
                   values for the standard spectral shape; if it is None then
                   the Sirius spectral template is used

    Returns
    -------

    factor:    The relative change in photon flux for the extinnction
               parameters, a float value >= 1.0.

    """
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
                   "Gaia:Grp", "Bessell:B", "Bessell:I", "Bessell:R",
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
    except ValueError:
        path = './'
    fpath = path + 'filter_subset/'
    if filter_name not in label_names:
        return 1.0
    if (wavelengths is None) or (spectrum is None):
        try:
            hdu1 = fits.open(path + 'sirius_mod_004.fits')
            table = hdu1[1].data
            wavelengths = numpy.copy(table['wavelength'])/10000.0
            spectrum = numpy.copy(table['flux'])*10.
            addwl = numpy.arange(300., 2010., 10.)
            wl1 = wavelengths[-3]
            fl1 = spectrum[-3]
            addfl = fl1 * numpy.power(wl1/addwl, 4.)
            spectrum = numpy.append(spectrum, addfl)
            wavelengths = numpy.append(wavelengths, addwl)
        except:
            print('Error getting the Sirius spectrum')
            return 1.0
    for loop in range(len(label_names)):
        if filter_name == label_names[loop]:
            tab1 = parse_single_table(fpath + profile_names[loop])
            filterwl = tab1.array['Wavelength'].data/10000.
            response = numpy.copy(tab1.array['Transmission'].data)
            midwl = numpy.mean(filterwl)
            minds1 = numpy.where(wavelengths > midwl)
            delwl1 = wavelengths[minds1[0][0]] - wavelengths[minds1[0][0]-1]
            minds2 = numpy.where(filterwl > midwl)
            delwl2 = filterwl[minds2[0][0]] - filterwl[minds2[0][0]-1]
            if delwl1 < delwl2:
                newwavelengths = numpy.copy(wavelengths)
                newspectrum = numpy.copy(spectrum)
                newresponse = numpy.interp(wavelengths, filterwl,
                                           response, left=0., right=0.)
            else:
                newresponse = response
                newspectrum = numpy.interp(filterwl, wavelengths, spectrum)
                newwavelengths = filterwl
            flux0 = numpy.trapz(newspectrum*newwavelengths*newresponse,
                                newwavelengths)
            if flag > 0:
                ext1 = get_extinction(newwavelengths, avpar, ebmvpar, flag)
            else:
                ext1, rlwavelm, rlvalues = get_extinction(
                    newwavelengths, avpar, ebmvpar, flag)
            extspectrum = apply_extinction(newspectrum, ext1, False)
            flux1 = numpy.trapz(extspectrum*newwavelengths*newresponse,
                                newwavelengths)
            if (flux1 <= 0.) or (flux0 <= 0.):
                return 1.0
            factor = flux0/flux1
            return factor
