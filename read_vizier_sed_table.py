#! /usr/bin/env python
#
"""
This routine allows one to read in a Vizer SED service output .vot file.

Routines defined:

arc_distance     calculate the arc-distance between two sky positions given
                 in degrees (nominally RA and Dec)

get_reference_position   determine reference position for a Vizier table,
                         read from the Target INFO field if possible

read_vizier_vot    the main routine to read data from a Vizier SED .vot file


"""
from astropy.io.votable import parse_single_table
import numpy
from astropy.coordinates import SkyCoord
import astropy.units as units


def arc_distance(ra1, dec1, ra2, dec2):
    """
    Calculate the arc-distance between two positions.

    This is a utility routine.  Given two sets of positions (RA, Dec)
    in degrees the code finds the separation, position angle, and RA/Dec
    offsets of the first position with respect to the second position.

    Parameters
    ----------
        ra1 :   a floating point value, the RA of position 1 in degrees

        dec1 :  a floating point value, the Dec of position 1 in degrees

        ra2 :   a floating point value, the RA of position 2 in degrees

        dec2 :  a floating point value, the Dec of position 2 in degrees

    Position 2 is taken to be reference position, and the distance from this
    point to position 1 is calculated

    Returns
    -------
        separation : a floating point value, the separation in arc-seconds

        angle :      a floating point value, the position angle in
                     degrees E of N

        delra :      a floating point value, the RA separation in arc-seconds

        deldec :     a floating point value, the Dec separation in arc-seconds

    """
    pos1 = SkyCoord(ra1*units.deg, dec1*units.deg)
    pos2 = SkyCoord(ra2*units.deg, dec2*units.deg)
    separation = pos2.separation(pos1)
    angle = pos2.position_angle(pos1)
    pos3 = SkyCoord(ra1*units.deg, dec2*units.deg)
    deldec = pos3.separation(pos1)
    pos4 = SkyCoord(ra2*units.deg, dec1*units.deg)
    ddec = deldec.arcsecond
    if (angle.deg > 90.) and (angle.deg < 270.):
        ddec = -1. * ddec
    delra = pos4.separation(pos1)
    dra = delra.arcsecond
    if angle.deg > 180.:
        dra = -1. * dra
    return separation.arcsecond, angle.deg, dra, ddec


def get_reference_position(filename):
    """
    Read a Vizier SED VOT file and parse the reference position.

    Parmaeters
    ----------
        filename : the name of a Vizier SED VOT file

    Returns
    -------
        ra0 : (float) the reference position RA in decimal degrees

        dec0 : (float) the reference position Dec in decimal degrees

        search_radius : (float) the search radius for the query in
                        arc-seconds

    """
    infile = open(filename)
    lines = infile.readlines()
    infile.close()
    try:
        for line in lines:
            line = line.strip('\n')
            if '<INFO ID="Target"' in line:
                values = line.split('"')
                str1 = values[-2]
                fields = str1.split(',')
                if '+' in fields[0]:
                    v1 = fields[0].split('+')
                    ra0 = float(v1[0])
                    dec0 = float(v1[1])
                else:
                    v1 = fields[0].split('-')
                    ra0 = float(v1[0])
                    dec0 = -1. * float(v1[1])
                search_radius = float(fields[1].replace('rs=', ''))
                return ra0, dec0, search_radius
        return None, None, None
    except:
        return None, None, None


def read_vizier_vot(filename):
    """
    Read a Vizier SED VOT file and parse the reference position.

    Parmaeters
    ----------
        filename : the name of a Vizier SED VOT file

    Returns
    -------
        photom_values :  a floating point numpy array of dimension [14,N]
                         where N is the number of photometry points in the
                         file; the values in the different indexes are

                         0  Wavelength of observation in microns

                         1  Frequency of observation in GHz

                         2  RA of observation in decimal degrees

                         3  Dec of observation in decimal degrees

                         4  Frequency flux density of observation in Jy

                         5  Uncertainty in frequency flux density (Jy)

                         6  Wavelength flux density of observation in
                            W/m^2/micron

                         7  Uncertainty in wavelength flux density
                            (W/m^2/micron)

                         8  wavelength * wavelength flux density of
                            observation (W/m^2)

                         9  uncertainty in wavelength * wavelength flux
                            density (W/m^2)

                        10  distance of observation from reference position
                            (arc-seconds)

                        11  angle of observation from reference position
                            (degrees E of N)

                        12  RA separation (arc-seconds)

                        13  Dec separation (arc-seconds)

        fnu_error_mask :  a numpy boolean array of length N, holding the
                          error mask value, True if the point has no
                          uncertainty value. In such cases the uncertainty
                          is set to zero

        filter_names :  a numpy string array of length N holding the filter
                        names for the observations

        references :  a numpy string array of length N holding the observation
                      reference code from Vizier

        refpos : a three element numpy float vector containing the reference
                 sky position (Ra, Dec) in degrees and the Vizier search
                 radius in arc-seconds; the last value is zero if it cannot
                 be determined from the input file

    If an error occurs the values are returned as None.
    """
    try:
        tab1 = parse_single_table(filename)
    except:
        return None, None, None, None, None
    vizier_fields = ['sed_freq', 'sed_flux', 'sed_eflux', 'sed_filter',
                     '_tabname', '_RAJ2000', '_DEJ2000', '_tabname', '_ID']
    for field in vizier_fields:
        try:
            values = tab1.array[field].data
        except:
            return None, None, None, None, None
    ra = numpy.copy(tab1.array['_RAJ2000'].data)
    dec = numpy.copy(tab1.array['_DEJ2000'].data)
    filters = numpy.copy(tab1.array['sed_filter'].data)
    ramean = numpy.mean(ra)
    decmean = numpy.mean(dec)
    ra0, dec0, search_radius = get_reference_position(filename)
    if ra0 is None:
        search_radius = 0.
        dmin = -1.
        for loop in range(len(filters)):
            if 'GAIA/GAIA2:G' == filters[loop].decode('UTF-8'):
                rain = ra[loop]
                decin = dec[loop]
                distance, angle, delra, deldec = arc_distance(
                    ra[loop], dec[loop], ramean, decmean)
                if (distance < dmin) or (dmin < 0.):
                    dmin = distance
                    ra0 = rain
                    dec0 = decin
    if ra0 is None:
        ra0 = ramean
        dec0 = decmean
        print('Warning:  using the mean RA/Dec as the reference point.')
    npoints = len(tab1.array['sed_freq'].data)
    c = 299792458.
    photom_values = numpy.zeros((14, npoints), dtype=numpy.float32)
    photom_values[1, :] = numpy.copy(tab1.array['sed_freq'].data)
    photom_values[0, :] = numpy.squeeze(
        1.e+06 * c / (photom_values[1, :] * 1.e+09))
    photom_values[4, :] = numpy.copy(tab1.array['sed_flux'].data)
    photom_values[5, :] = numpy.copy(tab1.array['sed_eflux'].data)
    photom_values[8, :] = numpy.squeeze(
        1.e-26 * photom_values[4, :] * photom_values[1, :] * 1.e+09)
    photom_values[6, :] = numpy.squeeze(
        photom_values[8, :] / photom_values[0, :])
    fnu_error_mask = numpy.copy(tab1.array['sed_eflux'].mask)
    photom_values[2, :] = ra
    photom_values[3, :] = dec
    filter_names = []
    references = []
    ids = numpy.copy(tab1.array['_ID'].data)
    for loop in range(npoints):
        filter_names.append(filters[loop].decode('UTF-8'))
        references.append(ids[loop].decode('UTF-8'))
        if fnu_error_mask[loop]:
            photom_values[5, loop] = 0.
            photom_values[7, loop] = 0.
            photom_values[9, loop] = 0.
        else:
            photom_values[9, loop] = (
                photom_values[5, loop] * 1.e-26 *
                (photom_values[1, loop] * 1.e+09))
            photom_values[7, loop] = (
                photom_values[9, loop] / photom_values[0, loop])
        distance, angle, delra, deldec = (
            arc_distance(ra[loop], dec[loop], ra0, dec0))
        photom_values[10, loop] = distance
        photom_values[11, loop] = angle
        photom_values[12, loop] = delra
        photom_values[13, loop] = deldec
    refpos = numpy.asarray([ra0, dec0, search_radius])
    return photom_values, fnu_error_mask, filter_names, references, refpos
