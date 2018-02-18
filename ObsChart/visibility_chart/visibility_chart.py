#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 17 14:19:24 2017

This script provides the option to plot a visibility chart.
It has a Location class which provides the information of the location,
a Telescope class which provides the information of the telescope and
a VisibilityChart class which does the actual job.
@author: Patrick Rauer
"""

from astropy.time import Time
from PyAstronomy import pyasl
from datetime import datetime, timedelta
import matplotlib.dates as mdates
import pylab as pl
import numpy as np


class Location:
    """
    Location is representer of a position on the earth.
    You need the location to create a visibility chart.
    """
    def __init__(self, latitude, longitude, altitude, east=False,
                 telescopes=None):
        # if the latitude is a string
        if type(latitude) == str:
            latitude = convert_to_deg(latitude)
        # if the longitude is a string
        if type(longitude) == str:
            # if the coordinates are in east direction
            if east:
                longitude = convert_to_deg(longitude)
            else:
                longitude = 360.-convert_to_deg(longitude)
        self.latitude = latitude
        self.longitude = longitude
        self.altitude = altitude
        # if the telescopes input is just one telescope
        if type(telescopes) != list:
            telescopes = [telescopes]
        self.telescopes = telescopes


class Telescope:
    """
    Telescope is a small class to represent a telescope in a visibility chart.
    It contains the name of the telescope and the limits of observabel 
    altitudes.
    """
    def __init__(self, name, min_altitude, max_altitude):
        self.min_altitude = min_altitude
        self.max_altitude = max_altitude
        self.name = name


def __add_twilights__(sp, civil, naut, astro):
    """
    Draws the lines of the twilight.
    Dotted line is civil twilight, dashed lines is nautic twilight and
    full line is astronomical twilight.

    :param sp: Set of axis
    :param civil: List of datetime with civil twilight
    :param naut: List of datetime with nautic twilight
    :param astro: List of datetime with astronomical twilight
    """
    # use the first and the second element
    for i in range(2):
        try:
            sp.plot_date([astro[i], astro[i]], [0, 90], fmt='-k')
            sp.plot_date([naut[i], naut[i]], [0, 90], fmt='--k')
            sp.plot_date([civil[i], civil[i]], [0, 90], fmt=':k')
        # if an index error appears (higher latitudes don't have civil
        # twilight in the summer)
        except IndexError:
            pass


def air_mass(alt):
    """
    Transformation from altitude to air mass.
    :param alt: The altitude
    :type alt: float, numpy.array
    :returns: The air mass of the input values (some order).
    :type: float, numpy.ndarray
    """
    # convert the altitude in a zenith angle in rad
    zenith_angle = np.deg2rad(90-alt)

    r = 12742
    h = 8.4
    # calculate the air mass
    air_mass_value = np.sqrt(np.square(r+h)-np.square(r*np.sin(zenith_angle)))
    air_mass_value -= r*np.cos(zenith_angle)
    air_mass_value /= h
    return air_mass_value


def __jds__(start, end, num):
    """
    Calculates the JDs between the start and the end datetime.

    :param start: The start of the JD list
    :type start: datetime.datetime
    :param end: The end of the JD list
    :type end: datetime.datetime
    :param num: The number of JDs between start and end
    :type num: int
    :returns: A list with the jds between start and end
    :rtype: numpy.ndarray
    """
    start_jd = Time(start, format='datetime').jd
    end_jd = Time(end, format='datetime').jd
    jds = np.linspace(start_jd, end_jd, num=num)
    return jds


def __adjust_dates__(start):
    """
    Redefines start in such a way that the complete night will be covered.
    :param start: The start date of the visibility chart.
    :type start: datetime.datetime
    :returns:
        The adjusted start date and the end date which is 20 hours later
    :rtype: two datetime.datetime items
    """
    start = start-timedelta(hours=start.hour-16, minutes=start.minute,
                            seconds=start.second)
    end = start+timedelta(hours=20)
    return start, end


def __set_axis__(sp, sp2, civil):
    """
    Set the axis labels and transforms the x axis to hours.

    :param sp: first set of axis
    :param sp2: second set of axis
    :param civil: List of datetime's with civil twilight
    """
    # set the limits to the y-axises
    sp.set_ylim(0, 90)
    sp2.set_ylim(5, 1)
    # set the limit to end of the civil twilight in the evening and
    # the start of the civil twilight in the morning
    sp.set_xlim(civil[0], civil[1])
    # convert the x-axis to the format HH:MM
    x_axis_format = mdates.DateFormatter('%H:%M')
    sp.xaxis.set_major_formatter(x_axis_format)
    # rotate the a-axis ticks by 45 deg
    pl.setp(sp.xaxis.get_majorticklabels(), rotation=45)
    # set the labels of the axis
    sp.set_ylabel("Altitude")
    sp.set_xlabel("UTC")
    # add the legend to the chart
    pl.legend(loc='best')
    # convert the second y-axis to the corresponding air mass
    sp2.set_yticklabels(np.round(air_mass(np.linspace(90, 0, 10)), 2))
    sp2.tick_params(axis='y', colors='gray')


def get_altitude(jds, ra, dec, lat, lon, alt):
    # get the altitude, azimuth and hour angle of the sun at the JDs
    alt, az, hour = pyasl.eq2hor(jds, ra, dec,
                                 lat=lat,
                                 lon=lon, 
                                 alt=alt)
    return alt, jds

class VisibilityChart:
    """
    This class creates a complete visibility chart for targets.
    """
    def __init__(self, location=None):
        # if there is no location, use La Palma as the default
        if location is None:
            self.location = LA_PALMA
        else:
            self.location = location

    def __twilight__(self, sun_alt):
        """
        Calculates the time when the sun is closest to the altitude.
        
        :param sun_alt: Searched altitude
        :type sun_alt: float
        :returns: The time when the sun is closest to the altitude
        :rtype: list of datetime
        """
        d_alt = np.abs(self.alt_sun-sun_alt)
        # select all points which are close to the altitude limit 
        d_alt = np.where(d_alt < np.min(d_alt)+0.5)[0]
        # select first element in the range and the last element in the range
        # and convert the JDs to datetime
        times = Time(self.jd_sun[[np.min(d_alt), np.max(d_alt)]],
                     format='jd').datetime
        return times
    
    def twilights(self):
        """
        Calculates the twilight's.
        
        :returns: civil, natuic and astronomical twilight
        """
        # twilight ends when sun is below -6 deg (civil), below -12 (nautic)
        # and below -18 (astronomical)
        civil_times = self.__twilight__(-6)
        naut_times = self.__twilight__(-12)
        astro_times = self.__twilight__(-18)
        return civil_times, naut_times, astro_times
        
    def __sun_pos__(self, jds):
        """
        Calculates the position of the sun for the jds.
        
        :param jds: The time of the sun positions in JD
        :type jds: numpy.ndarray
        """
        # get the sun positions for the jds
        jds, ra_sun, dec_sun = pyasl.sunpos(jds)
        # get the altitude, azimuth and hour angle of the sun at the JDs
        alt, az, hour = pyasl.eq2hor(jds[0], ra_sun[0], dec_sun[0],
                                     lat=self.location.latitude,
                                     lon=self.location.longitude, 
                                     alt=self.location.altitude)
        self.jd_sun = jds[0]
        self.alt_sun = alt

    def __add_telescope_limit__(self, sp, civil):
        """
        Adds the lower limits of the observable altitude of the telescopes
        to the chart.
        
        :param sp: Set of axis
        :param civil: 
            List of datetime's with the start and ending of the observation.
        """
        # if there are telescopes at this location
        if self.location.telescopes is not None:
            # use all telescopes
            for tele in self.location.telescopes:
                # draw the minimal altitude of the telescope
                sp.plot_date([civil[0], civil[1]], 
                             [tele.min_altitude,
                              tele.min_altitude], fmt=':',
                             color='gray')
                # draw the name of the telescope at the edge
                sp.text(civil[1], tele.min_altitude, tele.name)

    def __plot_visibility_curves__(self, sp, targets, jds):
        """
        Plots the visibility curves of the targets.
        
        :param sp: Set of axis
        :param targets: Table or numpy.ndarray with the targets
        :param jds: List with jds of the night (plot dates)
        """
        num = len(jds)
        if type(targets) is dict:
            # calculate the altitude, azimuth and hour angle of the target
            # at the JDs
            alt, az, hour = pyasl.eq2hor(jds,
                                         np.linspace(targets['ra'], targets['ra'], num=num),
                                         np.linspace(targets['dec'], targets['dec'], num=num),
                                         lat=self.location.latitude,
                                         lon=self.location.longitude,
                                         alt=self.location.altitude)
            # convert the JDs back to datetime objects
            times = Time(jds, format='jd').datetime
            # draw the visibility curve
            sp.plot_date(times, alt, label=targets['name'], fmt='-')
        # go through all targets
        for t in targets:
            # calculate the altitude, azimuth and hour angle of the target
            # at the JDs
            alt, az, hour = pyasl.eq2hor(jds,
                                         np.linspace(t['ra'], t['ra'], num=num),
                                         np.linspace(t['dec'], t['dec'], num=num), 
                                         lat=self.location.latitude,
                                         lon=self.location.longitude, 
                                         alt=self.location.altitude)
            # convert the JDs back to datetime objects
            times = Time(jds, format='jd').datetime
            # draw the visibility curve
            sp.plot_date(times, alt, label=t['name'], fmt='-')

    def create_chart(self, targets, start=None, num=500, path=None):
        """
        Creates the visibility chart.
        
        :param targets: 
            Table or numpy.structurearray with the targets. It must have the 
            columns 'name', 'ra' and 'dec'.
        :type targets: numpy.structurearray, astropy.table.Table
        :param start: 
            Start date of the visibility chart. Default is None which means
            today.
        :type start: datetime.datetime
        :param num: Number of points per visibility curve
        :type num: int
        :param path: 
            The path to the save place. Default is None which means no saving
        :type path: str
        """
        # if there is no start date, use the current date
        if start is None:
            start = datetime.utcnow()
        # adjust the start and the end date
        start, end = __adjust_dates__(start)
        # create the JD list
        jds = __jds__(start, end, num)
        
        # calculate the sun position
        self.__sun_pos__(jds)
        # calculate the twilight
        civil, naut, astro = self.twilights()
        
        pl.clf()
        sp = pl.subplot()
        sp2 = sp.twinx()
        # plot the visibility curves
        self.__plot_visibility_curves__(sp, targets, jds)
        # set the axises
        __set_axis__(sp, sp2, civil)
        # add the twilight to the chart
        __add_twilights__(sp, civil, naut, astro)
        # add the telescope limits to the chart
        self.__add_telescope_limit__(sp, civil)
        
        # if there is a path, save the chart at this place
        if path is not None:
            pl.savefig(path)
        # otherwise show the plot
        else:
            pl.show()
        return sp


def convert_to_deg(inp):
    """
    Converts dd:mm:ss.s to degree
    
    :param inp: String in the style dd:mm:ss.s
    :type inp: str
    :returns: degree
    :rtype: float
    """
    c = inp.split(':')
    deg = float(c[0])
    deg += float(c[1])/60
    deg += float(c[2])/3600
    return deg


# Isaac Newton telescope
INT = Telescope('INT', 20, 90)
# Whiliam Harschel Telescope
WHT = Telescope('WHT', 12, 90)
# La Palma
LA_PALMA = Location('28:45:38', '17:52:54', 2344, telescopes=[INT, WHT])
# Tautenburger Exoplanet Survey Telescope, Spain
TEST = Telescope('TEST', 20, 85)
# Alfred Jensch Telescope
AJT = Telescope('AJT', 25, 85)
# Tautenburg, Germany
TAUTENBURG = Location('50:58:48.3', '11:42:41.0', 341,
                      east=True, telescopes=[TEST, AJT])


def get_location_data(location):
    """
    Returns the location objects of the name
    :param location: The name of the location
    :type location: str
    :returns: Location object of these location
    :rtype: Location
    """
    if location == 'La Palma':
        return LA_PALMA
    elif location == 'Tautenburg':
        return TAUTENBURG
    
# direct run (testing)


#if __name__ == '__main__':
#    vc = VisibilityChart(location=TAUTENBURG)
#    from astropy.table import Table
#    
#    tab = Table()
#    tab['ra'] = [20.0, 50.0]
#    tab['dec'] = [25.0, 25.0]
#    tab['name'] = 'test'
#    start_date = datetime.utcnow()
#    end_date = start_date+timedelta(hours=12)
#    sp_chart = vc.create_chart(tab, start_date)
