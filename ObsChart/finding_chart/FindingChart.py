#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 09:55:45 2017

@author: Patrick Rauer
"""
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy import units as u
from astroquery.skyview import SkyView
from matplotlib.colors import LogNorm

import numpy as np
import pylab as pl


class FindingChartConfig:
    """
    Class with the finding chart configuration.
    """
    def __init__(self, 
                 width=5.*u.arcmin,
                 height=5.*u.arcmin,
                 survey='DSS2 Red',
                 cross=True,
                 color='r',
                 pixels=120):
        """
        A object of FindingChartConfig contains the global configuration
        of a finding chart

        :param width: The width of the image in deg, arcmin or arcsec
        :type width: astropy.units.quantity.Quantity or float in degrees
        :param height: The height of the image in deg, arcmin or arcsec
        :type height: astropy.units.quantity.Quantity or float in degrees
        :param survey: Name of the survey for the image
        :type survey: str
        :param cross: 
            True if you want to have a cross around the target, else False
        :type cross: bool
        :param color: 
            The color of the name, coordinates, cross, scale and
            coordinate arrow
        :type color: str
        :param pixels: The number of pixels in the image
        :type pixels: int
        """
        if type(width) == float:
            width = width*u.deg
        if type(height) == float:
            height = height*u.deg
        self.width = width
        self.height = height
        self.survey = survey
        self.cross = cross
        self.color = color
        self.pixels = str(int(pixels))


class FindingChart:
    def __init__(self, ra, dec, name='', config=FindingChartConfig()):
        if type(ra) == float or type(ra) == int:
            self.coord = SkyCoord(ra*u.deg, dec*u.deg)
        else:
            self.coord = SkyCoord(ra, dec)
        self.name = name
        self.config = config
        self.color = 'k'
        
    def __image__(self, img, sp):
        """
        Draws the image in a reverse gray scale.
        """
        shape = img.shape
        y = int(shape[0]/2)
        x = int(shape[1]/2)
        img -= np.min(img)
        vmin = np.median(img)
        vmax = np.max(img[y-10: y+10,
                          x-10: x+10])/2
        sp.imshow(img, cmap='gray_r',
                  norm=LogNorm(vmin=vmin,
                               vmax=vmax),
                  origin='lower')
        # name and coordinates of the target
        sp.text(x/20, img.shape[0]-y/5,
                '{}\n{}'.format(self.name, self.coord.to_string('hmsdms')),
                color=self.color, fontsize=13)
        
    def __cross__(self, sp, shape, center):
        """
        Draws the cross around the target.
        """
        x0 = center[0]
        y0 = center[1]
        sp.plot([x0, x0], [y0-y0/20, y0-y0/10], color=self.color)
        sp.plot([x0, x0], [y0+y0/20, y0+y0/10], color=self.color)
        sp.plot([x0-x0/20, x0-x0/10], [y0, y0], color=self.color)
        sp.plot([x0+x0/20, x0+x0/10], [y0, y0], color=self.color)
        
    def __arrows__(self, sp, shape):
        """
        Draws the coordinate arrows in the left bottom
        """
        x = shape[1]
        # axis arrows and labels
        dx = x/10
        dy = shape[0]/10
        sp.arrow(x-dx, dy, 0, 2*dy, color=self.color,
                 head_width=dx/4, linewidth=2)
        sp.arrow(x-dx, dy, -2*dx, 0, color=self.color,
                 head_width=dx/4, linewidth=2)
        sp.text(x-1.25*dx, 3.5*dy, 'N', color=self.color, fontsize=13)
        sp.text(x-4*dx, dy*0.75, 'E', color=self.color, fontsize=13)
        
    def __scale__(self, sp, wcs):
        """
        Draws the scaling line in the right bottom 
        """
        ra, dec = wcs.all_pix2world(20, 20, 0)
        x, y = wcs.all_world2pix(ra-1./60, dec, 0)
        sp.plot([20, x], [20, 20], color=self.color, linewidth=2)
        sp.text(x/2, 30, '1\'', color=self.color, fontsize=13)
        
    def create_chart(self, path=None):
        """
        Creates the finding chart 
        :param path: Path the save place of the finding chart
        :type path: str
        :return:
        """
        paths = SkyView.get_images(position=self.coord,
                                   survey=self.config.survey,
                                   width=self.config.width,
                                   height=self.config.height,
                                   pixels=self.config.pixels)[0][0]
        wcs = WCS(paths.header)
        
        shape = paths.shape
        
        fig = pl.figure(num=0)
        fig.clf()
        sp = fig.add_subplot(111, projection=wcs)
        self.__image__(paths.data.copy(), sp)
        
        # cross around the target
        if self.config.cross:
            center = wcs.all_world2pix(self.coord.ra.degree,
                                       self.coord.dec.degree,
                                       0)
            self.__cross__(sp, shape, center)
            
        self.__arrows__(sp, shape)
        
        self.__scale__(sp, wcs)
        
        if path is not None:
            fig.savefig(path)


ACAM_Config = FindingChartConfig(width=8.3*u.arcmin,
                                 height=8.3*u.arcmin)
if __name__ == '__main__':
    s = SkyCoord('11h40m15.0s -01d58m14.0s')
    ACAM_Config.survey = 'SDSSr'
#    ACAM_Config.width = 2.*u.arcmin
#    ACAM_Config.height = 2.*u.arcmin
    ACAM_Config.pixels = 360
    fc = FindingChart(s.ra.degree, s.dec.degree, 'HCSC candidate',
                      config=ACAM_Config)
    fc.create_chart()
