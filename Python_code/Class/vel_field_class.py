from __future__ import print_function
# class for carrying out analysis of 2D velocity fields created via the 
# cube class. i.e. the input data file should be a 2D field of velocity 
# measurements, which are made via the pipeline class initially 

# import the relevant modules
import os, sys, numpy as np, random, math
import pyraf
import numpy.polynomial.polynomial as poly
import lmfit
import scipy
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.mlab as mlab
from scipy import stats
from scipy import interpolate
from lmfit.models import GaussianModel, PolynomialModel
from scipy import optimize
from lmfit import Model
from scipy.optimize import minimize
from scipy.optimize import basinhopping
from astropy.io import fits
from pylab import *
from matplotlib.colors import LogNorm
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MaxNLocator
from scipy.spatial import distance
from copy import copy

####################################################################

class vel_field(object):

    """
    Def: 
    Class for housing analysis methods relevant to 2D velocity 
    field data computed elsewhere
    Input: 
    2D array of velocity field data (IMAGE)
    """

    # Initialiser creates an instance of the cube object
    # Input must be a combined data cube with two extensions - data and noise 

    def __init__(self, fileName, c_rot_x, c_rot_y):

        """
        Def:
        Initialiser method 
        Input: filename - name of file containing 1D spectrum 
                z - The redshift of the galaxy being manipulated
        """
        self.self = self
        # Initialise the fileName object 
        self.fileName = fileName
        # initialise the centre of rotation in the x and y
        self.c_rot_x = c_rot_x
        self.c_rot_y = c_rot_y
        # Variable containing all the fits extensions 
        self.Table = fits.open(fileName)
        # create an object to house the data 
        self.data = self.Table[0].data
        # set the lengths of the bins as x & y
        xbin = np.arange(0, self.data.shape[0], 1)
        ybin = np.arange(0, self.data.shape[1], 1)
        self.xbin, self.ybin = np.meshgrid(xbin, ybin)
        self.xbin = np.ravel(self.xbin)
        self.ybin = np.ravel(self.ybin)

        # initialise the flattened velocity array
        self.vel_flat = []

        # now that we have the xbins and ybins, these are the coordinates 
        # at which we want to evaluate the velocity in the 
        # fit_kinematic_pa method. We can loop around the 
        # the coordinates and create a flattened array 
        # so that the velocity data values correspond to the bins 

        for x, y in zip(self.xbin, self.ybin):
            # evaluate the velocity point 
            self.vel_flat.append(self.data[x][y])

        # make sure that vel_flat is a numpy array 
        self.vel_flat = np.array(self.vel_flat)

        # now can subtract the central positions from 
        # xbins and ybins 
        self.xbin = self.xbin - self.c_rot_x
        self.ybin = self.ybin - self.c_rot_y


        # print (self.vel_flat)

        # print (len(self.vel_flat))
        # print (len(self.xbin))
        # print (len(self.ybin))
        # print (len(self.xbin) * len(self.ybin))
        # sauron colour dictionary
        self._cdict = {'red':((0.000,   0.01,   0.01),
                             (0.170,   0.0,    0.0),
                             (0.336,   0.4,    0.4),
                             (0.414,   0.5,    0.5),
                             (0.463,   0.3,    0.3),
                             (0.502,   0.0,    0.0),
                             (0.541,   0.7,    0.7),
                             (0.590,   1.0,    1.0),
                             (0.668,   1.0,    1.0),
                             (0.834,   1.0,    1.0),
                             (1.000,   0.9,    0.9)),
                    'green':((0.000,   0.01,   0.01), 
                             (0.170,   0.0,    0.0),
                             (0.336,   0.85,   0.85),
                             (0.414,   1.0,    1.0),
                             (0.463,   1.0,    1.0),
                             (0.502,   0.9,    0.9),
                             (0.541,   1.0,    1.0),
                             (0.590,   1.0,    1.0),
                             (0.668,   0.85,   0.85),
                             (0.834,   0.0,    0.0),
                             (1.000,   0.9,    0.9)),
                     'blue':((0.000,   0.01,   0.01),
                             (0.170,   1.0,    1.0),
                             (0.336,   1.0,    1.0),
                             (0.414,   1.0,    1.0),
                             (0.463,   0.7,    0.7),
                             (0.502,   0.0,    0.0),
                             (0.541,   0.0,    0.0),
                             (0.590,   0.0,    0.0),
                             (0.668,   0.0,    0.0),
                             (0.834,   0.0,    0.0),
                             (1.000,   0.9,    0.9))
                      }

        self.sauron = colors.LinearSegmentedColormap('sauron', self._cdict)


    def _rotate_points(self, x, y, ang):
        """
        Rotates points counter-clockwise by an angle ANG in degrees.
        Michele cappellari, Paranal, 10 November 2013
        
        """
        theta = np.radians(ang - 90.)
        xNew = x*np.cos(theta) - y*np.sin(theta)
        yNew = x*np.sin(theta) + y*np.cos(theta)
        return xNew, yNew


    def symmetrize_velfield(self, xbin, ybin, velBin, sym=2, pa=90.):
        """
        This routine generates a bi-symmetric ('axisymmetric') 
        version of a given set of kinematical measurements.
        PA: is the angle in degrees, measured counter-clockwise,
          from the vertical axis (Y axis) to the galaxy major axis.
        SYM: by-simmetry: is 1 for (V,h3,h5) and 2 for (sigma,h4,h6)

        """

        xbin, ybin, velBin = map(np.asarray, [xbin, ybin, velBin])
        x, y = self._rotate_points(xbin, ybin, -pa)  # Negative PA for counter-clockwise
        
        xyIn = np.column_stack([x, y])
        xout = np.hstack([x,-x, x,-x])
        yout = np.hstack([y, y,-y,-y])
        xyOut = np.column_stack([xout, yout])
        velOut = interpolate.griddata(xyIn, velBin, xyOut)
        velOut = velOut.reshape(4, xbin.size)
        if sym == 1:
            velOut[[1, 3], :] *= -1.
        velSym = np.nanmean(velOut, axis=0)
        
        # print ('This is the symmetrised vel_field: %s' % np.nanmean(velSym))
        return velSym

    def plot_velfield(self, x, y, vel, vmin=None, vmax=None, ncolors=64, nodots=False,
                      colorbar=False, label=None, flux=None, fixpdf=False,
                      nticks=7, **kwargs):

        if vmin is None:
            vmin = -1000

        if vmax is None:
            vmax = 1000

        x, y, vel = map(np.ravel, [x, y, vel])
        levels = np.linspace(vmin, vmax, ncolors)

        ax = plt.gca()
        cs = ax.tricontourf(x, y, vel.clip(vmin, vmax), levels=levels,
                           cmap=kwargs.get("cmap", self.sauron))

        ax.axis('image')  # Equal axes and no rescaling
        ax.minorticks_on()
        ax.tick_params(length=10, which='major')
        ax.tick_params(length=5, which='minor')

        if flux is not None:
            ax.tricontour(x, y, -2.5*np.log10(flux/np.max(flux).ravel()),
                          levels=np.arange(20), colors='k') # 1 mag contours

        if fixpdf:  # remove white contour lines in PDF at expense of larger file size
            ax.tricontour(x, y, vel.clip(vmin, vmax), levels=levels, zorder=0,
                          cmap=kwargs.get("cmap", self.sauron))

        if not nodots:
            ax.plot(x, y, '.k', markersize=kwargs.get("markersize", 3))

        if colorbar:
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.1)
            ticks = MaxNLocator(nticks).tick_values(vmin, vmax)
            cbar = plt.colorbar(cs, cax=cax, ticks=ticks)
            if label:
                cbar.set_label(label)


        return cs

    def display_pixels(self, x, y, val, pixelsize=None, angle=None, **kwargs):
        """
        Display vectors of square pixels at coordinates (x,y) coloured with "val".
        An optional rotation around the origin can be applied to the whole image.

        This routine is designed to be fast even with large images and to produce
        minimal file sizes when the output is saved in a vector format like PDF.

        """
        if pixelsize is None:
            pixelsize = np.min(distance.pdist(np.column_stack([x, y])))

        xmin, xmax = np.min(x), np.max(x)
        ymin, ymax = np.min(y), np.max(y)
        nx = round((xmax - xmin)/pixelsize) + 1
        ny = round((ymax - ymin)/pixelsize) + 1
        j = np.round((x - xmin)/pixelsize).astype(int)
        k = np.round((y - ymin)/pixelsize).astype(int)
        mask = np.ones((nx, ny), dtype=bool)
        img = np.empty((nx, ny))
        mask[j, k] = 0
        img[j, k] = val
        img = np.ma.masked_array(img, mask)

        ax = plt.gca()

        if (angle is None) or (angle == 0):

            f = ax.imshow(np.rot90(img), interpolation='none',
                          cmap=kwargs.get("cmap", self.sauron),
                          extent=[xmin-pixelsize/2, xmax+pixelsize/2,
                                  ymin-pixelsize/2, ymax+pixelsize/2])

        else:

            x, y = np.ogrid[xmin-pixelsize/2 : xmax+pixelsize/2 : (nx+1)*1j,
                            ymin-pixelsize/2 : ymax+pixelsize/2 : (ny+1)*1j]
            ang = np.radians(angle)
            x, y = x*np.cos(ang) - y*np.sin(ang), x*np.sin(ang) + y*np.cos(ang)

            mask1 = np.ones_like(x, dtype=bool)
            mask1[:-1, :-1] *= mask  # Flag the four corners of the mesh
            mask1[:-1, 1:] *= mask
            mask1[1:, :-1] *= mask
            mask1[1:, 1:] *= mask
            x = np.ma.masked_array(x, mask1)  # Mask is used for proper plot range
            y = np.ma.masked_array(y, mask1)

            f = ax.pcolormesh(x, y, img, cmap=kwargs.get("cmap", self.sauron))
            ax.axis('image')

        ax.minorticks_on()
        ax.tick_params(length=10, width=1, which='major')
        ax.tick_params(length=5, width=1, which='minor')

        return f

    def display_bins(self, x, y, binNum, velBin):
        """
        NAME:
            display_bins()
            
        AUTHOR:
            Michele Cappellari, University of Oxford
            cappellari_at_astro.ox.ac.uk

        PURPOSE:
            This simple routine illustrates how to display a Voronoi binned map.
            
        INPUTS:
            (x, y): (length npix) Coordinates of the original spaxels before binning;
            binNum: (length npix) Bin number corresponding to each (x, y) pair,
                    as provided in output by the voronoi_2d_binning() routine;
            velBin: (length nbins) Quantity associated to each bin, resulting
                    e.g. from the kinematic extraction from the binned spectra.
                  
        MODIFICATION HISTORY:
            V1.0.0: Michele Cappellari, Oxford, 15 January 2015          
        
        """
        npix = len(binNum)
        if (npix != len(x)) or (npix != len(y)):
            raise ValueError('The vectors (x, y, binNum) must have the same size')
            
        f = self.display_pixels(x, y, velBin[binNum])
        
        return f

    def fit_kinematic_pa(self, debug=False, nsteps=361, 
                         quiet=False, plot=True, dvel=None):
        """
             NAME:
           FIT_KINEMATIC_PA

         PURPOSE:
           Determine the global kinematic position angle of a
           galaxy with the method described in Appendix C of
           Krajnovic, Cappellari, de Zeeuw, & Copin 2006, MNRAS, 366, 787


         INPUT PARAMETERS:
           XBIN, YBIN: vectors with the coordinates of the bins (or pixels)
               measured from the centre of rotation (typically the galaxy centre).
             - IMPORTANT: The routine will not give meaningful output unless 
               (X,Y)=(0,0) is an estimate of the centre of rotation.        
           VEL: measured velocity at the position (XBIN,YBIN). 
             - IMPORTANT: An estimate of the systemic velocity has to be already 
               subtracted from this velocity [e.g. VEL = VEL - median(VEL)]. 
               The routine will then provide in the output VELSYST a correction 
               to be added to the adopted systemic velocity.

         INPUT KEYWORDS:
           NSTEPS: number of steps along which the angle is sampled.
               Default is 361 steps which gives a 0.5 degr accuracy.
               Decrease this number to limit computation time during testing.

         OUTPUT PARAMETER:
           ANGLEBEST: kinematical PA. Note that this is the angle along which
               |Vel| is maximum (note modulus!). If one reverses the sense of
               rotation in a galaxy ANGLEBEST does not change. The sense of
               rotation can be trivially determined by looking at the map of Vel.
           ANGLEERROR: corresponding error to assign to ANGLEBEST.
           VELSYST: Best-fitting correction to the adopted systemic velocity 
               for the galaxy.
             - If the median was subtracted to the input velocity VEL before 
               the PA fit, then the corrected systemnic velocity will be 
               median(VEL)+VELSYST.

         REQUIRED ROUTINES:
           The following five additional routines are needed:
           - 1. CAP_SYMMETRIZE_VELFIELD and 2. CAP_RANGE: by Michele Cappellari
             (included in this FIT_KINEMATIC_PA distribution)
           - 3. SAURON_COLORMAP and 4. PLOT_VELFIELD: can be obtained from:
             http://purl.org/cappellari/idl#binning
           - 5. SIGRANGE: from IDL astro library http://idlastro.gsfc.nasa.gov/

        """
        vel = copy(self.vel_flat)
        if dvel is None:
            dvel = vel*0 + 10.0 # Adopt here constant 10 km/s errors!
        
        nbins = self.xbin.size
        n = nsteps
        angles = np.linspace(0, 180, n) # 0.5 degrees steps by default
        chi2 = np.empty_like(angles)
        for j, ang in enumerate(angles):
            velSym = self.symmetrize_velfield(self.xbin, self.ybin, vel, sym=1, pa=ang)
            chi2[j] = np.nansum(((vel-velSym)/dvel)**2)
            if debug:
                print('Ang, chi2/DOF:', ang, chi2[j]/nbins)
                self.plot_velfield(self.xbin, self.ybin, velSym)
                plt.pause(0.01)
        k = np.argmin(chi2)
        angBest = angles[k]
        
        # Compute fit at the best position
        #
        velSym = self.symmetrize_velfield(self.xbin, self.ybin, vel, sym=1, pa=angBest)
        if angBest < 0:
            angBest += 180
        
        # 3sigma confidence limit, including error on chi^2
        #
        f = chi2 - chi2[k] <= 9 + 3*np.sqrt(2*nbins)
        if f.sum():
            angErr = (np.max(angles[f]) - np.min(angles[f]))/2.0
            if angErr >= 45:
                good = np.degrees(np.arctan(np.tan(np.radians(angles[f]))))
                angErr = (np.max(good) - np.min(good))/2.0
        else:
            angErr = max(0.5, (angles[1]-angles[0])/2.0)
        
        # angErr = angErr.clip(max(0.5, (angles[1]-angles[0])/2.0)) # Force errors to be larger than 0.5 deg
        vSyst = np.nanmedian(vel - velSym)
        
        if not quiet:
            print('  Kin PA:', angBest, ' +/- ', angErr, ' (3*sigma error)')
            print('Velocity Offset:', vSyst)
        
        # Plot results
        #
        if plot:    
        
            mn, mx = stats.scoreatpercentile(velSym, [2.5, 97.5])
            mx = np.nanmin([mx, -mn])

            plt.subplot(121)
            self.plot_velfield(self.xbin, self.ybin, velSym, vmin=-mx, vmax=mx) 
            plt.title('Symmetrized')
            
            # debugging 
            print (velSym)
            print (np.nanmin(vel - vSyst), np.nanmax(vel - vSyst))

            plt.subplot(122)
            self.plot_velfield(self.xbin, self.ybin, velSym, vmin=-mx, vmax=mx) 
            plt.title('Data and best PA')
            rad = np.sqrt(np.max(self.xbin**2 + self.ybin**2))
            ang = [0,np.pi] + np.radians(angBest)
            plt.plot(rad*np.cos(ang), rad*np.sin(ang), '--', linewidth=3) # Zero-velocity line
            plt.plot(-rad*np.sin(ang), rad*np.cos(ang), linewidth=3) # Major axis PA

        plt.show()
        return angBest, angErr, vSyst