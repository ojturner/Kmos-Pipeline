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
import pickle
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
from photutils import CircularAperture
from photutils import aperture_photometry
import scipy.optimize as op
import emcee
import corner

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

        self.vel_data = self.Table[0].data

        # also load the associated velocity field errors

        try:

            self.er_table = fits.open('%serror_field.fits' % self.fileName[:-14])

            self.error_data = self.er_table[0].data

        except IOError:

            print 'No associated error field'

        try:

            self.sig_table = fits.open('%ssig_field.fits' % self.fileName[:-14])

            self.sig_data = self.sig_table[0].data

        except IOError:

            print 'No associated error field'

        try:

            self.sig_er_table = fits.open('%ssig_error_field.fits' % self.fileName[:-14])

            self.sig_error_data = self.sig_er_table[0].data

        except IOError:

            print 'No associated error field'

        # variables housing the pickled chain and lnp

        self.chain_name = self.fileName[:-5] + '_chain.obj'

        self.ln_p_name = self.fileName[:-5] + '_lnp.obj'

        self.param_file = self.fileName[:-5] + '_params.txt'

        #initialise x and y dimensions

        self.xpix = self.vel_data.shape[0]

        self.ypix = self.vel_data.shape[1]

        # set the lengths of the bins as x & y

        xbin = np.arange(0, self.vel_data.shape[0], 1)

        ybin = np.arange(0, self.vel_data.shape[1], 1)

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

            self.vel_flat.append(self.vel_data[x][y])

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

    def disk_function(self,
                      theta,
                      xpos,
                      ypos):
        """
        Def: Function to calculate disk velocity given input values.
        Note that all angles must be given in radians
        """
        # unpack the parameters

        xcen, ycen, inc, pa, rt, vasym = theta

        # look at the difference between central pixel and pixel
        # being modelled

        diff_x = (xcen - xpos) * 1.0

        diff_y = (ycen - ypos) * 1.0

        # print diff_x, diff_y

        # calculate the pixel angle

        if diff_y == 0 and diff_x != 0:

            pixel_angle = np.arctan(np.sign(diff_x)*np.inf)

            # print 'This is the pixel angle %s' % pixel_angle

        elif diff_y == 0 and diff_x == 0:

            # print 'In the middle'

            pixel_angle = 0.0

        else:

            # print 'computing pixel angle'

            pixel_angle = np.arctan(diff_x / diff_y)

            # print 'pixel angle %s' % (pixel_angle * 180 / np.pi)

        # work out phi which is the overall angle between
        # the spaxel being modelled and the central spaxel/position angle
        # this involves summing with a rotation angle which depends on
        # the spaxel quadrant

        if diff_x >= 0 and diff_y >= 0 and not(diff_x == 0 and diff_y == 0):

            # print 'top left'
            # we're in the upper left quadrant, want rot to be 270

            rot = 3 * np.pi / 2

        elif diff_x >= 0 and diff_y < 0:

            # print 'top right'

            # we're in the upper right quandrant, want rot to be 90

            rot = np.pi / 2

        elif diff_x < 0 and diff_y < 0:

            # print 'lower right'

            # we're in the lower right quadrant

            rot = np.pi / 2

        elif diff_x < 0 and diff_y >= 0:

            # print 'lower left'

            # we're in the lower left quadrant

            rot = 3 * np.pi / 2

        elif diff_x == 0 and diff_y == 0:

            # print 'middle'

            # we're in the middle

            rot = pa

        phi = pixel_angle - pa + rot

    #    print 'differences: %s %s' % (diff_x, diff_y)
    #    print 'pixel angle %s' % (pixel_angle * 180 / np.pi)
    #    print 'position angle %s' % (pa * 180 / np.pi)
    #    print 'rotation angle %s' % (rot * 180 / np.pi)
    #    print 'overall angle %s' % (phi * 180 / np.pi)
    #    print 'cosine of angle %s' % (np.cos(phi))

        r = np.sqrt(diff_x*diff_x + diff_y*diff_y)

        vel = np.cos(phi) * np.sin(inc) * (2 / np.pi) * vasym * np.arctan(r / rt)

        # print vel, xpix, ypix

        return vel

    def grid(self):

        """
        Def: return an empty grid with the specified dimensions
        """

        # create a 1D arrays of length dim_x * dim_y containing the 
        # spaxel coordinates

        xbin = np.arange(0, self.xpix, 1)

        ybin = np.arange(0, self.ypix, 1)

        ybin, xbin = np.meshgrid(ybin, xbin)

        xbin = np.ravel(xbin)

        ybin = np.ravel(ybin)

        return np.array(xbin) * 1.0, np.array(ybin) * 1.0

    def compute_model_grid(self,
                           theta):

        """
        Def:
        Use the grid function to construct a basis for the model.
        Then apply the disk function to each spaxel in the basis
        reshape back to 2d array and plot the model velocity
        """

        xbin, ybin = self.grid()

        # setup list to house the velocity measurements

        vel_array = []

        # compute the model at each spaxel location

        for xpos, ypos in zip(xbin, ybin):

            # run the disk function

            vel_array.append(self.disk_function(theta,
                                                xpos,
                                                ypos))

        # create numpy array from the vel_array list 

        vel_array = np.array(vel_array)

        # reshape back to the chosen grid dimensions

        vel_2d = vel_array.reshape((self.xpix, self.ypix))

        # plot as a 2d array

        #fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        #im = ax.imshow(vel_2d,
        #               cmap=plt.get_cmap('jet'),
        #               interpolation='nearest')
        # add colourbar to each plot
        #divider = make_axes_locatable(ax)
        #cax_new = divider.append_axes('right', size='10%', pad=0.05)
        #plt.colorbar(im, cax=cax_new)
        # set the title
        #ax.set_title('model velocity')
        # plt.show()
        plt.close('all')

        return vel_2d

    def lnlike(self, 
               theta):
        """
        Def: Return the log likelihood for the velocity field function.
        All that has to be done is to compute the model in a grid the same size
        as the data and then plug into the standard likelihood formula.

        Input:
                vel_data - the actual velocity field unicodedata
                vel_errors - the velocity field error grid
                theta - list of parameter values to be fed into the model

        Output:
                some single numerical value for the log likelihood
        """
        # sometimes nice to see what parameters are being tried in the
        # MCMC step

        print theta

        # compute the model grid

        model = self.compute_model_grid(theta)

        # find the grid of inverse sigma values

        inv_sigma2 = 1.0 / (self.error_data * self.error_data)

        ans = -0.5 * (np.nansum((self.vel_data - model)*(self.vel_data - model) *
                                inv_sigma2 - np.log(inv_sigma2)))

        print ans

        return ans

    def lnprior(self,
                theta):

        """
        Def:
        Set an uninformative prior distribution for the parameters in the model
        """

        xcen, ycen, inc, pa, rt, vasym = theta

        if 0 < xcen < 30.0 and \
           0 < ycen < 30.0 and \
           0 < inc < np.pi / 2 and \
           0 < pa < 2 * np.pi and \
           0 < rt < 20.0 and \
           0 < vasym < 300:

            return 0.0

        return -np.inf

    def lnprob(self,
               theta):

        lp = self.lnprior(theta)

        if not np.isfinite(lp):

            return -np.inf

        return lp + self.lnlike(theta)

    def run_emcee(self,
                  theta,
                  nsteps,
                  nwalkers,
                  burn_no):

        ndim = len(theta)

        pos = [theta + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

        sampler = emcee.EnsembleSampler(nwalkers,
                                        ndim,
                                        self.lnprob)

        sampler.run_mcmc(pos,
                         nsteps)

        samples = sampler.chain[:, burn_no:, :].reshape((-1, ndim))

        fig = corner.corner(samples,
                            labels=["$xcen$",
                                    "$ycen$",
                                    "$inc$",
                                    "$pa$",
                                    "rt",
                                    "vasym"],
                            truths=theta)

        fig.savefig('%s_corner_plot.png' % self.fileName[:-5])

        # plt.show()

        # print samples
        # going to save pickled versions of the chain and the lnprobability
        # so that these can be accessed again later if necessary

        if os.path.isfile(self.chain_name):

            os.system('rm %s' % self.chain_name)

        if os.path.isfile(self.ln_p_name):

            os.system('rm %s' % self.ln_p_name)

        chain_file = open(self.chain_name, 'w')

        pickle.dump(sampler.chain, chain_file)

        chain_file.close()

        ln_p_file = open(self.ln_p_name, 'w')

        pickle.dump(sampler.lnprobability, ln_p_file)

        ln_p_file.close()

        # now use the helper function below to open up
        # the pickled files and write to params file

        self.write_params(burn_no)

        # set a variable to the log probability value

    def write_params(self,
                     burn_no):

        """
        Def:
        Helper function to open up the pickled chain and lnp files
        and write the maximum likelihood, 50th percentile, 16th per and 84th
        per parameters to file for ease of application later on

        Input:
                burn_no - number of entries (steps) to burn from the chain
        """

        chain = pickle.load(open(self.chain_name, 'r'))

        lnp = pickle.load(open(self.ln_p_name, 'r'))

        samples = chain[:, burn_no:, :].reshape((-1, chain.shape[2]))

        # initialise the parameter names

        param_names = ['type',
                       'centre_x',
                       'centre_y',
                       'inclination',
                       'position_angle',
                       'Flattening_Radius',
                       'Flattening_Velocity']

        # find the max likelihood parameters

        max_p = np.unravel_index(lnp.argmax(), lnp.shape)

        max_params = chain[max_p[0], max_p[1], :]

        x_mcmc, y_mcmc, i_mcmc, \
            pa_mcmc, rt_mcmc, va_mcmc \
            = zip(*np.percentile(samples, [16, 50, 84],
                  axis=0))

        param_file = self.fileName[:-5] + '_params.txt'

        if os.path.isfile(param_file):

            os.system('rm %s' % param_file)

        # write all of these values to file

        with open(param_file, 'a') as f:

            for item in param_names:

                f.write('%s\t' % item)

            f.write('\nMAX_lnp:\t')

            for item in max_params:

                f.write('%s\t' % item)

            f.write('\n50th_lnp:\t%s\t%s\t%s\t%s\t%s\t%s\t' % (x_mcmc[1],
                                                               y_mcmc[1],
                                                               i_mcmc[1],
                                                               pa_mcmc[1],
                                                               rt_mcmc[1],
                                                               va_mcmc[1]))

            f.write('\n16th_lnp:\t%s\t%s\t%s\t%s\t%s\t%s\t' % (x_mcmc[0],
                                                               y_mcmc[0],
                                                               i_mcmc[0],
                                                               pa_mcmc[0],
                                                               rt_mcmc[0],
                                                               va_mcmc[0]))

            f.write('\n84th_lnp:\t%s\t%s\t%s\t%s\t%s\t%s\t' % (x_mcmc[2],
                                                               y_mcmc[2],
                                                               i_mcmc[2],
                                                               pa_mcmc[2],
                                                               rt_mcmc[2],
                                                               va_mcmc[2]))


    def plot_comparison(self):

        """
        Def:
        Plot the best fitting model alongside the original velocity field
        with position angle and morphological angle also plotted

        Input:
                theta - the now best fit set of parameters
                vel_data - the velocity field unicodedata
                vel_errors - the velocity field errors

        """

        # load in the file

        param_file = np.genfromtxt(self.param_file)

        theta_max = param_file[1][1:]

        theta_50 = param_file[2][1:]

        theta_16 = param_file[3][1:]

        theta_84 = param_file[4][1:]

        # compute the model grid with the specified parameters

        model_max = self.compute_model_grid(theta_max)

        model_50 = self.compute_model_grid(theta_50)

        model_16 = self.compute_model_grid(theta_16)

        model_84 = self.compute_model_grid(theta_84)

        # only want to see the evaluated model at the grid points
        # where the data is not nan. Loop round the data and create
        # a mask which can multiply the model

        mask_array = np.empty(shape=(self.xpix, self.ypix))

        for i in range(0, self.xpix):

            for j in range(0, self.ypix):

                if np.isnan(self.vel_data[i][j]):

                    mask_array[i][j] = np.nan

                else:

                    mask_array[i][j] = 1.0

        # take product of model and mask_array to return new data

        trunc_model_max = mask_array * model_max

        trunc_model_50 = mask_array * model_50

        trunc_model_16 = mask_array * model_16

        trunc_model_84 = mask_array * model_84

        # plot the results

        vel_min, vel_max = np.nanpercentile(self.vel_data,
                                            [5.0, 95.0])

        mod_min, mod_max = np.nanpercentile(trunc_model_max,
                                            [5.0, 95.0])

        plt.close('all')

        fig, ax = plt.subplots(1, 2, figsize=(12, 6))

        im = ax[0].imshow(self.vel_data,
                          cmap=plt.get_cmap('jet'),
                          vmin=mod_min,
                          vmax=mod_max,
                          interpolation='nearest')

        # add colourbar to each plot
        divider = make_axes_locatable(ax[0])
        cax_new = divider.append_axes('right', size='10%', pad=0.05)
        plt.colorbar(im, cax=cax_new)

        # set the title
        ax[0].set_title('[OIII] Velocity Data')

        im = ax[1].imshow(trunc_model_50,
                          cmap=plt.get_cmap('jet'),
                          vmin=mod_min,
                          vmax=mod_max,
                          interpolation='nearest')

        # add colourbar to each plot
        divider = make_axes_locatable(ax[1])
        cax_new = divider.append_axes('right', size='10%', pad=0.05)
        plt.colorbar(im, cax=cax_new)

        # set the title
        ax[1].set_title('[OIII] Velocity Model')

        # plt.show()

        fig.savefig('%s_model_comparison.png' % self.fileName[:-5])

        plt.close('all')

    def extract_in_apertures(self,
                             r_aper,
                             d_aper):

        """
        Def: Extract the velocity field along the kinematic axis returned by the
        model fitting in both the data and the model for comparison. The model
        will show a perfect arctangent function.

        Input:
                theta - array of best fitting model parameter values
                model_data - best fit model computed from the compute_model_grid
                vel_data - array containing the actual velocity data
                r_aper - aperture size in pixels to use for each aperture
                d_aper - distance spacing between apertures
        Output:
                1D arrays containing the extracted model and real velocity fields
                along the kinematic major axis
        """

        # assign the best fit parameters to variables from the theta array

        # load in the file

        param_file = np.genfromtxt(self.param_file)

        theta_max = param_file[1][1:]

        pa_max = theta_max[3]

        xcen_max = theta_max[0]

        ycen_max = theta_max[1]

        theta_50 = param_file[2][1:]

        pa_50 = theta_50[3]

        xcen_50 = theta_50[0]

        ycen_50 = theta_50[1]

        theta_16 = param_file[3][1:]

        pa_16 = theta_16[3]

        xcen_16 = theta_16[0]

        ycen_16 = theta_16[1]

        theta_84 = param_file[4][1:]

        pa_84 = theta_84[3]

        xcen_84 = theta_84[0]

        ycen_84 = theta_84[1]

        # compute the model grid with the specified parameters

        model_max = self.compute_model_grid(theta_max)

        model_50 = self.compute_model_grid(theta_50)

        model_16 = self.compute_model_grid(theta_16)

        model_84 = self.compute_model_grid(theta_84)

        # initialise the list of aperture positions with the xcen and ycen

        positions_max = []

        positions_50 = []

        positions_16 = []

        positions_84 = []

        # first job is to compute the central locations of the apertures
        # do this by fixing the distance along the KA between aperture centres

        xdim = self.xpix - 2

        ydim = self.ypix - 2

        # find the steps along the KA with which to increment

        x_inc_max = d_aper * abs(np.sin((np.pi / 2.0) - pa_max))

        y_inc_max = d_aper * abs(np.cos((np.pi / 2.0) - pa_max))

        x_inc_50 = d_aper * abs(np.sin((np.pi / 2.0) - pa_50))

        y_inc_50 = d_aper * abs(np.cos((np.pi / 2.0) - pa_50))

        x_inc_16 = d_aper * abs(np.sin((np.pi / 2.0) - pa_16))

        y_inc_16 = d_aper * abs(np.cos((np.pi / 2.0) - pa_16))

        x_inc_84 = d_aper * abs(np.sin((np.pi / 2.0) - pa_84))

        y_inc_84 = d_aper * abs(np.cos((np.pi / 2.0) - pa_84))

        # now find the sequence of aperture centres up until the boundaries
        # this is tricky - depending on the PA need to increase and decrease
        # both x and y together, or increase one and decrease the other

        # statements for the maximum likelihood position angle

        if 0 < pa_max < np.pi / 2.0 or np.pi < pa_max < 3 * np.pi / 2.0:

            print 'Top Right and Bottom Left'

            # need to increase x and decrease y and vice versa

            new_x_max = xcen_max + x_inc_max

            new_y_max = ycen_max - y_inc_max

            # while loop until xdim is breached or 0 is breached for y

            while new_x_max < xdim and new_y_max > 2:

                # append aperture centre to the positions array

                positions_max.append((new_y_max, new_x_max))

                new_x_max += x_inc_max

                new_y_max -= y_inc_max

                # print new_x_max, new_y_max

            # starting from the left so need to reverse list direction
            # and append the central point

            positions_max = positions_max[::-1]

            positions_max.append((ycen_max, xcen_max))

            # now go in the other direction

            new_x_max = xcen_max - x_inc_max

            new_y_max = ycen_max + y_inc_max

            # while loop until xdim is breached or 0 is breached for y

            while new_x_max > 2 and new_y_max < ydim:

                # append aperture centre to the positions_max array

                positions_max.append((new_y_max, new_x_max))

                new_x_max -= x_inc_max

                new_y_max += y_inc_max

                # print new_x, new_y_max

        # deal with the other cases of position angle

        else:

            print 'Top Left and Bottom Right'

            # need to increase x and increase y and vice versa

            new_x_max = xcen_max - x_inc_max

            new_y_max = ycen_max - y_inc_max

            # while loop until xdim is 2 or ydim is 2

            while new_x_max > 2 and new_y_max > 2:

                # append aperture centre to the positions_max array

                positions_max.append((new_y_max, new_x_max))

                new_x_max -= x_inc_max

                new_y_max -= y_inc_max

            # starting from the left so need to reverse list direction
            # and append the central point

            positions_max = positions_max[::-1]

            positions_max.append((ycen_max, xcen_max))

            # now go in the other direction

            new_x_max = xcen_max + x_inc_max

            new_y_max = ycen_max + y_inc_max

            # while loop until xdim is breached or ydim is breached

            while new_x_max < xdim and new_y_max < ydim:

                # append aperture centre to the positions_max array

                positions_max.append((new_y_max, new_x_max))

                new_x_max += x_inc_max

                new_y_max += y_inc_max

        # statements for the 50th percentile position angle

        if 0 < pa_50 < np.pi / 2.0 or np.pi < pa_50 < 3 * np.pi / 2.0:

            print 'Top Right and Bottom Left'

            # need to increase x and decrease y and vice versa

            new_x_50 = xcen_50 + x_inc_50

            new_y_50 = ycen_50 - y_inc_50

            # while loop until xdim is breached or 0 is breached for y

            while new_x_50 < xdim and new_y_50 > 2:

                # append aperture centre to the positions array

                positions_50.append((new_y_50, new_x_50))

                new_x_50 += x_inc_50

                new_y_50 -= y_inc_50

                # print new_x_50, new_y_50

            # starting from the left so need to reverse list direction
            # and append the central point

            positions_50 = positions_50[::-1]

            positions_50.append((ycen_50, xcen_50))

            # now go in the other direction

            new_x_50 = xcen_50 - x_inc_50

            new_y_50 = ycen_50 + y_inc_50

            # while loop until xdim is breached or 0 is breached for y

            while new_x_50 > 2 and new_y_50 < ydim:

                # append aperture centre to the positions_50 array

                positions_50.append((new_y_50, new_x_50))

                new_x_50 -= x_inc_50

                new_y_50 += y_inc_50

                # print new_x, new_y_50

        # deal with the other cases of position angle

        else:

            print 'Top Left and Bottom Right'

            # need to increase x and increase y and vice versa

            new_x_50 = xcen_50 - x_inc_50

            new_y_50 = ycen_50 - y_inc_50

            # while loop until xdim is 2 or ydim is 2

            while new_x_50 > 2 and new_y_50 > 2:

                # append aperture centre to the positions_50 array

                positions_50.append((new_y_50, new_x_50))

                new_x_50 -= x_inc_50

                new_y_50 -= y_inc_50

            # starting from the left so need to reverse list direction
            # and append the central point

            positions_50 = positions_50[::-1]

            positions_50.append((ycen_50, xcen_50))

            # now go in the other direction

            new_x_50 = xcen_50 + x_inc_50

            new_y_50 = ycen_50 + y_inc_50

            # while loop until xdim is breached or ydim is breached

            while new_x_50 < xdim and new_y_50 < ydim:

                # append aperture centre to the positions_50 array

                positions_50.append((new_y_50, new_x_50))

                new_x_50 += x_inc_50

                new_y_50 += y_inc_50

        # statements for the 16th percentile position angle

        if 0 < pa_16 < np.pi / 2.0 or np.pi < pa_16 < 3 * np.pi / 2.0:

            print 'Top Right and Bottom Left'

            # need to increase x and decrease y and vice versa

            new_x_16 = xcen_16 + x_inc_16

            new_y_16 = ycen_16 - y_inc_16

            # while loop until xdim is breached or 0 is breached for y

            while new_x_16 < xdim and new_y_16 > 2:

                # append aperture centre to the positions array

                positions_16.append((new_y_16, new_x_16))

                new_x_16 += x_inc_16

                new_y_16 -= y_inc_16

                # print new_x_16, new_y_16

            # starting from the left so need to reverse list direction
            # and append the central point

            positions_16 = positions_16[::-1]

            positions_16.append((ycen_16, xcen_16))

            # now go in the other direction

            new_x_16 = xcen_16 - x_inc_16

            new_y_16 = ycen_16 + y_inc_16

            # while loop until xdim is breached or 0 is breached for y

            while new_x_16 > 2 and new_y_16 < ydim:

                # append aperture centre to the positions_16 array

                positions_16.append((new_y_16, new_x_16))

                new_x_16 -= x_inc_16

                new_y_16 += y_inc_16

                # print new_x, new_y_16

        # deal with the other cases of position angle

        else:

            print 'Top Left and Bottom Right'

            # need to increase x and increase y and vice versa

            new_x_16 = xcen_16 - x_inc_16

            new_y_16 = ycen_16 - y_inc_16

            # while loop until xdim is 2 or ydim is 2

            while new_x_16 > 2 and new_y_16 > 2:

                # append aperture centre to the positions_16 array

                positions_16.append((new_y_16, new_x_16))

                new_x_16 -= x_inc_16

                new_y_16 -= y_inc_16

            # starting from the left so need to reverse list direction
            # and append the central point

            positions_16 = positions_16[::-1]

            positions_16.append((ycen_16, xcen_16))

            # now go in the other direction

            new_x_16 = xcen_16 + x_inc_16

            new_y_16 = ycen_16 + y_inc_16

            # while loop until xdim is breached or ydim is breached

            while new_x_16 < xdim and new_y_16 < ydim:

                # append aperture centre to the positions_16 array

                positions_16.append((new_y_16, new_x_16))

                new_x_16 += x_inc_16

                new_y_16 += y_inc_16

        # statements for the 84th percenntile position angle

        if 0 < pa_84 < np.pi / 2.0 or np.pi < pa_84 < 3 * np.pi / 2.0:

            print 'Top Right and Bottom Left'

            # need to increase x and decrease y and vice versa

            new_x_84 = xcen_84 + x_inc_84

            new_y_84 = ycen_84 - y_inc_84

            # while loop until xdim is breached or 0 is breached for y

            while new_x_84 < xdim and new_y_84 > 2:

                # append aperture centre to the positions array

                positions_84.append((new_y_84, new_x_84))

                new_x_84 += x_inc_84

                new_y_84 -= y_inc_84

                # print new_x_84, new_y_84

            # starting from the left so need to reverse list direction
            # and append the central point

            positions_84 = positions_84[::-1]

            positions_84.append((ycen_84, xcen_84))

            # now go in the other direction

            new_x_84 = xcen_84 - x_inc_84

            new_y_84 = ycen_84 + y_inc_84

            # while loop until xdim is breached or 0 is breached for y

            while new_x_84 > 2 and new_y_84 < ydim:

                # append aperture centre to the positions_84 array

                positions_84.append((new_y_84, new_x_84))

                new_x_84 -= x_inc_84

                new_y_84 += y_inc_84

                # print new_x, new_y_84

        # deal with the other cases of position angle

        else:

            print 'Top Left and Bottom Right'

            # need to increase x and increase y and vice versa

            new_x_84 = xcen_84 - x_inc_84

            new_y_84 = ycen_84 - y_inc_84

            # while loop until xdim is 2 or ydim is 2

            while new_x_84 > 2 and new_y_84 > 2:

                # append aperture centre to the positions_84 array

                positions_84.append((new_y_84, new_x_84))

                new_x_84 -= x_inc_84

                new_y_84 -= y_inc_84

            # starting from the left so need to reverse list direction
            # and append the central point

            positions_84 = positions_84[::-1]

            positions_84.append((ycen_84, xcen_84))

            # now go in the other direction

            new_x_84 = xcen_84 + x_inc_84

            new_y_84 = ycen_84 + y_inc_84

            # while loop until xdim is breached or ydim is breached

            while new_x_84 < xdim and new_y_84 < ydim:

                # append aperture centre to the positions_84 array

                positions_84.append((new_y_84, new_x_84))

                new_x_84 += x_inc_84

                new_y_84 += y_inc_84

        # positions array should now be populated with all of the apertures

        # print positions

        # now perform aperture photometry on the model data to check that this
        # actually works. Remember that the velocity computed for each
        # aperture will be the sum returned divided by the area

        pixel_area = np.pi * r_aper * r_aper

        # the max likelihood extraction parameters

        apertures_max = CircularAperture(positions_max, r=r_aper)

        mod_phot_table_max = aperture_photometry(model_max, apertures_max)

        real_phot_table_max = aperture_photometry(self.vel_data, apertures_max)

        real_error_table_max = aperture_photometry(self.error_data, apertures_max)

        sig_table_max = aperture_photometry(self.sig_data, apertures_max)

        sig_error_table_max = aperture_photometry(self.sig_error_data, apertures_max)

        mod_velocity_values_max = mod_phot_table_max['aperture_sum'] / pixel_area

        real_velocity_values_max = real_phot_table_max['aperture_sum'] / pixel_area

        real_error_values_max = real_error_table_max['aperture_sum'] / pixel_area

        sig_values_max = sig_table_max['aperture_sum'] / pixel_area

        sig_error_values_max = sig_error_table_max['aperture_sum'] / pixel_area

        x_max = np.arange(0, len(mod_phot_table_max['aperture_sum']), 1)

        # the 50th percentile extraction parameters

        apertures_50 = CircularAperture(positions_50, r=r_aper)

        mod_phot_table_50 = aperture_photometry(model_50, apertures_50)

        real_phot_table_50 = aperture_photometry(self.vel_data, apertures_50)

        real_error_table_50 = aperture_photometry(self.error_data, apertures_50)

        sig_table_50 = aperture_photometry(self.sig_data, apertures_50)

        sig_error_table_50 = aperture_photometry(self.sig_error_data, apertures_50)

        mod_velocity_values_50 = mod_phot_table_50['aperture_sum'] / pixel_area

        real_velocity_values_50 = real_phot_table_50['aperture_sum'] / pixel_area

        real_error_values_50 = real_error_table_50['aperture_sum'] / pixel_area

        sig_values_50 = sig_table_50['aperture_sum'] / pixel_area

        sig_error_values_50 = sig_error_table_50['aperture_sum'] / pixel_area

        x_50 = np.arange(0, len(mod_phot_table_50['aperture_sum']), 1)

        # the 16th percentile extraction parameters

        apertures_16 = CircularAperture(positions_16, r=r_aper)

        mod_phot_table_16 = aperture_photometry(model_16, apertures_16)

        real_phot_table_16 = aperture_photometry(self.vel_data, apertures_16)

        real_error_table_16 = aperture_photometry(self.error_data, apertures_16)

        sig_table_16 = aperture_photometry(self.sig_data, apertures_16)

        sig_error_table_16 = aperture_photometry(self.sig_error_data, apertures_16)

        mod_velocity_values_16 = mod_phot_table_16['aperture_sum'] / pixel_area

        real_velocity_values_16 = real_phot_table_16['aperture_sum'] / pixel_area

        real_error_values_16 = real_error_table_16['aperture_sum'] / pixel_area

        sig_values_16 = sig_table_16['aperture_sum'] / pixel_area

        sig_error_values_16 = sig_error_table_16['aperture_sum'] / pixel_area

        x_16 = np.arange(0, len(mod_phot_table_16['aperture_sum']), 1)

        # the 84th percentile extraction parameters

        apertures_84 = CircularAperture(positions_84, r=r_aper)

        mod_phot_table_84 = aperture_photometry(model_84, apertures_84)

        real_phot_table_84 = aperture_photometry(self.vel_data, apertures_84)

        real_error_table_84 = aperture_photometry(self.error_data, apertures_84)

        sig_table_84 = aperture_photometry(self.sig_data, apertures_84)

        sig_error_table_84 = aperture_photometry(self.sig_error_data, apertures_84)

        mod_velocity_values_84 = mod_phot_table_84['aperture_sum'] / pixel_area

        real_velocity_values_84 = real_phot_table_84['aperture_sum'] / pixel_area

        real_error_values_84 = real_error_table_84['aperture_sum'] / pixel_area

        sig_values_84 = sig_table_84['aperture_sum'] / pixel_area

        sig_error_values_84 = sig_error_table_84['aperture_sum'] / pixel_area

        x_84 = np.arange(0, len(mod_phot_table_84['aperture_sum']), 1)

        # plotting the model and extracted quantities

        fig, ax = plt.subplots(1, 1, figsize=(12, 6))

        ax.plot(x_max,
                mod_velocity_values_max,
                color='red',
                label='max_model')

        ax.errorbar(x_max,
                    real_velocity_values_max,
                    yerr=real_error_values_max,
                    fmt='o',
                    color='red',
                    label='max_data')

        ax.plot(x_50,
                mod_velocity_values_50,
                color='blue',
                label='50_model')

        ax.errorbar(x_50,
                    real_velocity_values_50,
                    yerr=real_error_values_50,
                    fmt='o',
                    color='blue',
                    label='50_data')

        ax.plot(x_16,
                mod_velocity_values_16,
                color='orange',
                linestyle='--',
                label='16_model')

        ax.plot(x_84,
                mod_velocity_values_84,
                color='purple',
                linestyle='--',
                label='84_model')

        ax.legend(prop={'size':10})

        ax.set_title('Model and Real Velocity')

        ax.set_ylabel('velocity (kms$^{-1}$)')

        ax.set_xlabel('aperture number')

        # plt.show()

        fig.savefig('%s_1d_velocity_plot.png' % self.fileName[:-5])

        plt.close('all')

        fig, ax = plt.subplots(1, 1, figsize=(12, 6))

        ax.errorbar(x_max,
                    sig_values_max,
                    yerr=sig_error_values_max,
                    fmt='o',
                    color='red',
                    label='max_data')

        ax.errorbar(x_50,
                    sig_values_50,
                    yerr=sig_error_values_50,
                    fmt='o',
                    color='blue',
                    label='50_data')

        ax.set_title('Velocity Dispersion')

        ax.set_ylabel('velocity (kms$^{-1}$)')

        ax.set_xlabel('aperture number')

        ax.legend(prop={'size':10})

        # plt.show()

        fig.savefig('%s_1d_dispersion_plot.png' % self.fileName[:-5])