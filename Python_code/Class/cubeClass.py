# This class houses the methods which are relevant to performing manipulations 
#on the reconstructed and combined data cubes. As such, the class will be a datacube  
#object 


# import the relevant modules
import os, sys, numpy as np, random, matplotlib.mlab as mlab, matplotlib.pyplot as plt, math, matplotlib.cm as cm
import pyraf
import numpy.polynomial.polynomial as poly
from numpy import poly1d
import lmfit
import scipy
import math
from lmfit.models import GaussianModel, ExponentialModel, LorentzianModel, VoigtModel, PolynomialModel
from scipy import stats
from scipy import optimize
from lmfit import Model
from scipy.optimize import minimize
from scipy.optimize import basinhopping
from astropy.io import fits
from matplotlib.backends.backend_pdf import PdfPages
from pylab import *
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable


####################################################################

class cubeOps(object):
    """
    Def: 
    Class mainly for combined data cubes output from the
    kmo_sci_red recipe. Contains a series of functions and 
    definitions useful for manipulating the datacubes. 
    Input: 
    sci_combined datacube
    """
    #Initialiser creates an instance of the cube object
    #Input must be a combined data cube with two extensions - data and noise 
    def __init__(self, fileName):
        self.self = self
        self.fileName = fileName
        #Variable containing all the fits extensions 
        self.Table = fits.open(fileName)
        #Variable housing the primary data cube
        self.data = self.Table[1].data
        #Collapse over the wavelength axis to get an image
        self.imData = np.nanmedian(self.data, axis=0)
        try:
            self.total_spec = np.nanmedian(np.nanmedian(self.data[:,4:len(self.data[0]),4:len(self.data[1])], axis=1), axis=1)
        except:
            print 'Cannot extract the total spectrum'
        #Create a plot of the image 
        #Variable housing the noise data cube
        #self.noise = self.Table[2].data
        #Primary Header 
        self.primHeader = self.Table[0].header
        #data Header 
        self.dataHeader = self.Table[1].header
        #noise Header 
        #self.noiseHeader = self.Table[2].header
        #Extract the wavelength calibration info
        try:
            self.start_L = self.dataHeader["CRVAL3"]
            self.dL = self.dataHeader["CDELT3"]
        except KeyError:
            print '[INFO]: This is a raw image'
        #Extract the IFU number from the data 
        #Cube may not have this keyword, try statement
        #Not sure if this is good practice - conditional attribute.
        try:
            self.IFUNR = self.dataHeader["HIERARCH ESO PRO IFUNR"]
            key = 'HIERARCH ESO OCS ARM' + str(self.IFUNR) + ' NAME'
            #print key
            self.IFUName = self.dataHeader[key]
        except KeyError:
            try:
                self.noise_header = self.Table[2].header
                self.IFUNR = self.noise_header["HIERARCH ESO PRO IFUNR"]
                key = 'HIERARCH ESO OCS ARM' + str(self.IFUNR) + ' NAME'
                #print key
                self.IFUName = self.noise_header[key]
            except:

                print '[INFO]: This is not a combined Frame, setting arm name...'
                try:
                    ext_name = self.dataHeader['EXTNAME']
                    num_string = ''
                    for s in ext_name:
                        if s.isdigit():
                            num_string += s
                    num_string = int(num_string)
                    self.IFUNR = copy(num_string)
                    self.IFUName = self.primHeader["HIERARCH ESO OCS ARM" + str(self.IFUNR) + " NAME"]
                    print '[INFO]: You have specified a reconstructed type'
                except KeyError:
                    print "[Warning]: not a datacube"

        #Set the RA and DEC positions of all the arms. These are in 
        #sexagesimal format - convert to degrees for the plot 
        self.raDict = {}
        self.decDict = {}
        self.combDict = {}
        self.offList = []
        for i in range(1, 25):
            # print "HIERARCH ESO OCS ARM" + str(i) + " NAME"
            # print "ARM" + str(i) + "_SCI"

            # if statement to cover the possibility of operational IFUs 
            # assigned to sky 
            if (self.primHeader["HIERARCH ESO OCS ARM" + str(i) + " NAME"]) == ("ARM" + str(i) + "_SCI"):
                # Add the IFU to the offlist
                print 'ARM%s is operational but on sky' % str(i)
                self.offList.append(i)
            else: 
                # Construct the list of combined science frames that will 
                # be thrown out by the science pipeline
                try:
                    nuName = "HIERARCH ESO OCS ARM" + str(i) + " NOTUSED"
                    temp = self.primHeader[nuName]
                except KeyError:
                    combKey = "HIERARCH ESO OCS ARM" + str(i) + " ORIGARM"
                    combName = "HIERARCH ESO OCS ARM" + str(i) + " NAME"
                    self.combDict[self.primHeader[combName]] = self.primHeader[combKey]

                try:      
                    raName = "HIERARCH ESO OCS ARM" + str(i) + " ALPHA"
                    DictName = 'Arm' + str(i)
                    decName = "HIERARCH ESO OCS ARM" + str(i) + " DELTA"
                    #print raName, decName
                    self.raDict[DictName] = self.raToDeg(self.primHeader[raName])
                    self.decDict[DictName] = self.decToDeg(self.primHeader[decName])

                except KeyError:
                    print '[INFO]: IFU %s Not in Use, or not pointing at at object' % DictName
                    self.offList.append(i)




        #Construct the list of combined science names separately
        #This is now in order of the IFU
        self.combNames = []
        for entry in self.combDict.keys():
            combinedName = 'sci_combined_' + entry + '__telluric_skytweak.fits'
            self.combNames.append(combinedName)

        #Also construct the list of kmo_combine recipe combined names 
        self.rec_combNames = []
        for entry in self.combDict.keys():
            combinedName = 'combine_sci_reconstructed_' + entry + '.fits'
            self.rec_combNames.append(combinedName)



        self.offList = np.array(self.offList)
        self.offList = self.offList - 1
        self.raArray = self.raDict.values()
        self.decArray = self.decDict.values()
        self.IFUArms = self.raDict.keys()
        self.xDit = self.primHeader['HIERARCH ESO OCS TARG DITHA']
        self.yDit = self.primHeader['HIERARCH ESO OCS TARG DITHD']

        #Find the pixel scale if this is a combined cube 
        try:
            self.pix_scale = self.primHeader['HIERARCH ESO PRO REC1 PARAM7 VALUE']
        except KeyError:
            print '[INFO]: Could not set pixel scale - not a datacube'
            self.pix_scale = 0

        #Create the wavelength array if this is a combined data type
        try:
            self.wave_array = self.start_L + np.arange(0, 2048*(self.dL), self.dL)
        except:
            print '[INFO]: cannot set wavelength array'
        #Can define all kinds of statistics from the data common to the methods
        #Find the brightest median pixel in the array 
        self.med_array = np.nanmedian(self.data, axis=0)
        self.num = np.nanargmax(self.med_array)
        self.ind1 = self.num / len(self.med_array)
        self.ind2 = self.num - (len(self.med_array) * self.ind1)



    def raToDeg(self, ra):
        """
        Def:
        Helper function - convert sexagesimal RA to degrees. 
        Needs to check number digits before the decimal point, 
        because the fits files doesn't always give 6 

        Input: Sexagesimal RA (HH:MM:SS.SS)
        Output: Ra in degrees 

        """
        ra = str(ra)
        #Figure out how many characters before the decimal point
        i = 0
        for char in ra:
            if char == '.':
                break
            else:
                i += 1
        #Conditionally convert, depending on whether i is 4 or 6
        #Also conditionally execute depending on number decimal places
        if i == 6:
            hours = int(ra[0:2])
            mins = int(ra[2:4])
            secs = float(ra[4:])
            raDeg = (hours * 15) + (mins * 0.25) + (secs * 1.0/240)
        elif i == 4:
            mins = int(ra[0:2])
            secs = float(ra[2:])
            raDeg = (mins * 0.25) + (secs * 1.0/240)            
        else:
            secs = float(ra)
            raDeg = (secs * 1.0/240)
        return raDeg


    def decToDeg(self, dec):
        """
        Def:
        Helper function - convert sexagesimal dec to degrees. 
        Needs to check number digits before the decimal point, 
        because the fits files doesn't always give 6 

        Input: Sexagesimal dec (DD:MM:SS.SS)
        Output: Dec in degrees 

        """     
        dec = str(dec)
        #Figure out how many characters before the decimal point
        i = 0
        if dec[0] == '-':

            for char in dec:
                if char == '.':
                    break
                else:
                    i += 1
            #Conditionally convert, depending on whether i is 4 or 6
            #Also conditionally execute depending on number decimal places
            if i == 7:
                deg = float(dec[1:3])
                mins = float(dec[3:5])
                secs = float(dec[5:])
                #Careful whether deg is negative or not
                #Becoming more positive if > 0
                decDeg = (deg * -1) - (mins * 1.0/60) - (secs * 1.0/3600)
            
            elif i == 5:
                mins = float(dec[1:3])
                secs = float(dec[3:])
                decDeg = (mins * -1.0/60) - (secs * 1.0/3600)                       
            else:
                secs = float(dec)
                decDeg = (secs * 1.0/3600)
            return decDeg

        else:
            
            for char in dec:
                if char == '.':
                    break
                else:
                    i += 1
            #Conditionally convert, depending on whether i is 4 or 6
            #Also conditionally execute depending on number decimal places
            #print i
            if i == 6:
                deg = float(dec[0:2])
                mins = float(dec[2:4])
                secs = float(dec[4:])
                #Careful whether deg is negative or not
                #Becoming more positive if > 0
                decDeg = deg + (mins * 1.0/60) + (secs * 1.0/3600)
            
            elif i == 4:
                mins = float(dec[0:2])
                secs = float(dec[2:])
                decDeg = (mins * 1.0/60) + (secs * 1.0/3600)                        
            else:
                secs = float(dec)
                decDeg = (secs * 1.0/3600)
            return decDeg


    def ifuCoordsPlot(self):
        """
        Def:
        Plots the already recorded IFU Positions on the sky 

        Inputs: None
        Output: Matplotlib plot of coordinates

        """

        plt.scatter(self.raArray, self.decArray)
        plt.show()
        plt.close('all')



    def specPlot(self, gridSize):
        """
        Def:
        Takes a data cube and creates a median stacked 1-D
        spectrum around the brightest pixel, with the size 
        of the median stack defined by gridSize

        Input: 
        gridSize - must be a positive int less than the number of 
        spatial pixels in the cube. 

        Output: 
        matplotlib plot of the 1D stacked spectrum 


        """
        #If choosing to construct the 1D spectrum from a single pixel:
        print '[INFO]: The Brightest Pixel is at: (%s, %s)' % (self.ind1, self.ind2)
        print self.IFUName
        print self.IFUNR
        if gridSize == 1:
            flux_array = self.data[:,self.ind1, self.ind2]
        else:
            lst = []
            for i in range(self.ind1 - (gridSize / 2), self.ind1 + (gridSize / 2)):
                for j in range(self.ind2 - (gridSize / 2), self.ind2 + (gridSize / 2)):
                    lst.append(self.data[:,i,j])
            flux_array = np.nanmedian(lst, axis=0)

        #Now make very basic plot at the moment 
        fig, ax = plt.subplots(1, 1, figsize=(18.0,12.0))
        ax.plot(self.wave_array, flux_array)
        ax.set_title('Flux vs. Wavelength')
        ax.set_xlabel('Wavelength ($\mu m$)')
        ax.set_ylim(0,100)
        saveName = (self.fileName)[:-5] + '.png'
        fig.savefig(saveName)
        #plt.show()     
        plt.close('all')
        return flux_array

    def centralSpec(self):
        """
        Def:
        Takes a data cube and creates a median stacked 1-D
        spectrum in a 3x3 grid around the central pixel 

        Output: 
        matplotlib plot of the 1D stacked spectrum 


        """
        lst = []
        for i in range((len(self.data[0]) / 2) - 1, (len(self.data[0]) / 2) + 1):
            for j in range((len(self.data[0]) / 2) - 1, (len(self.data[0]) / 2) + 1):
                lst.append(self.data[:,i,j])
        flux_array = np.nanmedian(lst, axis=0)
        return flux_array               

    def singlePixelExtract(self, centre_x, centre_y):
        """
        Def: 
        Extracts a 1-D spectrum at the central x and y locations provided
        Input - centre_x: the central location of the galaxy on the 2D image in the x direction 
                centre_y: the central location of the galaxy on the 2D image in the y direction
        Output - FluxArray: Extracted 1D flux spectrum from the object at the chosen location 
        """
        #Already have Data defined - want to collapse this down to a 1D array at the chosen x-y location 
        return self.data[:,centre_y, centre_x]

    def specPlot2D(self, orientation):
        """
        Def:
        Takes a data cube and creates a median stacked 2-D
        spectrum across either the horizontal row or vertical column

        Input: 
        orientation - either 'vertical' or 'horizontal' (default down)

        Output: 
        2D array specifying the 2D spectrum 


        """

        #Check the orientation input to see what has been specified#
        try:
            #If 'up' median stack across the rows
            if orientation == 'vertical':
                plot_vec = np.nanmedian(self.data, axis=1)
            elif orientation == 'horizontal':
                plot_vec = np.nanmedian(self.data, axis=2)
            else:
                raise ValueError('Choose either vertical or horizontal for the orientation')
        except ValueError:
            print 'check input for Orientation'
            raise   

        #Now have the plot vector, plot it.
        plt.imshow(plot_vec)
        plt.savefig('test.png')     
        plt.close('all')

    def optimalSpec(self):

        """
        Def: 
        Optimally extract the spectrum of the object from the whole image. 
        Use the PSF of the object to get the weighting for the extraction. 
        Do I just sum over the axis? 
        Input: None - everything already defined
        """
        #Fit a gaussian to the fully combined science cube
        #to determine the optimal extraction profile
        print '[INFO]: Fitting the optimal spectrum for object: %s' % self.IFUName
        params, psfMask, fwhm, offList = self.psfMask()

        #Multiply the cube data by the psfMask
        modCube = psfMask * self.data

        #Recover the width of the gaussian 
        width = fwhm / 2.3548
        width = int(np.round(width))
        #Recover the central value
        x = params[2]
        y = params[1]

        #Set the upper and lower limits for optimal spectrum extraction
        x_upper = int(x + (1.5*width))
        if x_upper > len(self.data[0]):
            x_upper = len(self.data[0])
        x_lower = int(x - (1.5*width))
        if x_lower < 0:
            x_lower = 0 

        y_upper = int(y + (1.5*width))
        if y_upper > len(self.data[0]):
            y_upper = len(self.data[0])
        y_lower = int(y - (1.5*width))
        if y_lower < 0:
            y_lower = 0

        print x_lower, x_upper, y_lower, y_upper        

        #Set all values greater than 2sigma from the centre = 0 
        #Don't want to include these in the final spectrum 
        modCube[:,0:y_lower, :] = 0 
        modCube[:,y_upper:len(self.data[0]), :] = 0
        modCube[:,: , 0:x_lower] = 0
        modCube[:, :, x_upper: len(self.data[0])] = 0
        #Sum over each spatial dimension to get the spectrum 
        first_sum = np.nansum(modCube, axis=1)
        spectrum = np.nansum(first_sum, axis=1)
        return spectrum

    def optimalSpecFromProfile(self, profile, fwhm, centre_x, centre_y):

        """
        Def: 
        Optimally extract the spectrum of the object from the whole image. 
        Use the PSF of the object to get the weighting for the extraction. 
        Do I just sum over the axis? 
        Input: Profile - a specified 2D normalised profile for extraction
                fwhm - the fwhm of the tracked star
        """
        try:
            print '[INFO]: Fitting the optimal spectrum for object: %s' % self.IFUName
        except AttributeError:
            print '[INFO]: Fitting the optimal spectrum'
        #Multiply the cube data by the psfMask
        modCube = profile * self.data

        #Recover the width of the gaussian 
        width = fwhm / 2.3548
        #Recover the central value
        x = copy(centre_x)
        y = copy(centre_y)
        print '[INFO]: The central values of the Gaussian are: %s %s' % (x, y)
        print '[INFO]: And the width is: %s' % width

        #Set the upper and lower limits for optimal spectrum extraction
        x_upper = int(np.round((x + (2.0*width))))
        if x_upper > len(self.data[0]):
            x_upper = len(self.data[0])
        x_lower = int(np.round((x - (2.0*width))))
        if x_lower < 0:
            x_lower = 0 

        y_upper = int(np.round((y + (2.0*width))))
        if y_upper > len(self.data[0]):
            y_upper = len(self.data[0])
        y_lower = int(np.round((y - (2.0*width))))
        if y_lower < 0:
            y_lower = 0

        print x_lower, x_upper, y_lower, y_upper        

        #Set all values greater than 2sigma from the centre = 0 
        #Don't want to include these in the final spectrum 
        modCube[:,0:y_lower, :] = 0 
        modCube[:,y_upper:len(self.data[0]), :] = 0
        modCube[:,: , 0:x_lower] = 0
        modCube[:, :, x_upper: len(self.data[0])] = 0   

        imModCube = copy(modCube)
        imModCube = np.nanmedian(modCube, axis=0)
        #Check to see that the gaussian and shifted profile align
        colFig, colAx = plt.subplots(1,1, figsize=(12.0,12.0))
        colCax = colAx.imshow(imModCube, interpolation='bicubic')
        colFig.colorbar(colCax)
        plt.show()
        plt.close('all')
        #Sum over each spatial dimension to get the spectrum 
        first_sum = np.nansum(modCube, axis=1)
        spectrum = np.nansum(first_sum, axis=1)
        return spectrum


    #Create a gaussian function for use with lmfit  
    def gaussian(self, x1, x2, pedastal, height, center_x, center_y, width_x, width_y):
        #make sure we have floating point values 
        width_x = float(width_x)
        width_y = float(width_y)
        #Specify the gaussian function here 
        func = pedastal + height*exp(-(((center_x-x1)/width_x)**2+((center_y-x2)/width_y)**2)/2)
        return func

    #Create a gaussian function for use with the integral
    def gaussianLam(self, pedastal, height, center_x, center_y, width_x, width_y):
        """Returns a gaussian function with the given parameters"""
        width_x = float(width_x)
        width_y = float(width_y)
        return lambda x,y: pedastal + height*exp(
                    -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)


    def gauss2dMod(self):
        mod = Model(self.gaussian, independent_vars=['x1', 'x2'], param_names=\
            ['pedastal','height', 'center_x', 'center_y', 'width_x', 'width_y'], missing='drop')
        #print mod.independent_vars
        #print mod.param_names
        
        return mod

    def moments(self, data):
        """Returns (pedastal, height, center_x, center_y, width_x, width_y)
        the gaussian parameters of a 2D distribution by calculating its
        moments """
        #First set all np.nan values in data to 0 
        #And all the negative values to 0 
        #These shouldn't influence the moment calculation
        #data[np.isnan(data)] = 0
        #data[data < 0] = 0
        #total = np.nansum(data)
        #print 'The sum over the data is: %s' % total
        #X, Y = indices(data.shape)
        #print 'The Indices are: %s, %s' % (X, Y)
        #x = np.nansum((X*data))/total
        #y = np.nansum((Y*data))/total
        #print x, y
        #col = data[:, int(y)]
        #if col.sum() == 0:
        #   width_x = sqrt(abs((arange(col.size)-y)**2*1.0).sum()/1.0)
        #else:
        #   width_x = sqrt(abs((arange(col.size)-y)**2*col).sum()/col.sum())
        #row = data[int(x), :]

        #if row.sum() == 0:
        #   width_y = sqrt(abs((arange(row.size)-x)**2*1.0).sum()/1.0)
        #else:
        #   width_y = sqrt(abs((arange(row.size)-x)**2*row).sum()/row.sum())
        #height = np.nanmax(data)
        #print '[INFO]: The Initial Guess at the G.Params;\nArea:%s\nCentre_y:%s\nCentre_x:%s\nWidth_y:%s\nWidth_x:%s' % (height, x, y, width_x, width_y)
        #return height, x, y, width_x, width_y, pedastal

        pedastal = np.nanmedian(data)
        height = np.nanmax(data)
        data[np.isnan(data)] = 0
        data[data < 0] = 0
        total = np.nansum(data)
        #print 'The sum over the data is: %s' % total
        X, Y = indices(data.shape)
        #print 'The Indices are: %s, %s' % (X, Y)
        center_x = np.nansum((X*data))/total
        center_y = np.nansum((Y*data))/total
        width_x = 1.2
        width_y = 1.2
        #print '[INFO]: The Initial Guess at the G.Params;\nArea:%s\nCentre_y:%s\nCentre_x:%s\nWidth_y:%s\nWidth_x:%s' % (height, x, y, width_x, width_y)
        return [height, center_x, center_y, width_x, width_y, pedastal] 

    def fitgaussian(self, data):
        """Returns (height, x, y, width_x, width_y)
        the gaussian parameters of a 2D distribution found by a fit"""
        #params = self.moments(data)
        #print '[INFO]: The Data has shape: %s %s' % (data.shape[0], data.shape[1])
        #errorfunction = lambda p: ravel(self.gaussian(*p)(*indices(data.shape)) - data)
        #p, success = optimize.leastsq(errorfunction, params)
        #print '[INFO]: The Gaussian fitting parameters are;\nArea:%s\nCentre_y:%s\nCentre_x:%s\nWidth_y:%s\nWidth_x:%s' % (p[0], p[1], p[2], p[3], p[4])
        #return p
        #At the moment will assume that the data is imagedata which needs flattened 
        pars = self.moments(data)
        flat_data = np.ndarray.flatten(data)
        print 'This is the flattened data: %s' % flat_data
        #Get the gaussian model 
        mod = self.gauss2dMod()
        #Set the parameter hints from the initialPars method 
        mod.set_param_hint('height', value=pars[0])
        mod.set_param_hint('center_x', value=pars[1])
        mod.set_param_hint('center_y', value=pars[2])
        mod.set_param_hint('width_x', value=pars[3])
        mod.set_param_hint('width_y', value=pars[4])
        mod.set_param_hint('pedastal', value=pars[5])
        #Initialise a parameters object to use in the fit
        fit_pars = mod.make_params()
        #Guess isn't implemented for this model 
        #Need to pass independent variables for the fit. these come from 
        #flattening the indices of data.shape
        x1 = np.ndarray.flatten(indices(data.shape)[0])
        #print 'The first independent variable: %s %s' % (x1, type(x1))
        x2 = np.ndarray.flatten(indices(data.shape)[1])
        #print 'The second independent variable: %s' % x2
        mod_fit = mod.fit(flat_data, x1=x1, x2=x2, params=fit_pars)
        #print mod_fit.best_values['pedastal']
        return mod_fit, x1, x2

    def psfMask(self):
        """Returns (FWHM, psfMask) which are the FWHM of the 2D gaussian 
        fit to the collapsed object image and the mask of values found 
        after integrating the gaussian function over all the pixels and 
        normalising by this value.  
        """
        #Find the FWHM and the masking profile of a given datacube

        #Step 1 - perform least squares minimisation to find the parameters  
        mod_fit, x1, x2 = self.fitgaussian(self.imData)
#       #Check to find sigma 
#       if (np.isnan(params[3]) and np.isnan(params[4])):
#           sigma = 3.0
#       elif (np.isnan(params[3]) and not np.isnan(params[4])):
#           sigma = 3.0
#       elif (not np.isnan(params[3]) and np.isnan(params[4])):
#           sigma = 3.0
#       else:
#           sigma = 0.5*(params[3] + params[4])
        #Set the params variable as the best fit attribute 
        params = mod_fit.best_values
        sigma = 0.5*(params['width_x'] + params['width_y'])
        FWHM = 2.3548 * sigma
        try: 
            print '[INFO]: The FWHM of object %s is: %s' % (self.IFUName, FWHM)
        except AttributeError:
            print '[INFO]: The FWHM is: %s' % FWHM


        #Create an instance of the gaussian for integrating 
        fit = self.gaussianLam(params['pedastal'], params['height'],\
         params['center_x'], params['center_y'], params['width_x'], params['width_y'])
        #Evaluate the gaussian with the best fit parameters  
        mod_eval = mod_fit.eval(x1=x1, x2=x2)
        #Step 3 - reshape back to 2D array 
        gEval = np.reshape(mod_eval, self.imData.shape)
        #This is the initial grid of values, but need to normalise to 1 
        #Integrate the gaussian using double quadrature 
        integral = scipy.integrate.dblquad(fit, a=0, b=self.imData.shape[1],\
            gfun=lambda x: 0 , hfun=lambda x: self.imData.shape[1])
        #Plot the image and the fit 
        colFig, colAx = plt.subplots(1,1, figsize=(14.0,14.0))
        colCax = colAx.imshow(self.imData, interpolation='bicubic')
        colAx.contour(gEval)
        colFig.colorbar(colCax)
        saveName = (self.fileName)[:-5] + '_gauss.png'
        colFig.savefig(saveName)
        #plt.show()     
        plt.close('all')
        #return the FWHM and the masked profile 
        return params, (gEval / integral[0]), FWHM, self.offList


    def plot_HK_sn_map(self, redshift, savefig=False):
        """
        Def: 
        Check the signal to noise of the emission lines over the face of a cube 
        with known redshift 
        Input: redshift - redshift of the galaxy in the cube 
               savefig - whether or not to save the figures 
        """

        fig, axes = plt.subplots(figsize=(14, 4), nrows=1, ncols=3)
        fig.subplots_adjust(right=0.83)

        cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])

        # open the data
        data = self.data
        noise = self.Table[2].data

        # get the wavelegnth index of the oiii5007 line:
        wl_0 = self.Table[1].header['CRVAL3']
        dwl = self.Table[1].header['CDELT3']
        wl_n = wl_0 + (data.shape[0] * dwl)

        wl = np.linspace(wl_0, wl_n, data.shape[0])

        # create a sn dictionary to house the line sn maps
        sn_dict = {}

        for line, ax in zip(['[OII]', 'Hb', '[OIII]5007'], axes.flatten()):

            ax.minorticks_on()

            if line == '[OIII]5007':
                oiii5007_wl = 0.500824 * (1. + redshift)
                line_idx = np.argmin(np.abs(wl - oiii5007_wl))
            elif line == 'Hb':
                hb_wl = 0.486268 * (1. + redshift)
                line_idx = np.argmin(np.abs(wl - hb_wl))
            elif line == '[OII]':
                oii_wl = 0.3729875 * (1. + redshift)
                line_idx = np.argmin(np.abs(wl - oii_wl))

            # the shape of the data is (spectrum, xpixel, ypixel)
            # loop through each x and y pixel and get the OIII5007 S/N
            xpixs = data.shape[1]
            ypixs = data.shape[2]

            sn_array = np.empty(shape=(xpixs, ypixs))

            for i, xpix in enumerate(np.arange(0, xpixs, 1)):

                for j, ypix in enumerate(np.arange(0, ypixs, 1)):

                    spaxel_spec = data[:, i, j]
                    spaxel_noise = noise[:, i, j]

                    line_counts = np.median(spaxel_spec[line_idx - 3:
                                                        line_idx + 3])

                    line_noise = np.median(spaxel_noise[line_idx - 3:
                                                        line_idx + 3])

                    line_sn = line_counts / line_noise

                    if np.isnan(line_sn):
                        sn_array[i, j] = -99.
                    else:
                        sn_array[i, j] = line_sn

            # print max(sn_array.flatten())
            #add the result to the sn_dict
            sn_dict[line] = sn_array

            im = ax.imshow(sn_array, aspect='auto', vmin=0.,
                           vmax=3.,
                           cmap=plt.get_cmap('hot'))

            ax.set_title('%s' % line)

        fig.colorbar(im, cax=cbar_ax)

        #plt.tight_layout()
        # plt.show()

        if savefig:
            fig.savefig('%s_sn_map.pdf' % self.fileName[:-5])
        # return the dictionary containing the noise values
        return sn_dict


    def plot_K_sn_map(self, redshift, savefig=False):

        fig, axes = plt.subplots(figsize=(10, 4), nrows=1, ncols=2)
        fig.subplots_adjust(right=0.83)

        cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])

        # open the data
        data = self.data
        noise = self.Table[2].data

        # get the wavelegnth index of the oiii5007 line:
        wl_0 = self.Table[1].header['CRVAL3']
        dwl = self.Table[1].header['CDELT3']
        wl_n = wl_0 + (data.shape[0] * dwl)

        wl = np.linspace(wl_0, wl_n, data.shape[0])
        # Create a dictionary to house the sn_arrays
        sn_dict = {}

        for line, ax in zip(['Hb', '[OIII]5007'], axes.flatten()):

            ax.minorticks_on()

            if line == '[OIII]5007':
                oiii5007_wl = 0.500824 * (1. + redshift)
                line_idx = np.argmin(np.abs(wl - oiii5007_wl))
            elif line == 'Hb':
                hb_wl = 0.486268 * (1. + redshift)
                line_idx = np.argmin(np.abs(wl - hb_wl))

            # the shape of the data is (spectrum, xpixel, ypixel)
            # loop through each x and y pixel and get the OIII5007 S/N
            xpixs = data.shape[1]
            ypixs = data.shape[2]

            sn_array = np.empty(shape=(xpixs, ypixs))

            for i, xpix in enumerate(np.arange(0, xpixs, 1)):

                for j, ypix in enumerate(np.arange(0, ypixs, 1)):

                    spaxel_spec = data[:, i, j]
                    spaxel_noise = noise[:, i, j]

                    line_counts = np.median(spaxel_spec[line_idx - 3:
                                                        line_idx + 3])

                    line_noise = np.median(spaxel_noise[line_idx - 3:
                                                        line_idx + 3])

                    line_sn = line_counts / line_noise

                    if np.isnan(line_sn):
                        sn_array[i, j] = -99.
                    else:
                        sn_array[i, j] = line_sn

            # print max(sn_array.flatten())
            # add the result to the sn dictionary
            sn_dict[line] = sn_array

            im = ax.imshow(sn_array, aspect='auto', vmin=0.,
                           vmax=3.,
                           cmap=plt.get_cmap('hot'))

            ax.set_title('%s' % line)

        fig.colorbar(im, cax=cbar_ax)

        #plt.tight_layout()
        # plt.show()

        if savefig:
            fig.savefig('%s_sn_map.pdf' % self.fileName[:-5])
        return sn_dict


    def plot_HK_image(self, redshift, savefig=False):
        """
        Def: 
        Check the signal to noise of the emission lines over the face of a cube 
        with known redshift 
        Input: redshift - redshift of the galaxy in the cube 
               savefig - whether or not to save the figures 
        """

        fig, axes = plt.subplots(figsize=(14, 4), nrows=1, ncols=3)
        fig.subplots_adjust(right=0.83)

        cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])

        # open the data
        data = self.data
        noise = self.Table[2].data

        # get the wavelegnth index of the oiii5007 line:
        wl_0 = self.Table[1].header['CRVAL3']
        dwl = self.Table[1].header['CDELT3']
        wl_n = wl_0 + (data.shape[0] * dwl)

        wl = np.linspace(wl_0, wl_n, data.shape[0])

        for line, ax in zip(['[OII]', 'Hb', '[OIII]5007'], axes.flatten()):

            ax.minorticks_on()

            # the shape of the data is (spectrum, xpixel, ypixel)
            # loop through each x and y pixel and get the OIII5007 S/N
            xpixs = data.shape[1]
            ypixs = data.shape[2]

            if line == '[OIII]5007':
                oiii5007_wl = 0.500824 * (1. + redshift)
                line_idx = np.argmin(np.abs(wl - oiii5007_wl))
                met_array_OIII = np.empty(shape=(xpixs, ypixs))  
            elif line == 'Hb':
                hb_wl = 0.486268 * (1. + redshift)
                line_idx = np.argmin(np.abs(wl - hb_wl))
                met_array_Hb = np.empty(shape=(xpixs, ypixs))
            elif line == '[OII]':
                oii_wl = 0.3729875 * (1. + redshift)
                line_idx = np.argmin(np.abs(wl - oii_wl))
                met_array_OII = np.empty(shape=(xpixs, ypixs))



            sn_array = np.empty(shape=(xpixs, ypixs))

            for i, xpix in enumerate(np.arange(0, xpixs, 1)):

                for j, ypix in enumerate(np.arange(0, ypixs, 1)):

                    spaxel_spec = data[:, i, j]
                    spaxel_noise = noise[:, i, j]

                    line_counts = np.median(spaxel_spec[line_idx - 3:
                                                        line_idx + 3])

                    line_noise = np.median(spaxel_noise[line_idx - 3:
                                                        line_idx + 3])

                    line_sn = line_counts / line_noise

                    if np.isnan(line_sn):
                        sn_array[i, j] = -99
                    else:
                        sn_array[i, j] = line_counts

                    if line == '[OIII]5007':
                        if line_sn < 0.8:
                            met_array_OIII[i, j] = np.nan
                        else:
                            met_array_OIII[i, j] = line_counts

                    if line == 'Hb':
                        if line_sn < 0.8:
                            met_array_Hb[i, j] = np.nan
                        else:
                            met_array_Hb[i, j] = line_counts

                    if line == '[OII]':
                        if line_sn < 0.8:
                            met_array_OII[i, j] = np.nan
                        else:
                            met_array_OII[i, j] = line_counts

            # print max(sn_array.flatten())



            im = ax.imshow(sn_array, aspect='auto', vmin=0.,
                           vmax=3.,
                           cmap=plt.get_cmap('hot'))

            ax.set_title('%s' % line)

        fig.colorbar(im, cax=cbar_ax)

        # plt.tight_layout()
        # plt.show()
        if savefig:
            fig.savefig('%s_images.pdf' % self.fileName[:-5])
        plt.close('all')    

#        # before creating the three plots create plots for each graph 
#        fig, ax = plt.subplots(1, figsize=(10, 10))
#        im = ax.imshow(met_array_OIII, aspect='auto', vmin=0.,
#                       vmax=3.,
#                       cmap=plt.get_cmap('hot'))
#        ax.set_title('[OIII]')
#        fig.colorbar(im)#

#        plt.show()
#        plt.close('all')#

#        fig, ax = plt.subplots(1, figsize=(10, 10))
#        im = ax.imshow(met_array_Hb, aspect='auto', vmin=0.,
#                       vmax=3.,
#                       cmap=plt.get_cmap('hot'))
#        ax.set_title('Hb')
#        fig.colorbar(im)#

#        plt.show()
#        plt.close('all')#

#        fig, ax = plt.subplots(1, figsize=(10, 10))
#        im = ax.imshow(met_array_OII, aspect='auto', vmin=-3.0,
#                       vmax=3.,
#                       cmap=plt.get_cmap('hot'))
#        ax.set_title('[OII]')
#        fig.colorbar(im)#

#        plt.show()
#        plt.close('all')    

        # now should also have the Hb and OIII metallicity maps
        # divide the two and plot the result 
        overall_met = met_array_OIII / met_array_Hb
        overall_met_OII = met_array_OIII / met_array_OII

        # now for each of these convert to metallicity using the Maiolino 
        # relations. The problem here is with which root of the polynomial 
        # to take. Different roots should be applicable in the high and low 
        # metallicity intervals 
        # First the Hb ratio, set up a new array to house the results 

        x_shape = overall_met.shape[0]
        y_shape = overall_met.shape[1]

        Hb_met_array = np.empty(shape=(x_shape, y_shape))


        # initialise the coefficients, given in Maiolino 2008 
        c_0_Hb = 0.1549
        c_1_Hb = -1.5031
        c_2_Hb = -0.9790
        c_3_Hb = -0.0297

        for i, xpix in enumerate(np.arange(0, x_shape, 1)):

            for j, ypix in enumerate(np.arange(0, y_shape, 1)):
                # print 'This is the number: %s' % overall_met[i, j]

                # if the number is nan, leave it as nan 

                if np.isnan(overall_met[i, j]) \
                   or np.isinf(overall_met[i, j]) \
                   or (overall_met[i, j]) < 0:

                    Hb_met_array[i, j] = np.nan

                # else subtract the log10(number) from
                # c_0_Hb and set up the polynomial from poly1D

                else:

                    c_0_Hb_new = c_0_Hb - np.log10(overall_met[i, j])

                    p = poly1d([c_3_Hb, c_2_Hb, c_1_Hb, c_0_Hb_new])
                    # print p.r
                    # the roots of the polynomial are given in units 
                    # of metallicity relative to solar. add 8.69 
                    # met_value = p.r[0] + 8.69
                    # if the root has an imaginary component, just take 
                    # the real part 
                    Hb_met_array[i, j] = p.r[2].real + 8.69

        # Next the OII ratio, set up a new array to house the results 

        x_shape = overall_met_OII.shape[0]
        y_shape = overall_met_OII.shape[1]

        OII_met_array = np.empty(shape=(x_shape, y_shape))

        # initialise the coefficients, given in Maiolino 2008 
        c_0_OII = -0.2839
        c_1_OII = -1.3881
        c_2_OII = -0.3172

        for i, xpix in enumerate(np.arange(0, x_shape, 1)):

            for j, ypix in enumerate(np.arange(0, y_shape, 1)):

                # if the number is nan, leave it as nan 
                if np.isnan(overall_met_OII[i, j]) \
                   or np.isinf(overall_met_OII[i, j]) \
                   or (overall_met_OII[i, j]) < 0:

                    OII_met_array[i, j] = np.nan

                # else subtract the number from
                # c_0_OII and set up the polynomial from poly1D

                else:
                    # print 'This is the number: %s' % overall_met_OII[i, j]
                    c_0_OII_new = c_0_OII - np.log10(overall_met_OII[i, j])

                    p = poly1d([c_2_OII, c_1_OII, c_0_OII_new])
                    # print p.r
                    # the roots of the polynomial are given in units 
                    # of metallicity relative to solar. add 8.69 
                    # met_value = p.r[0] + 8.69
                    # if the root has an imaginary component, just take 
                    # the real part 
                    if np.isreal(p.r[1]):
                        OII_met_array[i, j] = p.r[1] + 8.69
                    else:
                        OII_met_array[i, j] = -100

        fig, ax = plt.subplots(1, figsize=(10, 10))

        im = ax.imshow(Hb_met_array, aspect='auto', 
                       vmin=7.5, vmax=9.0, cmap=plt.get_cmap('jet'))

        ax.set_title('[OIII] / Hb')

        fig.colorbar(im) 

        # plt.show()
        if savefig:
            fig.savefig('%s_OIII_Hb.pdf' % self.fileName[:-5])
        plt.close('all')

        fig, ax = plt.subplots(1, figsize=(10, 10))

        im = ax.imshow(OII_met_array, aspect='auto', 
                       vmin=7.5, vmax=9.0, cmap=plt.get_cmap('jet'))

        ax.set_title('[OIII] / [OII]')

        fig.colorbar(im) 
        # plt.show()
        if savefig:
            fig.savefig('%s_OIII_OII.pdf' % self.fileName[:-5])
        plt.close('all')
        return Hb_met_array, OII_met_array


    def plot_K_image(self, redshift, savefig=False):

        fig, axes = plt.subplots(figsize=(10, 4), nrows=1, ncols=2)
        fig.subplots_adjust(right=0.83)

        cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])

        # open the data
        data = self.data
        noise = self.Table[2].data

        # get the wavelegnth index of the oiii5007 line:
        wl_0 = self.Table[1].header['CRVAL3']
        dwl = self.Table[1].header['CDELT3']
        wl_n = wl_0 + (data.shape[0] * dwl)

        wl = np.linspace(wl_0, wl_n, data.shape[0])

        for line, ax in zip(['Hb', '[OIII]5007'], axes.flatten()):

            ax.minorticks_on()

            # the shape of the data is (spectrum, xpixel, ypixel)
            # loop through each x and y pixel and get the OIII5007 S/N
            xpixs = data.shape[1]
            ypixs = data.shape[2]

            if line == '[OIII]5007':
                oiii5007_wl = 0.500824 * (1. + redshift)
                line_idx = np.argmin(np.abs(wl - oiii5007_wl))
                met_array_OIII = np.empty(shape=(xpixs, ypixs))  
            elif line == 'Hb':
                hb_wl = 0.486268 * (1. + redshift)
                line_idx = np.argmin(np.abs(wl - hb_wl))
                met_array_Hb = np.empty(shape=(xpixs, ypixs))


            sn_array = np.empty(shape=(xpixs, ypixs))

            for i, xpix in enumerate(np.arange(0, xpixs, 1)):

                for j, ypix in enumerate(np.arange(0, ypixs, 1)):

                    spaxel_spec = data[:, i, j]
                    spaxel_noise = noise[:, i, j]

                    line_counts = np.median(spaxel_spec[line_idx - 3:
                                                        line_idx + 3])

                    line_noise = np.median(spaxel_noise[line_idx - 3:
                                                        line_idx + 3])

                    line_sn = line_counts / line_noise

                    if np.isnan(line_sn):
                        sn_array[i, j] = -99
                    else:
                        sn_array[i, j] = line_counts

                    if line == '[OIII]5007':
                        if line_sn < 1.0:
                            met_array_OIII[i, j] = np.nan
                        else:
                            met_array_OIII[i, j] = line_counts

                    if line == 'Hb':
                        if line_sn < 1.0:
                            met_array_Hb[i, j] = np.nan
                        else:
                            met_array_Hb[i, j] = line_counts

            # print max(sn_array.flatten())

            im = ax.imshow(sn_array, aspect='auto', vmin=0.,
                           vmax=3.,
                           cmap=plt.get_cmap('hot'))

            ax.set_title('%s' % line)

        fig.colorbar(im, cax=cbar_ax)

        # plt.tight_layout()
        # plt.show()
        if savefig:
            fig.savefig('%s_images.pdf' % self.fileName[:-5])
        plt.close('all')        

        # now should also have the Hb and OIII metallicity maps
        # divide the two and plot the result 
        overall_met = met_array_OIII / met_array_Hb

        # now for each of these convert to metallicity using the Maiolino 
        # relations. The problem here is with which root of the polynomial 
        # to take. Different roots should be applicable in the high and low 
        # metallicity intervals 
        # First the Hb ratio, set up a new array to house the results 

        x_shape = overall_met.shape[0]
        y_shape = overall_met.shape[1]

        Hb_met_array = np.empty(shape=(x_shape, y_shape))


        # initialise the coefficients, given in Maiolino 2008 
        c_0_Hb = 0.1549
        c_1_Hb = -1.5031
        c_2_Hb = -0.9790
        c_3_Hb = -0.0297

        for i, xpix in enumerate(np.arange(0, x_shape, 1)):

            for j, ypix in enumerate(np.arange(0, y_shape, 1)):
                # print 'This is the number: %s' % overall_met[i, j]

                # if the number is nan, leave it as nan 

                if np.isnan(overall_met[i, j]) \
                   or np.isinf(overall_met[i, j]) \
                   or (overall_met[i, j]) < 0:

                    Hb_met_array[i, j] = np.nan

                # else subtract the log10(number) from
                # c_0_Hb and set up the polynomial from poly1D

                else:

                    c_0_Hb_new = c_0_Hb - np.log10(overall_met[i, j])

                    p = poly1d([c_3_Hb, c_2_Hb, c_1_Hb, c_0_Hb_new])
                    # print p.r
                    # the roots of the polynomial are given in units 
                    # of metallicity relative to solar. add 8.69 
                    # met_value = p.r[0] + 8.69
                    # if the root has an imaginary component, just take 
                    # the real part 
                    Hb_met_array[i, j] = p.r[2].real + 8.69

        fig, ax = plt.subplots(1, figsize=(14, 14))

        im = ax.imshow(Hb_met_array, aspect='auto', 
                       vmin=7.5, vmax=9.0, cmap=plt.get_cmap('jet'))

        ax.set_title('log([OIII] / Hb)')

        fig.colorbar(im) 
        # plt.show()
        if savefig:
            fig.savefig('%s_OIII_Hb.pdf' % self.fileName[:-5])
        plt.close('all')
        return Hb_met_array

    def spaxel_binning(self, data, xbin, ybin, interp='median'):

        """
        Def: bins spaxels in xbin and ybin shaped chunks. Sticking with the 
        convention that xbin will always refer to the first index.
        Input: xbin - number of spaxels to go into 1 along x direction
               ybin - number of spaxels to go into 1 along y direction
        """
        # the data is a 3D cube - need to preserve this 
        # use a for loop with a step size 

        # initially set the dimensions of the data cube
        xlength = data.shape[1]
        ylength = data.shape[2]

        print 'The original cube dimensions are: (%s, %s)' % (xlength, ylength)

        # calculate what the shape of the final array will be 
        # this is tricky if the bin sizes do not match the shape of 
        # the array. In this case take the modulus result and use that 
        # as the final bin width (more often than not this will be 1)
        # for the shape of the final array this means that it is given 
        # by the initial shape / binsize + 1 (if % != 0)

        if xlength % xbin == 0:
            new_xlength = xlength / xbin

        else:
            # account for the additional bin at the edge
            new_xlength = (xlength / xbin) + 1

        if ylength % ybin == 0:
            new_ylength = ylength / ybin

        else:
            # account for the additional bin at the edge
            new_ylength = (ylength / ybin) + 1

        # create the new array
        new_data = np.empty(shape=(data.shape[0], new_xlength, new_ylength))

        # loop round and create the new spaxels 
        # during each loop component need to check 
        # if the indices are reaching the original data size 
        # and if so create the final bin using modulus 
        # then save each 1D array at the appropriate location 
        # in the new_data cube

        # set counters to record the position to store the new spaxel 
        # in the new_data array 
        x_cube_counter = 0
        
        for x in range(0, xlength, xbin):
            y_cube_counter = 0
            # check if the xlength has been reached or exceeded 
            if x + xbin >= xlength:

                # need a different y for loop that uses the end point 
                # limits in the x-direction. First find what these are 
                modulus_x = xlength % xbin
                start_x = (xlength - modulus_x) + 1

                # initiate for loop for this scenario
                for y in range(0, ylength, ybin):

                    # this configuration means we will first 
                    # be looping down the way, for each row 
                    if y + ybin >= ylength:

                        # we've exceeded the original spaxel limit 
                        # meaning that indicing will fail. create the final bin
                        modulus_y = ylength % ybin
                        start_y = (ylength - modulus_y) + 1

                        # note the + 1 is required for proper indexing

                        # now take into account the chosen interpolation type
                        if interp == 'sum':
                            # print 'Sum interpolation chosen' 
                            new_spaxel = np.nansum(\
                                         np.nansum(data[:, start_x:xlength - 1,\
                                         start_y:ylength - 1], axis=1), axis=1)

                        elif interp == 'mean':
                            # print 'Mean interpolation chosen' 
                            new_spaxel = np.nanmean(\
                                         np.nanmean(data[:, start_x:xlength - 1,\
                                         start_y:ylength - 1], axis=1), axis=1)

                        # default value of median
                        else: 
                            new_spaxel = np.nanmedian(\
                                         np.nanmedian(data[:, start_x:xlength - 1,\
                                         start_y:ylength - 1], axis=1), axis=1)

                    # everything is okay, limit not exceeded
                    else:

                        if interp == 'sum':
                            # print 'Sum interpolation chosen' 
                            new_spaxel = np.nansum(\
                                         np.nansum(data[:, start_x:xlength - 1,\
                                         y:y + ybin], axis=1), axis=1)

                        elif interp == 'mean':
                            # print 'Mean interpolation chosen'
                            new_spaxel = np.nanmean(\
                                         np.nanmean(data[:, start_x:xlength - 1,\
                                         y:y + ybin], axis=1), axis=1)
                                                                                             
                        # default value of median
                        else: 
                            new_spaxel = np.nanmedian(\
                                         np.nanmedian(data[:, start_x:xlength - 1,\
                                         y:y + ybin], axis=1), axis=1)
                        
                    # add the new spaxel to the new_data
                    # cube in the correct position
                    new_data[:, x_cube_counter, y_cube_counter] = new_spaxel

                    # increment both the x and y cube counters
                    y_cube_counter += 1

            else:

                # everything is okay and the xlimit has not been reached
                for y in range(0, ylength, ybin):

                    # this configuration means we will first 
                    # be looping down the way, for each row 
                    if y + ybin >= ylength:

                        # we've exceeded the original spaxel limit 
                        # meaning that indicing will fail. create the final bin
                        modulus_y = ylength % ybin
                        start_y = (ylength - modulus_y) + 1

                        if interp == 'sum':
                            # print 'Sum interpolation chosen' 
                            new_spaxel = np.nansum(\
                                         np.nansum(data[:, x:x + xbin,\
                                         start_y:ylength - 1], axis=1), axis=1)

                        elif interp == 'mean':
                            # print 'Mean interpolation chosen' 
                            new_spaxel = np.nanmean(\
                                         np.nanmean(data[:, x:x + xbin,\
                                         start_y:ylength - 1], axis=1), axis=1)
                                                                                             
                        # default value of median
                        else: 
                            new_spaxel = np.nanmedian(\
                                         np.nanmedian(data[:, x:x + xbin,\
                                         start_y:ylength - 1], axis=1), axis=1)

                    # everything is okay, limit not exceeded
                    else:

                        if interp == 'sum': 
                            # print 'Sum interpolation chosen'
                            new_spaxel = np.nansum(\
                                         np.nansum(data[:, x:x + xbin,\
                                         y:y + ybin], axis=1), axis=1)

                        elif interp == 'mean':
                            # print 'Mean interpolation chosen' 
                            new_spaxel = np.nanmean(\
                                         np.nanmean(data[:, x:x + xbin,\
                                         y:y + ybin], axis=1), axis=1)
                                                                                             
                        # default value of median
                        else: 
                            new_spaxel = np.nanmedian(\
                                         np.nanmedian(data[:, x:x + xbin,\
                                         y:y + ybin], axis=1), axis=1)

                    # add the new spaxel to the new_data 
                    # cube in the correct position
                    new_data[:, x_cube_counter, y_cube_counter] = new_spaxel

                    # increment both the x and y cube counters
                    y_cube_counter += 1
            x_cube_counter += 1

        # return the new_data 
        return new_data


    def OIII_vel_map(self, redshift, savefig=False, binning=False, **kwargs):

        #TODO: CHECK THAT THE kwargs values xbin and ybin are ints < 10
        # open the data
        data = self.data
        noise = self.Table[2].data

        # get the wavelegnth index of the oiii5007 line:
        wl_0 = self.Table[1].header['CRVAL3']
        dwl = self.Table[1].header['CDELT3']
        wl_n = wl_0 + (data.shape[0] * dwl)

        wl = np.linspace(wl_0, wl_n, data.shape[0])

        # if binning is true, take the median of adjacent spaxels
        # this uses the spaxel_binning method which can bin in any 
        # different combination of shapes 
        if binning:

            # redefine the data as the binned spaxel equivalent
            if (not kwargs['xbin'] or not kwargs['ybin']):

                # haven't specified the bins 
                raise ValueError("Missing keyword arguments for binsize")

            # have specified the bins - take the values 
            else:

                # check that both bins are integers less than 10

                if (np.equal(np.mod(kwargs['xbin'], 1), 0) 
                    and kwargs['xbin'] < 10.0 
                    and np.equal(np.mod(kwargs['ybin'], 1), 0) 
                    and kwargs['ybin'] < 10.0): 

                    xbin = kwargs['xbin']
                    ybin = kwargs['ybin']

                else: 
                    raise ValueError("Non-integer or binsize too large")

            # check for an interpolation keyword 
            if kwargs['interp'] and kwargs['interp'] == 'sum':

                data = self.spaxel_binning(data, xbin, ybin, interp='sum')
                noise = self.spaxel_binning(noise, xbin, ybin, interp='sum')

            elif kwargs['interp'] and kwargs['interp'] == 'mean':

                data = self.spaxel_binning(data, xbin, ybin, interp='mean')
                noise = self.spaxel_binning(noise, xbin, ybin, interp='mean')

            # default median value chosen
            else:                                    
            # important that the data and noise have the same binning
                data = self.spaxel_binning(data, xbin, ybin)
                noise = self.spaxel_binning(noise, xbin, ybin)


        # the shape of the data is (spectrum, xpixel, ypixel)
        # loop through each x and y pixel and get the OIII5007 S/N
        xpixs = data.shape[1]
        ypixs = data.shape[2]

        # set the central wavelength of the OIII line 
        oiii5007_wl = 0.500824 * (1. + redshift)
        line_idx = np.argmin(np.abs(wl - oiii5007_wl))

        # initialise the empty velocity array
        OIII_vel_array = np.empty(shape=(xpixs, ypixs))  
        OIII_sigma_array = np.empty(shape=(xpixs, ypixs))

        for i, xpix in enumerate(np.arange(0, xpixs, 1)):

            for j, ypix in enumerate(np.arange(0, ypixs, 1)):

                spaxel_spec = data[:, i, j]
                spaxel_noise = noise[:, i, j]

                line_counts = np.median(spaxel_spec[line_idx - 3:
                                                    line_idx + 3])

                line_noise = np.median(spaxel_noise[line_idx - 3:
                                                    line_idx + 3])

                line_sn = line_counts / line_noise

                # check for nan, inf, poor s/n
                if np.isnan(line_sn):
                    OIII_vel_array[i, j] = np.nan
                    OIII_sigma_array[i, j] = np.nan

                elif np.isinf(line_sn):
                    OIII_vel_array[i, j] = np.nan
                    OIII_sigma_array[i, j] = np.nan

                elif line_sn < 1.0:
                    OIII_vel_array[i, j] = np.nan
                    OIII_sigma_array[i, j] = np.nan

                # now the condition where we have good s/n
                # can fit a gaussian to the data in each spaxel

                else:

                    # isolate the flux and wavelength data 
                    # to be used in the gaussian fit
                    # print 'Gaussian fitting spaxel [%s,%s]' % (i, j)

                    fit_wl = wl[line_idx - 6: line_idx + 6]
                    fit_flux = spaxel_spec[line_idx - 6: line_idx + 6]


                    
                    # construct gaussian model using lmfit
                    gmod = GaussianModel()
                    # set the initial parameter values 
                    pars = gmod.make_params()

                    pars['center'].set(value=oiii5007_wl, 
                                       min=oiii5007_wl - 0.0015, 
                                       max=oiii5007_wl + 0.0015)

                    pars['sigma'].set(0.0004)
                    pars['amplitude'].set(0.001)

                    # perform the fit 
                    out = gmod.fit(fit_flux, pars, x=fit_wl)

                    # assuming that the redshift measured in qfits is the 
                    # correct one - subtract the fitted centre and convert 
                    # to kms-1
                    c = 2.99792458E5
                    OIII_vel = c * ((out.best_values['center'] - oiii5007_wl) / oiii5007_wl)         
                    OIII_sig = c * ((out.best_values['sigma']) / oiii5007_wl)
                    # add this result to the velocity array 
                    OIII_vel_array[i, j] = OIII_vel
                    OIII_sigma_array[i, j] = OIII_sig

        # create a plot of the velocity field

        vel_fig, vel_ax = plt.subplots(figsize=(14, 6), nrows=1, ncols=2)
        # vel_fig.subplots_adjust(right=0.83)
        # cbar_ax = vel_fig.add_axes([0.85, 0.15, 0.02, 0.7])
        vel_ax[0].minorticks_on()
        vel_ax[1].minorticks_on()

        vel_min, vel_max = np.nanpercentile(OIII_vel_array, [2.5, 97.5])
        sig_min, sig_max = np.nanpercentile(OIII_sigma_array, [2.5, 97.5])

        im_vel = vel_ax[0].imshow(OIII_vel_array, aspect='auto', 
                           vmin=vel_min,
                           vmax=vel_max,
                           interpolation='nearest',
                           cmap=plt.get_cmap('jet'))

        vel_ax[0].set_title('[OIII] velocity')

        # add colourbar to each plot 
        divider_vel = make_axes_locatable(vel_ax[0])
        cax_vel = divider_vel.append_axes('right', size='10%', pad=0.05)
        plt.colorbar(im_vel, cax=cax_vel)

        im_sig = vel_ax[1].imshow(OIII_sigma_array, aspect='auto', 
                           vmin=sig_min,
                           vmax=sig_max,
                           interpolation='nearest',
                           cmap=plt.get_cmap('jet'))

        vel_ax[1].set_title('[OIII] Dispersion')

        # add colourbar to each plot 
        divider_sig = make_axes_locatable(vel_ax[1])
        cax_sig = divider_sig.append_axes('right', size='10%', pad=0.05)
        plt.colorbar(im_sig, cax=cax_sig)


        # vel_fig.colorbar(im)

        # plt.tight_layout()
        # plt.show()
        if savefig:
            if binning:
                vel_fig.savefig('%s_velocity_OIII_binned.pdf'\
                                % self.fileName[:-5])
            else:
                vel_fig.savefig('%s_velocity_OIII.pdf' % self.fileName[:-5])
        plt.close('all')

        # also write out the velocity array to a fits image file 
        # will use a very simple format now with no header and 
        # only a single primary extension 

        hdu = fits.PrimaryHDU(OIII_vel_array)
        hdu.writeto('%s_velocity_map.fits' % self.fileName[:-5], clobber=True)

        # return the velocity array  
        return OIII_vel_array, OIII_sigma_array


    def OII_vel_map(self, redshift, savefig=False, binning=False, **kwargs):

        #TODO: CHECK THAT THE kwargs values xbin and ybin are ints < 10
        # open the data
        data = self.data
        noise = self.Table[2].data

        # get the wavelegnth index of the oiii5007 line:
        wl_0 = self.Table[1].header['CRVAL3']
        dwl = self.Table[1].header['CDELT3']
        wl_n = wl_0 + (data.shape[0] * dwl)

        wl = np.linspace(wl_0, wl_n, data.shape[0])

        # if binning is true, take the median of adjacent spaxels
        # this uses the spaxel_binning method which can bin in any 
        # different combination of shapes 
        if binning:

            # redefine the data as the binned spaxel equivalent
            if (not kwargs['xbin'] or not kwargs['ybin']):

                # haven't specified the bins 
                raise ValueError("Missing keyword arguments for binsize")

            # have specified the bins - take the values 
            else:

                # check that both bins are integers less than 10

                if (np.equal(np.mod(kwargs['xbin'], 1), 0) 
                    and kwargs['xbin'] < 10.0 
                    and np.equal(np.mod(kwargs['ybin'], 1), 0) 
                    and kwargs['ybin'] < 10.0): 

                    xbin = kwargs['xbin']
                    ybin = kwargs['ybin']

                else: 
                    raise ValueError("Non-integer or binsize too large")

            # check for an interpolation keyword 
            if kwargs['interp'] and kwargs['interp'] == 'sum':

                data = self.spaxel_binning(data, xbin, ybin, interp='sum')
                noise = self.spaxel_binning(noise, xbin, ybin, interp='sum')

            elif kwargs['interp'] and kwargs['interp'] == 'mean':

                data = self.spaxel_binning(data, xbin, ybin, interp='mean')
                noise = self.spaxel_binning(noise, xbin, ybin, interp='mean')

            # default median value chosen
            else:                                    
            # important that the data and noise have the same binning
                data = self.spaxel_binning(data, xbin, ybin)
                noise = self.spaxel_binning(noise, xbin, ybin)


        # the shape of the data is (spectrum, xpixel, ypixel)
        # loop through each x and y pixel and get the OIII5007 S/N
        xpixs = data.shape[1]
        ypixs = data.shape[2]

        # set the central wavelength of the OIII line 
        oii_wl = 0.3729875 * (1. + redshift)
        line_idx = np.argmin(np.abs(wl - oii_wl))

        # initialise the empty velocity array
        OII_vel_array = np.empty(shape=(xpixs, ypixs))  
        OII_sigma_array = np.empty(shape=(xpixs, ypixs))

        for i, xpix in enumerate(np.arange(0, xpixs, 1)):

            for j, ypix in enumerate(np.arange(0, ypixs, 1)):

                spaxel_spec = data[:, i, j]
                spaxel_noise = noise[:, i, j]

                line_counts = np.median(spaxel_spec[line_idx - 3:
                                                    line_idx + 3])

                line_noise = np.median(spaxel_noise[line_idx - 3:
                                                    line_idx + 3])

                line_sn = line_counts / line_noise

                # check for nan, inf, poor s/n
                if np.isnan(line_sn):
                    OII_vel_array[i, j] = np.nan
                    OII_sigma_array[i, j] = np.nan

                elif np.isinf(line_sn):
                    OII_vel_array[i, j] = np.nan
                    OII_sigma_array[i, j] = np.nan

                elif line_sn < 1.0:
                    OII_vel_array[i, j] = np.nan
                    OII_sigma_array[i, j] = np.nan

                # now the condition where we have good s/n
                # can fit a gaussian to the data in each spaxel

                else:

                    # isolate the flux and wavelength data 
                    # to be used in the gaussian fit
                    # print 'Gaussian fitting spaxel [%s,%s]' % (i, j)

                    fit_wl = wl[line_idx - 6: line_idx + 6]
                    fit_flux = spaxel_spec[line_idx - 6: line_idx + 6]


                    
                    # construct gaussian model using lmfit
                    gmod = GaussianModel()
                    # set the initial parameter values 
                    pars = gmod.make_params()

                    pars['center'].set(value=oii_wl, 
                                       min=oii_wl - 0.0015, 
                                       max=oii_wl + 0.0015)

                    pars['sigma'].set(0.0004)
                    pars['amplitude'].set(0.001)

                    # perform the fit 
                    out = gmod.fit(fit_flux, pars, x=fit_wl)

                    # assuming that the redshift measured in qfits is the 
                    # correct one - subtract the fitted centre and convert 
                    # to kms-1
                    c = 2.99792458E5
                    OII_vel = c * ((out.best_values['center'] - oii_wl) / oii_wl)         
                    OII_sig = c * ((out.best_values['sigma']) / oii_wl)
                    # add this result to the velocity array 
                    OII_vel_array[i, j] = OII_vel
                    OII_sigma_array[i, j] = OII_sig

        # create a plot of the velocity field

        vel_fig, vel_ax = plt.subplots(figsize=(14, 6), nrows=1, ncols=2)
        # vel_fig.subplots_adjust(right=0.83)
        # cbar_ax = vel_fig.add_axes([0.85, 0.15, 0.02, 0.7])
        vel_ax[0].minorticks_on()
        vel_ax[1].minorticks_on()

        vel_min, vel_max = np.nanpercentile(OII_vel_array, [2.5, 97.5])
        sig_min, sig_max = np.nanpercentile(OII_sigma_array, [2.5, 97.5])

        im_vel = vel_ax[0].imshow(OII_vel_array, aspect='auto', 
                           vmin=vel_min,
                           vmax=vel_max,
                           interpolation='nearest',
                           cmap=plt.get_cmap('jet'))

        vel_ax[0].set_title('[OII] velocity')

        # add colourbar to each plot 
        divider_vel = make_axes_locatable(vel_ax[0])
        cax_vel = divider_vel.append_axes('right', size='10%', pad=0.05)
        plt.colorbar(im_vel, cax=cax_vel)

        im_sig = vel_ax[1].imshow(OII_sigma_array, aspect='auto', 
                           vmin=sig_min,
                           vmax=sig_max,
                           interpolation='nearest',
                           cmap=plt.get_cmap('jet'))

        vel_ax[1].set_title('[OII] Dispersion')

        # add colourbar to each plot 
        divider_sig = make_axes_locatable(vel_ax[1])
        cax_sig = divider_sig.append_axes('right', size='10%', pad=0.05)
        plt.colorbar(im_sig, cax=cax_sig)


        # vel_fig.colorbar(im)

        # plt.tight_layout()
        # plt.show()
        if savefig:
            if binning:
                vel_fig.savefig('%s_velocity_OII_binned.pdf'\
                                % self.fileName[:-5])
            else:
                vel_fig.savefig('%s_velocity_OII.pdf' % self.fileName[:-5])
        plt.close('all')

        # also write out the velocity array to a fits image file 
        # will use a very simple format now with no header and 
        # only a single primary extension 

        hdu = fits.PrimaryHDU(OII_vel_array)
        hdu.writeto('%s_velocity_map.fits' % self.fileName[:-5], clobber=True)

        # return the velocity array  
        return OII_vel_array, OII_sigma_array

    def Hb_vel_map(self, redshift, savefig=False, binning=False, **kwargs):

        #TODO: CHECK THAT THE kwargs values xbin and ybin are ints < 10
        # open the data
        data = self.data
        noise = self.Table[2].data

        # get the wavelegnth index of the oiii5007 line:
        wl_0 = self.Table[1].header['CRVAL3']
        dwl = self.Table[1].header['CDELT3']
        wl_n = wl_0 + (data.shape[0] * dwl)

        wl = np.linspace(wl_0, wl_n, data.shape[0])

        # if binning is true, take the median of adjacent spaxels
        # this uses the spaxel_binning method which can bin in any 
        # different combination of shapes 
        if binning:

            # redefine the data as the binned spaxel equivalent
            if (not kwargs['xbin'] or not kwargs['ybin']):

                # haven't specified the bins 
                raise ValueError("Missing keyword arguments for binsize")

            # have specified the bins - take the values 
            else:

                # check that both bins are integers less than 10

                if (np.equal(np.mod(kwargs['xbin'], 1), 0) 
                    and kwargs['xbin'] < 10.0 
                    and np.equal(np.mod(kwargs['ybin'], 1), 0) 
                    and kwargs['ybin'] < 10.0): 

                    xbin = kwargs['xbin']
                    ybin = kwargs['ybin']

                else: 
                    raise ValueError("Non-integer or binsize too large")

            # check for an interpolation keyword 
            if kwargs['interp'] and kwargs['interp'] == 'sum':

                data = self.spaxel_binning(data, xbin, ybin, interp='sum')
                noise = self.spaxel_binning(noise, xbin, ybin, interp='sum')

            elif kwargs['interp'] and kwargs['interp'] == 'mean':

                data = self.spaxel_binning(data, xbin, ybin, interp='mean')
                noise = self.spaxel_binning(noise, xbin, ybin, interp='mean')

            # default median value chosen
            else:                                    
            # important that the data and noise have the same binning
                data = self.spaxel_binning(data, xbin, ybin)
                noise = self.spaxel_binning(noise, xbin, ybin)


        # the shape of the data is (spectrum, xpixel, ypixel)
        # loop through each x and y pixel and get the OIII5007 S/N
        xpixs = data.shape[1]
        ypixs = data.shape[2]

        # set the central wavelength of the OIII line 
        hb_wl = 0.486268 * (1. + redshift)
        line_idx = np.argmin(np.abs(wl - hb_wl))

        # initialise the empty velocity array
        Hb_vel_array = np.empty(shape=(xpixs, ypixs))  
        Hb_sigma_array = np.empty(shape=(xpixs, ypixs))

        for i, xpix in enumerate(np.arange(0, xpixs, 1)):

            for j, ypix in enumerate(np.arange(0, ypixs, 1)):

                spaxel_spec = data[:, i, j]
                spaxel_noise = noise[:, i, j]

                line_counts = np.median(spaxel_spec[line_idx - 3:
                                                    line_idx + 3])

                line_noise = np.median(spaxel_noise[line_idx - 3:
                                                    line_idx + 3])

                line_sn = line_counts / line_noise

                # check for nan, inf, poor s/n
                if np.isnan(line_sn):
                    Hb_vel_array[i, j] = np.nan
                    Hb_sigma_array[i, j] = np.nan

                elif np.isinf(line_sn):
                    Hb_vel_array[i, j] = np.nan
                    Hb_sigma_array[i, j] = np.nan

                elif line_sn < 1.0:
                    Hb_vel_array[i, j] = np.nan
                    Hb_sigma_array[i, j] = np.nan

                # now the condition where we have good s/n
                # can fit a gaussian to the data in each spaxel

                else:

                    # isolate the flux and wavelength data 
                    # to be used in the gaussian fit
                    # print 'Gaussian fitting spaxel [%s,%s]' % (i, j)

                    fit_wl = wl[line_idx - 6: line_idx + 6]
                    fit_flux = spaxel_spec[line_idx - 6: line_idx + 6]


                    
                    # construct gaussian model using lmfit
                    gmod = GaussianModel()
                    # set the initial parameter values 
                    pars = gmod.make_params()

                    pars['center'].set(value=hb_wl, 
                                       min=hb_wl - 0.0015, 
                                       max=hb_wl + 0.0015)

                    pars['sigma'].set(0.0004)
                    pars['amplitude'].set(0.001)

                    # perform the fit 
                    out = gmod.fit(fit_flux, pars, x=fit_wl)

                    # assuming that the redshift measured in qfits is the 
                    # correct one - subtract the fitted centre and convert 
                    # to kms-1
                    c = 2.99792458E5
                    Hb_vel = c * ((out.best_values['center'] - hb_wl) / hb_wl)         
                    Hb_sig = c * ((out.best_values['sigma']) / hb_wl)
                    # add this result to the velocity array 
                    Hb_vel_array[i, j] = Hb_vel
                    Hb_sigma_array[i, j] = Hb_sig

        # create a plot of the velocity field

        vel_fig, vel_ax = plt.subplots(figsize=(14, 6), nrows=1, ncols=2)
        # vel_fig.subplots_adjust(right=0.83)
        # cbar_ax = vel_fig.add_axes([0.85, 0.15, 0.02, 0.7])
        vel_ax[0].minorticks_on()
        vel_ax[1].minorticks_on()

        vel_min, vel_max = np.nanpercentile(Hb_vel_array, [2.5, 97.5])
        sig_min, sig_max = np.nanpercentile(Hb_sigma_array, [2.5, 97.5])

        im_vel = vel_ax[0].imshow(Hb_vel_array, aspect='auto', 
                           vmin=vel_min,
                           vmax=vel_max,
                           interpolation='nearest',
                           cmap=plt.get_cmap('jet'))

        vel_ax[0].set_title('[Hb] velocity')

        # add colourbar to each plot 
        divider_vel = make_axes_locatable(vel_ax[0])
        cax_vel = divider_vel.append_axes('right', size='10%', pad=0.05)
        plt.colorbar(im_vel, cax=cax_vel)

        im_sig = vel_ax[1].imshow(Hb_sigma_array, aspect='auto', 
                           vmin=sig_min,
                           vmax=sig_max,
                           interpolation='nearest',
                           cmap=plt.get_cmap('jet'))

        vel_ax[1].set_title('[Hb] Dispersion')

        # add colourbar to each plot 
        divider_sig = make_axes_locatable(vel_ax[1])
        cax_sig = divider_sig.append_axes('right', size='10%', pad=0.05)
        plt.colorbar(im_sig, cax=cax_sig)


        # vel_fig.colorbar(im)

        # plt.tight_layout()
        # plt.show()
        if savefig:
            if binning:
                vel_fig.savefig('%s_velocity_Hb_binned.pdf'\
                                % self.fileName[:-5])
            else:
                vel_fig.savefig('%s_velocity_Hb.pdf' % self.fileName[:-5])
        plt.close('all')

        # also write out the velocity array to a fits image file 
        # will use a very simple format now with no header and 
        # only a single primary extension 

        hdu = fits.PrimaryHDU(Hb_vel_array)
        hdu.writeto('%s_velocity_map.fits' % self.fileName[:-5], clobber=True)

        # return the velocity array  
        return Hb_vel_array, Hb_sigma_array
##############################################################################
#Uncomment to create test instance of class and try out the methods###########
##############################################################################
#   def testMeth(self):
#       print self.data 

#create class instance 
#cube.specPlot2D(orientation='vertical')
##############################################################################


#data = np.genfromtxt('15names.txt', dtype='str')
#Save the names and types as lists 
#print data
#for name in data:
#   cube = cubeOps(name)
#   cube.specPlot(1)

##Attempting to fit 2D gaussian
#cube = cubeOps('sci_combined_n55_19__skytweak.fits')
#params = cube.fitgaussian(cube.imData)
#fit = cube.gaussian(*params)
#print params
##indices creates the grid of x,y pairs to evaluate the gaussian at. 
#cube.colAx.contour(fit(*indices(cube.imData.shape)))#

##print fit(*indices(cube.imData.shape))
##Now try evaluating the gaussian at each point
#integral = scipy.integrate.dblquad(fit, a=0, b=cube.imData.shape[1],\
# gfun=lambda x: 0 , hfun=lambda x: cube.imData.shape[1])
#print integral[0]
##Now just need to fold all of this into a function in pipeline class  
##to recover the FWHM of the gaussian in each case and plot this in the 
##same way as before. Each analysis gives a full set of IFU combined frames 
##Each of those gives a FWHM. 







