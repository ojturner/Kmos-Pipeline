#This class houses the methods which are relevant to performing manipulations 
#on the reconstructed and combined data cubes. As such, the class will be a datacube  
#object 


#import the relevant modules
import os, sys, numpy as np, random, matplotlib.mlab as mlab, matplotlib.pyplot as plt, math, matplotlib.cm as cm
import pyraf
import numpy.polynomial.polynomial as poly
import lmfit
import scipy
import math
from lmfit.models import GaussianModel, ExponentialModel, LorentzianModel, VoigtModel, PolynomialModel
from scipy import stats
from scipy import optimize
from scipy.optimize import minimize
from scipy.optimize import basinhopping
from astropy.io import fits
from matplotlib.backends.backend_pdf import PdfPages
from pylab import *
from matplotlib.colors import LogNorm


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
		self.imData = np.median(self.data, axis=0)
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
			print 'This is a raw image'
		#Extract the IFU number from the data 
		#Cube may not have this keyword, try statement
		#Not sure if this is good practice - conditional attribute.
		try:
			self.IFUNR = self.dataHeader["HIERARCH ESO PRO IFUNR"]
			key = 'HIERARCH ESO OCS ARM' + str(self.IFUNR) + ' NAME'
			#print key
			self.IFUName = self.dataHeader[key]
		except KeyError:
			print 'This is not a combined Frame, setting arm name...'
			try:
				self.IFUNR = self.dataHeader["HIERARCH ESO OCS ARM1 NAME"]
				self.IFUName = self.primHeader["HIERARCH ESO OCS ARM" + str(self.IFUNR) + " NAME"]
				print 'You have specified a reconstructed type'
			except KeyError:
				print("Warning: not a datacube")

		#Set the RA and DEC positions of all the arms. These are in 
		#sexagesimal format - convert to degrees for the plot 
		self.raDict = {}
		self.decDict = {}
		self.combDict = {}
		self.offList = []
		for i in range(1, 25):

			#Construct the list of combined science frames that will 
			#be thrown out by the science pipeline
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
				print 'IFU %s Not in Use' % DictName
				self.offList.append(i)

		#Construct the list of combined science names separately
		#This is now in order of the IFU
		self.combNames = []
		for entry in self.combDict.keys():
			combinedName = 'sci_combined_' + entry + '__skytweak.fits'
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

		#Find the pixel scale if this is a combined cube 
		try:
			self.pix_scale = self.primHeader['HIERARCH ESO PRO REC1 PARAM7 VALUE']
		except KeyError:
			print 'Could not set pixel scale - not a datacube'
			self.pix_scale = 0

		#Create the wavelength array if this is a combined data type
		try:
			self.wave_array = self.start_L + np.arange(0, 2048*(self.dL), self.dL)
		except:
			print 'cannot set wavelength array'
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
			print i
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
		print 'The Brightest Pixel is at: (%s, %s)' % (self.ind1, self.ind2)
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
		x_upper = int(x + (2*width))
		if x_upper > len(self.data[0]):
			x_upper = len(self.data[0])
		x_lower = int(x - (2*width))
		if x_lower < 0:
			x_lower = 0	

		y_upper = int(y + (2*width))
		if y_upper > len(self.data[0]):
			y_upper = len(self.data[0])
		y_lower = int(y - (2*width))
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

		#Multiply the cube data by the psfMask
		modCube = profile * self.data

		#Recover the width of the gaussian 
		width = fwhm / 2.3548
		#Recover the central value
		x = copy(centre_x)
		y = copy(centre_y)
		print 'The central values of the Gaussian are: %s %s' % (x, y)
		print 'And the width is: %s' % width

		#Set the upper and lower limits for optimal spectrum extraction
		x_upper = int(np.round((x + (1.5*width))))
		if x_upper > len(self.data[0]):
			x_upper = len(self.data[0])
		x_lower = int(np.round((x - (1.5*width))))
		if x_lower < 0:
			x_lower = 0	

		y_upper = int(np.round((y + (1.5*width))))
		if y_upper > len(self.data[0]):
			y_upper = len(self.data[0])
		y_lower = int(np.round((y - (1.5*width))))
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


	#Attempt to use someone elses code to fit a 2D gaussian to the data	
	def gaussian(self, height, center_x, center_y, width_x, width_y):
	    """Returns a gaussian function with the given parameters"""
	    width_x = float(width_x)
	    width_y = float(width_y)
	    return lambda x,y: height*exp(
	                -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

	def moments(self, data):
	    """Returns (height, center_x, center_y, width_x, width_y)
	    the gaussian parameters of a 2D distribution by calculating its
	    moments """
	    #First set all np.nan values in data to 0 
	    #And all the negative values to 0 
	    #These shouldn't influence the moment calculation
	    data[np.isnan(data)] = 0
	    data[data < 0] = 0
	    total = np.nansum(data)
	    print 'The sum over the data is: %s' % total
	    X, Y = indices(data.shape)
	    #print 'The Indices are: %s, %s' % (X, Y)
	    x = np.nansum((X*data))/total
	    y = np.nansum((Y*data))/total
	    #print x, y
	    col = data[:, int(y)]
	    width_x = sqrt(abs((arange(col.size)-y)**2*col).sum()/col.sum())
	    row = data[int(x), :]
	    width_y = sqrt(abs((arange(row.size)-x)**2*row).sum()/row.sum())
	    height = data.max()
	    return height, x, y, width_x, width_y

	def fitgaussian(self, data):
	    """Returns (height, x, y, width_x, width_y)
	    the gaussian parameters of a 2D distribution found by a fit"""
	    params = self.moments(data)
	    errorfunction = lambda p: ravel(self.gaussian(*p)(*indices(data.shape)) -
	                                 data)
	    p, success = optimize.leastsq(errorfunction, params)
	    print p
	    return p

	def psfMask(self):
		"""Returns (FWHM, psfMask) which are the FWHM of the 2D gaussian 
		fit to the collapsed object image and the mask of values found 
		after integrating the gaussian function over all the pixels and 
		normalising by this value. 	
		"""
		#Find the FWHM and the masking profile of a given datacube

		#Step 1 - perform least squares minimisation to find the parameters  
		params = self.fitgaussian(self.imData)
		sigma = params[3]
		FWHM = 2.3548 * sigma
		print 'The FWHM is: %s' % FWHM
		#Step 2 - Return a gaussian function with the fit parameters 
		fit = self.gaussian(*params)
		#Step 3 - Evaluate the gaussian over the pixel range 
		gEval = fit(*indices(self.imData.shape))
		#This is the initial grid of values, but need to normalise to 1 
		#Integrate the gaussian using double quadrature 
		integral = scipy.integrate.dblquad(fit, a=0, b=self.imData.shape[1],\
 			gfun=lambda x: 0 , hfun=lambda x: self.imData.shape[1])
		#Plot the image and the fit 
		colFig, colAx = plt.subplots(1,1, figsize=(12.0,12.0))
		colCax = colAx.imshow(self.imData, interpolation='bicubic')
		colAx.contour(fit(*indices(self.imData.shape)))
		colFig.colorbar(colCax)
		saveName = (self.fileName)[:-5] + '_gauss.png'
		colFig.savefig(saveName)
		plt.show()		
		plt.close('all')
		#return the FWHM and the masked profile 
		return params, (gEval / integral[0]), FWHM, self.offList


##############################################################################
#Uncomment to create test instance of class and try out the methods###########
##############################################################################
#	def testMeth(self):
#		print self.data	

#create class instance 
#cube.specPlot2D(orientation='vertical')
##############################################################################

#data = np.genfromtxt('15names.txt', dtype='str')
#Save the names and types as lists 
#print data
#for name in data:
#	cube = cubeOps(name)
#	cube.specPlot(1)

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








