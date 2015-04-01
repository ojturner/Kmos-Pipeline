#This class houses the methods which are relevant to performing manipulations 
#on the reconstructed and combined data cubes. As such, the class will be a datacube  
#object 


#import the relevant modules
import os, sys, numpy as np, random, matplotlib.mlab as mlab, matplotlib.pyplot as plt, math, matplotlib.cm as cm
import pyraf
import numpy.polynomial.polynomial as poly
import lmfit
from lmfit.models import GaussianModel, ExponentialModel, LorentzianModel, VoigtModel, PolynomialModel
from scipy import stats
from scipy.optimize import minimize
from scipy.optimize import basinhopping
from astropy.io import fits
from matplotlib.backends.backend_pdf import PdfPages
from pylab import *
from matplotlib.colors import LogNorm


####################################################################

class cubeOps(object):
	#Initialiser creates an instance of the cube object
	#Input must be a combined data cube with two extensions - data and noise 
	def __init__(self, fileName):
		self.self = self
		self.fileName = fileName
		#Variable containing all the fits extensions 
		self.Table = fits.open(fileName)
		#Variable housing the primary data cube
		self.data = self.Table[1].data
		#Variable housing the noise data cube
		self.noise = self.Table[2].data
		#Primary Header 
		self.primHeader = self.Table[0].header
		#data Header 
		self.dataHeader = self.Table[1].header
		#noise Header 
		self.noiseHeader = self.Table[2].header
		#Extract the wavelength calibration info 
		self.start_L = self.dataHeader["CRVAL3"]
		self.dL = self.dataHeader["CDELT3"]
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
				self.IFUName = copy(self.IFUNR)
				print 'You have specified a reconstructed type'
			except KeyError:
				raise KeyError("Please Check the type of fits file fed to the class")
		#Create the wavelength array  
		self.wave_array = self.start_L + np.arange(0, 2048*(self.dL), self.dL)
		#Can define all kinds of statistics from the data common to the methods
		#Find the brightest median pixel in the array 
		self.med_array = np.nanmedian(self.data, axis=0)
		self.num = np.nanargmax(self.med_array)
		self.ind1 = self.num / len(self.med_array)
		self.ind2 = self.num - (len(self.med_array) * self.ind1)





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
		plt.plot(self.wave_array, flux_array)
		plt.ylim(0,100)
		saveName = (self.fileName)[:-5] + '.png'
		plt.savefig(saveName)
		plt.show()		
		plt.close('all')
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

		#Check the orientation input to see what has been specified
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













##############################################################################
#Uncomment to create test instance of class and try out the methods###########
##############################################################################
#	def testMeth(self):
#		print self.data	

#create class instance 
#cube = cubeOps('/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/Pipeline_Execution/16-3-15_Min_11Seg/Science_Output/sci_combined_n55_19__telluric.fits')
#cube.specPlot2D(orientation='vertical')
##############################################################################

#data = np.genfromtxt('15names.txt', dtype='str')
#Save the names and types as lists 
#print data
#for name in data:
#	cube = cubeOps(name)
#	cube.specPlot(1)











