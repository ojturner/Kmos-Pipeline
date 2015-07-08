#This class houses the methods which are relevant to performing manipulations 
#on the 1D spectra which have been extracted from the reconstructed data cubes. 



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
from lmfit import Model
from scipy.optimize import minimize
from scipy.optimize import basinhopping
from astropy.io import fits
from matplotlib.backends.backend_pdf import PdfPages
from pylab import *
from matplotlib.colors import LogNorm


####################################################################

class galPhys(object):
	"""
	Def: 
	Class for housing methods relevant to the extraction 
	of physical properties from KMOS galaxy spectra 
	Input: 
	1D spectrum from optimal galaxy extraction
	"""
	#Initialiser creates an instance of the cube object
	#Input must be a combined data cube with two extensions - data and noise 
	def __init__(self, fileName, z):
		"""
		Def:
		Initialiser method 
		Input: filename - name of file containing 1D spectrum 
				z - The redshift of the galaxy being manipulated
		"""
		self.self = self
		#Initialise the fileName object 
		self.fileName = fileName
		self.z = z
		#Variable containing all the fits extensions 
		self.Table = fits.open(fileName)
		#store the flux column
		self.flux = self.Table[1].data['Flux']
		#Store the wavelength column
		self.wavelength = self.Table[1].data['Wavelength']
		#####################################
		#HARDWIRED EMISSION LINES IN MICRONS#
		#####################################
		#Define the wavelength values of the relevant emission lines
		self.OII3727 = 0.3727092
		self.OII3729 = 0.3729875
		self.H_beta = 0.4861363
		self.OIII4959 = 0.4958911
		self.OIII5007 = 0.5006843
		self.H_alpha = 0.6562801
		self.NII6585 = 0.654805
		self.SII6718 = 0.671644
		self.SII6732 = 0.673082
		#Now apply the redshift formula to find where this will be observed
		#Note that for these SDSS spectra the OII doublet is not in range
		self.OII3727_shifted = self.OII3727 * (1 + z)
		self.OII3729_shifted = self.OII3729 * (1 + z)
		self.H_beta_shifted = self.H_beta * (1 + z)
		self.OIII4959_shifted = self.OIII4959 * (1 + z)
		self.OIII5007_shifted = self.OIII5007 * (1 + z)
		self.H_alpha_shifted = self.H_alpha * (1 + z)
		self.NII6585_shifted = self.NII6585 * (1 + z)
		self.SII6718_shifted = self.SII6718 * (1 + z)
		self.SII6732_shifted = self.SII6732 * (1 + z)
		##############################################
		##############################################



	def printProps(self):
		"""
		Def:
		Prints the values of some emission lines in microns 
		and the values of the object flux and wavelength 
		"""
		print 'The wavelength array is: %s' % self.wavelength
		print 'The flux array is: %s' % self.flux
		print 'The shifted OIII Wavelength is: %s' % self.OIII5007_shifted

	def plotSpec(self):
		"""
		Def:
		Quickly plot the spectrum without saving
		"""
		f, ax1 = plt.subplots(1, 1, sharex=True, figsize=(18.0, 10.0))
		ax1.plot(self.wavelength, self.flux, color='b')
		ax1.set_title('Object Spectrum', fontsize=30)
		ax1.set_ylim(0, 4)
		ax1.tick_params(axis='y', which='major', labelsize=15)
		ax1.set_xlabel(r'Wavelength ($\mu m$)', fontsize=24)
		ax1.set_ylabel(r'Flux', fontsize=24)
		f.tight_layout()
		plt.show()

	def fitHbandOIII(self):
		"""
		Def:
		Constructing a composite model for fitting the emission lines 
		OIII doublet and Hb in the K-band filter of high redshift galaxies
		Output: Composite model with fitted parameters 
		"""

		#Need to first fit and subtract the continuum with a polynomial
		poly_mod = PolynomialModel(6)
		#Have an initial guess at the model parameters 
		pars = poly_mod.guess(self.flux, x=self.wavelength)
		#Use the parameters for the full poly_model fit 
		out  = poly_mod.fit(self.flux, pars, x=self.wavelength)
		#The output of the model is the fitted continuum
		continuum = out.best_fit
		#Define the composite model as the sum of three gaussians 
		#First define the individual models, set the parameters 
		OIII5007_gauss = GaussianModel(prefix='OIII5007_')
		pars = OIII5007_gauss.make_params()
		pars['OIII5007_center'].set(self.OIII5007_shifted)
		pars['OIII5007_center'].set(vary=False)
		pars['OIII5007_sigma'].set(0.0008)
		pars['OIII5007_amplitude'].set(0.004)

		OIII4959_gauss = GaussianModel(prefix='OIII4959_')
		pars.update(OIII4959_gauss.make_params())
		pars['OIII4959_center'].set(self.OIII4959_shifted)
		pars['OIII4959_center'].set(vary=False)
		pars['OIII4959_sigma'].set(0.0008)
		pars['OIII4959_amplitude'].set(0.002)

		H_beta_gauss = GaussianModel(prefix='H_beta_')
		pars.update(H_beta_gauss.make_params())
		pars['H_beta_center'].set(self.H_beta_shifted)
		pars['H_beta_center'].set(vary=False)
		pars['H_beta_sigma'].set(0.0008)
		pars['H_beta_amplitude'].set(0.0005)

		#create the composite model as the sum of these three 
		comp_mod = OIII5007_gauss + OIII4959_gauss + H_beta_gauss

		init = comp_mod.eval(pars, x=self.wavelength)
		out = comp_mod.fit(self.flux - continuum, pars, x=self.wavelength)
		print out.fit_report()
		#plot an initial evaluation of the model on top of the spectrum 
		f, ax1 = plt.subplots(1, 1, sharex=True, figsize=(18.0, 10.0))
		ax1.plot(self.wavelength, self.flux - continuum, color='b')
		ax1.plot(self.wavelength, init, 'k--')
		ax1.plot(self.wavelength, out.best_fit, 'r-')
		ax1.set_title('Object Spectrum', fontsize=30)
		ax1.set_ylim(0, 4)
		ax1.tick_params(axis='y', which='major', labelsize=15)
		ax1.set_xlabel(r'Wavelength ($\mu m$)', fontsize=24)
		ax1.set_ylabel(r'Flux', fontsize=24)
		f.tight_layout()
		plt.show()