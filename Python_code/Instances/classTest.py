#This class houses the methods which are relevant 
#For fiting gaussians to a galactic spectra, 
#fitting template spectra to an observed spectrum 
#and estimating the galactic physical properties 

#import the relevant modules
import os, sys, numpy as np, random, matplotlib.mlab as mlab, matplotlib.pyplot as plt, math, matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages
from pylab import *
from matplotlib.colors import LogNorm
import numpy.polynomial.polynomial as poly
import lmfit
from lmfit.models import GaussianModel, ExponentialModel, LorentzianModel, VoigtModel, PolynomialModel

#add my toolkit.py file to the PYTHONPATH
sys.path.append('/Users/owenturner/Documents/PhD/KMOS/Analysis_Pipeline/Python_code/Class')

from astropy.io import fits
from pipelineClass import pipelineOps

objFile = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/NGC55/14-9-2014/KMOS_SPEC_OBS258_0007.fits'
skyFile = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/NGC55/14-9-2014/KMOS_SPEC_OBS258_0008.fits'
badPMap = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/8_12_14Cal_products/badpixel_dark.fits'
lcalMap = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/8_12_14Cal_products/lcal_YJYJYJ.fits'

pipe_methods = pipelineOps()


pipe_methods.computeOffsetSegments(objFile, skyFile, badPMap, lcalMap)

newFile = 'newCorrection.fits'

#Now try the subtraction method 
#pipe_methods.subFrames(newFile, skyFile)
print 'all good'
print 'Hello'
#pipe_methods.pixelHistogram('KMOS_SPEC_OBS258_0007_m_8_raw.fits', 'topFour7_m_8.fits', 780, 1270)

#pipe_methods.stackLcal(lcalMap)
#pipe_methods.computeOffsetTopFour('KMOS_SPEC_OBS258_0007_m_8_raw.fits', objFile)
#pipe_methods.subFrames('topFourCorrected7.fits', skyFile)