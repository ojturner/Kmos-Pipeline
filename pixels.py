#import the relevant modules
import os, sys, numpy as np, random, matplotlib.mlab as mlab, matplotlib.pyplot as plt, math, matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages
from pylab import *
from matplotlib.colors import LogNorm
import numpy.polynomial.polynomial as poly
import lmfit
from lmfit.models import GaussianModel, ExponentialModel, LorentzianModel, VoigtModel, PolynomialModel

#add my toolkit.py file to the PYTHONPATH
sys.path.append('/Users/owenturner/Documents/PhD/SUPA_Courses/AdvancedDataAnalysis/Homeworks')

#and import the toolkit
import toolkit
from astropy.io import fits

#Ultimately want to be looping round pairs of input files 
#And the different extensions within these, will have to try and code clever 

#First open the fits files containing the sky image and the object detection 
#Will ultimately have to check from a list of input fits files which are which 
#And do this procedure for every pair of input files 
#At the moment will try to write something that works and can be scaled up later

#Read in the object and sky fits images
#Will hard wire in the names at the moment
#Access the primary extension, later this will be looped  

table_o = fits.open("KMOS_SPEC_OBS258_0009.fits")
data_o = table_o[1].data

table_s = fits.open("KMOS_SPEC_OBS258_0008.fits")
data_s = table_s[1].data

bad_pixel_table = fits.open('badpixel_dark.fits')
bad_pixel_data = bad_pixel_table[1].data

lcal_table = fits.open('lcal_YJYJYJ.fits')
lcal_data = lcal_table[8].data

#Now we have both the object and sky data these are 2D arrays, 
#essentially a matrix where each number represents a pixel flux, 
#And the location of the number in the matrix represents the 
#Pixel position on the detector, which in turn corresponds to the 
#objects position on the sky. Need to slice this 2D array into a 
#1D array of 2D arrays, each of which is 64 pixels wide 
#So that I can examine these in turn and loop over them 

#Counters for the slicing 
x = 0
y = 64

#1D arrays to host the data 
skyArray = []
objArray = []
badPArray = []
lcalArray = []

for i in range(32):
   
   #Slice each of the data files into 32 columns of 64 pixels width
   newObjectArray = data_o[:,x:y]
   newSkyArray = data_s[:,x:y]
   newPArray = bad_pixel_data[:,x:y]
   newCalArray = lcal_data[:,x:y]
   objArray.append(newObjectArray)
   skyArray.append(newSkyArray)
   badPArray.append(newPArray)
   lcalArray.append(newCalArray)

   #Add 64 to the counters each time to create the slices
   x += 64
   y += 64

#Now wanto to read in the lcal and the bad pixel map, 
#To get a list of the pixel indices which should be averaged
#First need to slice this in the same way as the other file 
#Then simply find the bad pixel and the lcal positions, hide 
#these as nan and compute the median in each 64 pixel column 
#from everything that isn't nan, then apply the correction to the 
#data before stitching everything back together as a data file
#And feeding back into the pipeline at the appropriate location

#Redefine objArray and sky array, do manipulations to temp
objTemp = copy(objArray)
skyTemp = copy(skyArray)
#Start the loop for all the columns in the Array vectors 
for num in range(len(objArray)):


	#Find the coordinates of the bad pixels and the slitlets 
	bad_pixel_coords = np.where(badPArray[num] == 0)
	lcal_pixel_coords = np.where(np.isnan(lcalArray[num]))

	#Loop through the first 2048 x 64 data and sky columns and mask off these coords
	#Then compute the median in both columns. If median sky > median obj, add the difference 
	#between the median values to the obj. If the other way around, decrement
	#Need first to create a temporary vector for each of the arrays:


	#Loop around the bad pixel locations
	for i in range(len(bad_pixel_coords[0])):
		#Because of the way np.where works, need to define the x and y coords in this way
		xcoord = bad_pixel_coords[0][i]
		ycoord = bad_pixel_coords[1][i]
		#Now set all positions where there is a dead pixel to np.nan in the object and sky
		objTemp[num][xcoord][ycoord] = np.nan
		skyTemp[num][xcoord][ycoord] = np.nan

	#Loop around the slitlet positions
	for i in range(len(lcal_pixel_coords[0])):
		#Do the same, this time for the slitlet positions (substantially more will be nan)
		xcoord = lcal_pixel_coords[0][i]
		ycoord = lcal_pixel_coords[1][i]
		#Set all of these locations to nan 
		objTemp[num][xcoord][ycoord] = np.nan
		skyTemp[num][xcoord][ycoord] = np.nan


	#hdu = fits.PrimaryHDU(objTemp)
	#hdu.writeto('masked.fits')
	#Now all the pixels that shouldn't be included in the median 
	#Should have value nan. Can then just do np.nanmean(objTemp) which will ignore nan's
	#Then repeat the process for the sky, compare the mean's, compute and apply the offset.

	obj_mean = np.nanmean(objTemp[num])
	sky_mean = np.nanmean(skyTemp[num])
	print obj_mean
	print sky_mean

	#Need to compare the two medians to see how to apply the offset.
	#If the sky is brighter, add the difference to the object image

	if sky_mean > obj_mean:
		objArray[num] += abs(sky_mean - obj_mean)
	elif obj_mean > sky_mean:
		objArray[num] -= abs(obj_mean - sky_mean)	

#Have now made the correction to all 32 of the 64 pixel width columns.
#All that is left to do is stitch the objArray arrays back together to 
#give a single 2048x2048 array.
newObjData = np.hstack(objArray)
#hdu = fits.PrimaryHDU(newObjData)
#hdu.writeto('corrected_object9.fits')

#Now perform the subtraction and see if this looks better
objMinusSky = newObjData - data_s
hdu = fits.PrimaryHDU(objMinusSky)
hdu.writeto('objMinusSkyCorrected9.fits')



	
