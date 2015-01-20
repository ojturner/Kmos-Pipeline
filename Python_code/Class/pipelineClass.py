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
sys.path.append('/Users/owenturner/Documents/PhD/SUPA_Courses/AdvancedDataAnalysis/Homeworks')

#and import the toolkit
import toolkit
from astropy.io import fits

####################################################################

class pipelineOps(object):
	#Initialiser creates an instance of the spectrumFit object 
	def __init__(self):
		self.self = self

############################################################################################
#Current functions in this file: 
#computeOffset - from a sky file, object file, bad pixel map and full lcal file 
#                compute the column correction for the 2048x64 pixel chunks 
#
#computeOffsetTwo - Same as above, but uses the output of the stackLcal function 
#                   which are a more conservative set of stacked lcal maps for each detector
#
#computeOffsetTopFour - Use the top four pixels on the raw subtracted frame to compute the offset
#
#subFrames       - Take one fits image with three data extensions and subtract another 
#
#pixelHistogram  - select a chunk of pixels and compute and draw a histogram of their values
#
#stackLcal       - stack the lcal frames from the six different rotation angles. Makes sure 
#                   that the median pixel value will be more reliable 
#
############################################################################################

######################################################################################
#MODULE: computeOffset
#
#PURPOSE:
#Take the object image after calibration, along with a sky frame, bad pixel frame and 
#the lcal frame and homogenise the readout columns, so that after subtraction gives 0
#
#INPUTS:
#
#			objectFile: The object image to be corrected 	
#			skyFile: The sky image taken directly before or after the object frame
#			badPMap: The bad pixel frame generated by the pipeline
#			lcalMap: The calibration frame generated by the pipeline
#
#OUTPUTS: 	newObjData: The corrected 2048x2048 object arrray which can then be saved
#						
#
#USAGE: 	newObjData = computeOffset(objectFile, skyFile, badPMap, lcalMap)
#######################################################################################	

	
	#Access the primary extension, later this will be looped  
	def computeOffset(self, objectFile, skyFile, badPMap, lcalMap): 
		#Set up vector to house the corrected extensions
		correctedExtensions=[]
		#Read in the tables of data
		table_o = fits.open(objectFile)
		fitsHeader = table_o[0].header
		print fitsHeader
		table_s = fits.open(skyFile)
		bad_pixel_table = fits.open(badPMap)

		#Now choose the correct rotation angle 
		lcal_table = fits.open(lcalMap)
		#This is a list of all possible rotation angles 
		angleList = np.array([0, 60, 120, 180, 240, 300])
		#Select the ocs.rot.naangle keyword 
		obsAngle = table_o[0].header["HIERARCH ESO OCS ROT NAANGLE"]
		#Find where the difference between the observed and idealised angle is minimum
		newAngleList = obsAngle - angleList
		n = newAngleList.argmin()
		obsAngleNew = angleList[n]


		#Find the extension to which this corresponds
		val = 0
		if obsAngleNew == 0:
			val = 1 
		elif obsAngleNew == 60:
			val = 4
		elif obsAngleNew == 120:
			val = 7
		elif obsAngleNew == 180: 
			val = 10
		elif obsAngleNew == 240:
			val = 13
		elif obsAngleNew == 300:
			val = 16		
		print val	


		#Loop over the fits image extensions, do the same each time
		for count in range(1,4):
			print val
			data_o = table_o[count].data
			data_s = table_s[count].data
			bad_pixel_data = bad_pixel_table[count].data			
			lcal_data = lcal_table[val].data
			#Debug to see if mask is being applied properly
			test_array = np.zeros(shape=[2048, 2048])
			tempName = 'lcal' + str(count) + '.fits'
			#Loop around the bad pixel locations

			#Now we have both the object and sky data these are 2D arrays, 
			#essentially a matrix where each number represents a pixel flux, 
			#And the location of the number in the matrix represents the 
			#Pixel position on the detector, which in turn corresponds to the 
			#objects position on the sky. Need to slice this 2D array into a 
			#1D array of 2D arrays, each of which is 64 pixels wide 
			#so that I can examine these in turn and loop over them 

			#Counters for the slicing 
			x = 0
			y = 64

			#1D arrays to host the data 
			skyArray = []
			objArray = []
			badPArray = []
			lcalArray = []
			testArray = []

			for i in range(32):
			   
			   #Slice each of the data files into 32 columns of 64 pixels width
			   newObjectArray = data_o[:,x:y]
			   newSkyArray = data_s[:,x:y]
			   newPArray = bad_pixel_data[:,x:y]
			   newCalArray = lcal_data[:,x:y]
			   newTestArray = test_array[:,x:y]
			   testArray.append(newTestArray)
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
				lcal_pixel_coords = np.where(lcalArray[num] > 0)


				for i in range(len(bad_pixel_coords[0])):
					#Because of the way np.where works, need to define the x and y coords in this way
					xcoord = bad_pixel_coords[0][i]
					ycoord = bad_pixel_coords[1][i]
					#Now set all positions where there is a dead pixel to np.nan in the object and sky
					testArray[num][xcoord][ycoord] = np.nan

				for i in range(len(lcal_pixel_coords[0])):
					#Because of the way np.where works, need to define the x and y coords in this way
					xcoord = lcal_pixel_coords[0][i]
					ycoord = lcal_pixel_coords[1][i]
					#Now set all positions where there is a dead pixel to np.nan in the object and sky
					testArray[num][xcoord][ycoord] = np.nan

				#Save test_array as a fits file 
				fits.writeto(tempName, data=test_array, clobber=True)	
				#Loop through the first 2048 x 64 data and sky columns and mask off these coords
				#Then compute the median in both columns. If median sky > median obj, add the difference 
				#between the median values to the obj. If the other way around, decrement

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
					#Do the same, this time for the slitlet positions (substantially more will have a value)
					xcoord = lcal_pixel_coords[0][i]
					ycoord = lcal_pixel_coords[1][i]
					#Set all of these locations to nan 
					objTemp[num][xcoord][ycoord] = np.nan
					skyTemp[num][xcoord][ycoord] = np.nan


				#Now all the pixels that shouldn't be included in the median 
				#have value nan. Can then just do np.nanmean(objTemp) which will ignore nan's
				#then repeat the process for the sky, compare the mean's, compute and apply the offset.

				obj_mean = np.nanmedian(objTemp[num])
				sky_mean = np.nanmedian(skyTemp[num])
				print sky_mean
				print obj_mean

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
			correctedExtensions.append(newObjData)
			if count == 1:
				print 'Computed First Correction'
			elif count == 2:
				print 'Computed Second Correction'
			elif count == 3: 
				print 'Computed Third Correction'			
			val += 1				
		#Create the object fits file with the three corrected extensions
		fileName = raw_input('Enter a name for the corrected fits file: ') + '.fits'
		#Note that the readout complained about the header not being 
		#in the correct fits format
		fits.writeto(fileName, data = [], header=fitsHeader, clobber=True)
		fits.append(fileName, data=correctedExtensions[0])	
		fits.append(fileName, data=correctedExtensions[1])	
		fits.append(fileName, data=correctedExtensions[2])

	#Second compute offset method, this time with the refined lcal
	#Access the primary extension, later this will be looped  
	def computeOffsetTwo(self, objectFile, skyFile, badPMap, lcal1, lcal2, lcal3): 
		#Set up vector to house the corrected extensions
		correctedExtensions=[]
		#Read in the tables of data
		table_o = fits.open(objectFile)
		fitsHeader = table_o[0].header
		print fitsHeader
		table_s = fits.open(skyFile)
		bad_pixel_table = fits.open(badPMap)

		#Now read in the refined lcal maps 
		lcal1 = fits.open(lcal1)
		lcal2 = fits.open(lcal2)
		lcal3 = fits.open(lcal3)

		#Save to a dictionary 
		d = {}
		d[1] = lcal1[0].data
		d[2] = lcal2[0].data
		d[3] = lcal3[0].data

		#Loop over the fits image extensions, do the same each time
		for count in range(1,4):
			data_o = table_o[count].data
			data_s = table_s[count].data
			bad_pixel_data = bad_pixel_table[count].data			
			lcal_data = d[count]

			#Loop around the bad pixel locations
			#Now we have both the object and sky data these are 2D arrays, 
			#essentially a matrix where each number represents a pixel flux, 
			#And the location of the number in the matrix represents the 
			#Pixel position on the detector, which in turn corresponds to the 
			#objects position on the sky. Need to slice this 2D array into a 
			#1D array of 2D arrays, each of which is 64 pixels wide 
			#so that I can examine these in turn and loop over them 

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
				lcal_pixel_coords = np.where(lcalArray[num] == 0)


				#Loop through the first 2048 x 64 data and sky columns and mask off these coords
				#Then compute the median in both columns. If median sky > median obj, add the difference 
				#between the median values to the obj. If the other way around, decrement

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
					#Do the same, this time for the slitlet positions (substantially more will have a value)
					xcoord = lcal_pixel_coords[0][i]
					ycoord = lcal_pixel_coords[1][i]
					#Set all of these locations to nan 
					objTemp[num][xcoord][ycoord] = np.nan
					skyTemp[num][xcoord][ycoord] = np.nan


				#Now all the pixels that shouldn't be included in the median 
				#have value nan. Can then just do np.nanmean(objTemp) which will ignore nan's
				#then repeat the process for the sky, compare the mean's, compute and apply the offset.

				obj_mean = np.nanmedian(objTemp[num])
				sky_mean = np.nanmedian(skyTemp[num])
				print sky_mean
				print obj_mean

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
			correctedExtensions.append(newObjData)
			if count == 1:
				print 'Computed First Correction'
			elif count == 2:
				print 'Computed Second Correction'
			elif count == 3: 
				print 'Computed Third Correction'			
		#Create the object fits file with the three corrected extensions
		fileName = raw_input('Enter a name for the corrected fits file: ') + '.fits'
		#Note that the readout complained about the header not being 
		#in the correct fits format
		fits.writeto(fileName, data = [], header=fitsHeader, clobber=True)
		fits.append(fileName, data=correctedExtensions[0])	
		fits.append(fileName, data=correctedExtensions[1])	
		fits.append(fileName, data=correctedExtensions[2])	

	def computeOffsetTopFour(self, rawSubFile, objectFile):

		correctedExtensions = []
		#Read in the tables of data
		table_o = fits.open(objectFile)
		fitsHeader = table_o[0].header
		print fitsHeader
		table_s = fits.open(rawSubFile)	


		#Loop over the fits image extensions, do the same each time
		for count in range(1,4):
			data_o = table_o[count].data
			data_s = table_s[count].data

			#Counters for the slicing 
			x = 0
			y = 64

			#1D arrays to host the data 
			subArray = []
			objArray = []

			for i in range(32):
			   
			   #Slice each of the data files into 32 columns of 64 pixels width
			   newObjectArray = data_o[:,x:y]
			   newSubArray = data_s[2044:2048,x:y]

			   objArray.append(newObjectArray)
			   subArray.append(newSubArray)

			   #Add 64 to the counters each time to create the slices
			   x += 64
			   y += 64		

			testData = np.hstack(subArray)
			fileName = 'testTopFour' + str(count) + '.fits'
			fits.writeto(fileName, data=testData, clobber=True)   

			for num in range(len(objArray)):
				correctionMedian = np.median(subArray[num])
				objArray[num] -= correctionMedian
				print correctionMedian
				print subArray[num][:,0:1]
			#Have now made the correction to all 32 of the 64 pixel width columns.
			#All that is left to do is stitch the objArray arrays back together to 
			#give a single 2048x2048 array.
			newObjData = np.hstack(objArray)
			correctedExtensions.append(newObjData)
			if count == 1:
				print 'Computed First Correction'
			elif count == 2:
				print 'Computed Second Correction'
			elif count == 3: 
				print 'Computed Third Correction'			
		#Create the object fits file with the three corrected extensions
		fileName = raw_input('Enter a name for the corrected fits file: ') + '.fits'
		#Note that the readout complained about the header not being 
		#in the correct fits format
		fits.writeto(fileName, data = [], header=fitsHeader, clobber=True)
		fits.append(fileName, data=correctedExtensions[0])	
		fits.append(fileName, data=correctedExtensions[1])	
		fits.append(fileName, data=correctedExtensions[2])				
			


	def subFrames(self, objectFile, skyFile):	

		#Read in the object and sky files 
		objData = fits.open(objectFile)
		skyData = fits.open(skyFile)

		#Find the header and extensions of the new fits file
		header = objData[0].header
		print header
		ext1 = objData[1].data - skyData[1].data
		ext2 = objData[2].data - skyData[2].data
		ext3 = objData[3].data - skyData[3].data

		#Write out to a different fits file, with name user specified
		nameOfFile = raw_input('Enter the name of the subtracted file: ')
		nameOfFile = nameOfFile + '.fits'
		fits.writeto(nameOfFile, data=[], header=header, clobber=True)
		fits.append(nameOfFile, data=ext1)	
		fits.append(nameOfFile, data=ext2)	
		fits.append(nameOfFile, data=ext3)

	def pixelHistogram(self, subFile, subCorFile, x1, x2):

		#Create a histogram of pixel values on the 
		#subtracted frame before and after correction	

		#First read in the files 
		subData = fits.open(subFile)
		subCorData = fits.open(subCorFile)

		#At the moment we'll just consider the first extension for our data
		subData = subData[1].data
		subCorData = subCorData[1].data

		#The input numbers define the left and right edges of the pixel section
		subData = subData[:,x1:x2]
		subCorData = subCorData[:,x1:x2]
		
		#This gives an array of arrays, we just care about the numbers, not spatial info
		#Use ravel() to convert these into lists for the histogram
		subData = subData.ravel()
		subCorData = subCorData.ravel()	
		print len(subData)
		print subData
		print np.median(subData)
		print np.median(subCorData)

		#Create the bins array for both histograms 
		bins = np.arange(-15, 15, 1)

		plt.close('all')
		#Plot the histograms 
		n1, bins1, patches1 = plt.hist(subData, bins=bins, histtype='step',\
		 color='green', linewidth=3, label='Before Correction') 

		n2, bins2, patches2 = plt.hist(subCorData, bins=bins, histtype='step',\
		 color='blue', linewidth=3, alpha=0.5, label='After Correction')

		#Now want to fit the gaussians to the histograms using lmfit gaussian models 
		#gaussian model number 1
		mod1 = GaussianModel()

		#Create a new bins vector for the fit 
		fitBins = np.arange(-14.5, 14.5, 1)
		print len(fitBins)

		#Take an initial guess at what the model parameters are 
		#In this case the gaussian model has three parameters, 
		#Which are amplitude, center and sigma
		pars1 = mod1.guess(n1, x=fitBins)
		#Perform the actual fit 
		out1  = mod1.fit(n1, pars1, x=fitBins)
		#Now want to add this curve to our plot 
		plt.plot(fitBins, out1.best_fit, linewidth=2.0, label='b.c. model', color='green')

		#Repeat for the corrected data
		mod2 = GaussianModel()
		pars2 = mod2.guess(n2, x=fitBins)
		out2  = mod2.fit(n2, pars2, x=fitBins)
		plt.plot(fitBins, out2.best_fit, linewidth=2.0, label='a.c. model', color='blue')
		plt.xlabel('Counts per pixel')
		plt.ylabel('Number per bin')
		plt.title('Subtracted Frame, Improvement after Column Correction')
		plt.legend(loc='upper left', fontsize='small')
		#plt.savefig(raw_input('Enter the plot name: '), papertype='a4', orientation='landscape')
		plt.show()

		#Print out the fit reports to look at the centre at S.D. of each model
		print out1.fit_report() 
		print out2.fit_report() 

	#def topFourCorrection()	
	def stackLcal(self, lcalFile):
		#Read in the file 
		lcal = fits.open(lcalFile)
		#Want to have it so that only nan values survive. Can do this 
		#By substituting 'False' values for the numbers and true values for nan
		#So that when multiplied only True * True all the way will survive 

		#Loop round all the extensions and change the np.nan values to True 
		#and the pixels with values to False 
		d = {}
		for i in range(1, 19):
			data = lcal[i].data 

			#Define the coordinates where there is a value 
			value_coords = np.where(data > 0)

			#Define where there is np.nan
			nan_coords = np.where(np.isnan(data))
			temp = np.empty(shape=[2048,2048], dtype=float)
			#loop over the pixel values in data and change to either True or False 
			for j in range(len(value_coords[0])):
				xcoord = value_coords[0][j]
				ycoord = value_coords[1][j]
				temp[xcoord][ycoord] = 0

			for j in range(len(nan_coords[0])):
				xcoord = nan_coords[0][j]
				ycoord = nan_coords[1][j]
				temp[xcoord][ycoord] = 1

			d[i] = temp
			print 'done' + str(i)
		print len(d)
		#Check to see if we're creating the right array, which we are	
		#fits.writeto(filename='test.fits', data=d[1], clobber=True)	

		#Now create the individual lcal frames by stacking the extensions together
		#This gives a more conservative estimate of the pixels we should and should 
		#not be using throughout the data analysis stage 

		lcal1 = d[1] * d[4] * d[7] * d[10] * d[13] * d[16]
		lcal2 = d[2] * d[5] * d[8] * d[11] * d[14] * d[17]
		lcal3 = d[3] * d[6] * d[9] * d[12] * d[15] * d[18]

		fits.writeto(filename='lcal1.fits', data=lcal1, clobber=True)
		fits.writeto(filename='lcal2.fits', data=lcal2, clobber=True)
		fits.writeto(filename='lcal3.fits', data=lcal3, clobber=True)
