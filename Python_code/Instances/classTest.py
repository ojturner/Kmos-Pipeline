#Don't need to import all of the same things which are imported in the class. Done already.

#import the relevant modules
import os, sys


#add the class file to the PYTHONPATH
sys.path.append('/Users/owenturner/Documents/PhD/KMOS/Analysis_Pipeline/Python_code/Class')

#import the class 
from pipelineClass import pipelineOps

objFile = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/NGC55/14-9-2014/KMOS_SPEC_OBS258_0009_Corrected.fits'
skyFile = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/NGC55/14-9-2014/KMOS_SPEC_OBS258_0008.fits'
badPMap = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/Pipeline_Execution/8-12-14/Calibration_Files/badpixel_dark.fits'
lcalMap = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/Pipeline_Execution/26-1-15/Calibration_Files/lcal_YJYJYJ.fits'

#Create an instance of the class
pipe_methods = pipelineOps()

#pipe_methods.computeOffsetSegments(objFile, skyFile, badPMap, lcalMap)
#The following are examples of using the functions within the class
#pipe_methods.computeOffsetTopFour('KMOS_SPEC_OBS258_0001_m_2_raw.fits', objFile)

newFile = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/NGC55/14-9-2014/KMOS_SPEC_OBS258_0003_Corrected.fits'

#Now try the subtraction method 
#pipe_methods.subFrames('topFour_1_Corrected.fits', skyFile)
#pipe_methods.subFrames('KMOS_SPEC_OBS258_0001_Corrected.fits', skyFile)
#print 'all good'
#print 'Hello'
#pipe_methods.pixelHistogram('KMOS_SPEC_OBS258_0003_m_2_raw.fits', 'KMOS_SPEC_OBS258_0003_m_2.fits', 780, 1270)

#pipe_methods.stackLcal(lcalMap)
#pipe_methods.computeOffsetTopFour('KMOS_SPEC_OBS258_0007_m_8_raw.fits', objFile)
#pipe_methods.subFrames('topFourCorrected7.fits', skyFile)
#pipe_methods.applyCorrection('NGC55_14-9-2014_fileNames.txt', badPMap, lcalMap)
#pipe_methods.plotMedian('KMOS_SPEC_OBS258_0009_m_8_raw.fits', 'KMOS_SPEC_OBS258_0009_m_8.fits', \
#	'segmentsSubtracted_9_m_8_128.fits', 'topFour_9_m_8.fits', 1000, 1200, 800, 1270)
#pipe_methods.plotMedian('KMOS_SPEC_OBS258_0009_m_8_raw.fits', 'KMOS_SPEC_OBS258_0009_m_8.fits', \
#	'segmentsSubtracted_9_m_8_128.fits', 'topFour_9_m_8.fits', 1900, 2100, 800, 1270)	
#pipe_methods.plotMedian('KMOS_SPEC_OBS258_0009_m_8_raw.fits', 'KMOS_SPEC_OBS258_0009_m_8.fits', \
# 'segmentsSubtracted_9_m_8_128.fits', 'topFour_9_m_8.fits', 400, 600, 800, 1270)


rho = pipe_methods.crossCorr(ext=1, objFile=objFile, skyFile=skyFile, badpmap=badPMap, y1=100, y2=1800, x1=100, x2=1800)
print rho
#pipe_methods.shiftImage(ext=1, infile=objFile, skyfile=skyFile, interp_type = 'linear',stepsize=0.03, xmin=-0.2, xmax=0.2, ymin=-0.2, ymax=0.2)
#pipe_methods.rotateImage(ext=1, infile=objFile, skyfile=skyFile, interp_type = 'linear', stepsize=0.002, minAngle=-0.1, maxAngle=0.1)
#array = pipe_methods.imSplit(ext=1, infile=objFile, vertSegments=4, horSegments=4)
#pipe_methods.shiftImageSegments(ext=1, infile=objFile, skyfile=skyFile,\
#  	 vertSegments=2, horSegments=2, interp_type='spline3', \
#  	 stepsize=0.02, xmin=-0.1, xmax=0.1, ymin=-0.1, ymax=0.1)












