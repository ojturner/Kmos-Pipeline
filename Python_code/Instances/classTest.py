#Don't need to import all of the same things which are imported in the class. Done already.

#import the relevant modules
import os, sys


#add the class file to the PYTHONPATH
sys.path.append('/Users/owenturner/Documents/PhD/KMOS/Analysis_Pipeline/Python_code/Class')

#import the class 
from pipelineClass import pipelineOps

objFile = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/NGC55/14-9-2014/KMOS_SPEC_OBS258_0009.fits'
skyFile = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/NGC55/14-9-2014/KMOS_SPEC_OBS258_0008.fits'
badPMap = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/8_12_14Cal_products/badpixel_dark.fits'
lcalMap = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/8_12_14Cal_products/lcal_YJYJYJ.fits'

#Create an instance of the class
pipe_methods = pipelineOps()

#The following are examples of using the functions within the class
#pipe_methods.computeOffsetSegments(objFile, skyFile, badPMap, lcalMap)

newFile = 'newCorrection.fits'

#Now try the subtraction method 
#pipe_methods.subFrames('segmentsCorrected_9_128.fits', skyFile)
print 'all good'
print 'Hello'
pipe_methods.pixelHistogram('KMOS_SPEC_OBS258_0009_m_8_raw.fits', 'segmentsSubtracted_9_m_8_128.fits', 780, 1270)

#pipe_methods.stackLcal(lcalMap)
#pipe_methods.computeOffsetTopFour('KMOS_SPEC_OBS258_0007_m_8_raw.fits', objFile)
#pipe_methods.subFrames('topFourCorrected7.fits', skyFile)