#Don't need to import all of the same things which are imported in the class. Done already.

#import the relevant modules
import os, sys, numpy as np, matplotlib.pyplot as plt
from astropy.io import fits


#add the class file to the PYTHONPATH
sys.path.append('/Users/owenturner/Documents/PhD/KMOS/Analysis_Pipeline/Python_code/Class')

#import the class 
from pipelineClass import pipelineOps
from cubeClass import cubeOps


oldbadP = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/Pipeline_Execution/8-12-14/Calibration_Files/badpixel_dark.fits'
oldlcal = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/Pipeline_Execution/8-12-14/Calibration_Files/lcal_YJYJYJ.fits'

objFile = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/NGCLEE/K-band/KMOS.2014-08-05T03_14_39.263.fits'
skyFile = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/NGCLEE/K-band/KMOS.2014-08-05T03_15_35.245.fits'
badPMap14 = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/Pipeline_Execution/16-3-15_Min_11Seg/Calibration_Files/badpixel_dark_Added.fits'
badPMap15 = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/Pipeline_Execution/S24-3-15/Calibration_Files/badpixel_dark_Added.fits'
lcalMap14 = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/Pipeline_Execution/16-3-15_Min_11Seg/Calibration_Files/lcal_YJYJYJ.fits'
lcalMap15 = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/Pipeline_Execution/S24-3-15/Calibration_Files/lcal_YJYJYJ.fits'

hskyCube = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/NGCLEE/H-band/Science/cube_science.fits'
kskyCube = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/NGCLEE/K-band/Science/cube_science.fits'


#Create an instance of the class
pipe_methods = pipelineOps()
sky_cube = cubeOps(kskyCube)


#pipe_methods.computeOffsetSegments(objFile, skyFile, badPMap, lcalMap)
#The following are examples of using the functions within the class
#pipe_methods.computeOffsetTopFour('KMOS_SPEC_OBS258_0001_m_2_raw.fits', objFile)
k_names = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/NGCLEE/K-band/Calibrations/corrected_object_names.txt'
h_names = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/NGCLEE/H-band/Calibrations/corrected_object_names.txt'


#Permanent Names
raw_14 = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/NGC55/14-9-2014/NGC55_14_Names.txt'
names_14 = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/NGC55/14-9-2014/NGC55_14_Corrected_Names.txt'
names_14_short = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/NGC55/14-9-2014/NGC55_14_Corrected_Names_short.txt'
names_14_shifted = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/NGC55/14-9-2014/NGC55_14_Corrected_Names_shifted.txt'

raw_15 = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/NGC55/15-09-2014/NGC55_15_Names.txt'
names_15 = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/NGC55/15-09-2014/NGC55_15_Corrected_Names.txt'
names_15_short = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/NGC55/15-09-2014/NGC55_15_Corrected_Names_short.txt'
names_15_shifted = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/NGC55/15-09-2014/NGC55_15_Corrected_Names_shifted.txt'


noTel_names_15 = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/NGC55/15-09-2014/NGC55_15_Corrected_Names_noTel.txt'
noTel_names_14 = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/NGC55/14-9-2014/NGC55_14_Corrected_Names_noTel.txt'

#Changes depending on reduction process
frame_check_names = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/Pipeline_Execution/frameCheck/combNames.txt'
sci_names_14 = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/Pipeline_Execution/S24-3-15/5sig_Science_Output'
sci_names_15 = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/Pipeline_Execution/S24-3-15/5sig_Science_Output'
sci_names_14_noTel = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/Pipeline_Execution/10-4-15_Pairs_14/all_but_9_ws/sci_names_noTel.txt'
sci_names_14_noTel_1 = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/Pipeline_Execution/10-4-15_Pairs_14/all_but_9/sci_names_noTel.txt'
sci_names_15_noTel = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/Pipeline_Execution/10-4-15_250-1500_15/Science_Output/sci_names_noTel.txt'


newFile2 = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/NGC55/15-09-2014/KMOS_SPEC_OBS259_0014_Corrected_22_spline3_Shifted.fits'
newFile3 = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/NGC55/14-9-2014/KMOS_SPEC_OBS258_0009_Corrected_Subtracted.fits'

namesOfFile = np.genfromtxt(sci_names_14_noTel, dtype='str')
namesOfFile_1 = np.genfromtxt(sci_names_14_noTel_1, dtype='str')

objCube = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/Pipeline_Execution/10-4-15_Pairs_14/all_but_9/sci_combined_n55_19__skytweak.fits'
objSpec2 = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/Pipeline_Execution/10-4-15_Pairs_14/all_but_9/sci_combined_n55_19__skytweak_spectrum.fits'
objSpec = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/NGCLEE/H-band/Science/Best_sci_combined_P108__skytweak_spectrum.fits'
skySpec = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/NGCLEE/H-band/Science/cubesky_spectrum.fits'

hobjframe = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/NGCLEE/H-band/raw_frames/KMOS.2014-08-03T00:05:24.218_Corrected_11_spline3_Shifted.fits'
hskyframe = '/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/NGCLEE/H-band/raw_frames/KMOS.2014-08-03T00:03:33.904.fits'

#pipe_methods.saveSpec(hskyCube)
#pipe_methods.plotSpecs(objSpec, skySpec, 1)

#pipe_methods.multiExtractSpec(skyCube=hskyCube, frameNames=h_names)
#Now try the subtraction method 
#pipe_methods.subFrames('/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/NGC55/14-9-2014/KMOS_SPEC_OBS258_0009_Corrected_11_spline3_Shifted.fits', skyFile)
#pipe_methods.subFrames('KMOS_SPEC_OBS258_0001_Corrected.fits', skyFile)
#print 'all good'
#print 'Hello'
#pipe_methods.pixelHistogram('KMOS_SPEC_OBS258_0001_m_2_raw.fits', 'KMOS_SPEC_OBS258_0001_m_2.fits', 780, 1270)

#pipe_methods.stackLcal(lcalMap)
#pipe_methods.computeOffsetTopFour('KMOS_SPEC_OBS258_0007_m_8_raw.fits', objFile)
#pipe_methods.subFrames(objFile, skyFile)
#pipe_methods.applyCorrection(raw_14, badPMap14, lcalMap14)
#pipe_methods.plotMedian('KMOS_SPEC_OBS258_0009_m_8_raw.fits', 'KMOS_SPEC_OBS258_0009_m_8.fits', \
#	'segmentsSubtracted_9_m_8_128.fits', 'topFour_9_m_8.fits', 1000, 1200, 800, 1270)
#pipe_methods.plotMedian('KMOS_SPEC_OBS258_0009_m_8_raw.fits', 'KMOS_SPEC_OBS258_0009_m_8.fits', \
#	'segmentsSubtracted_9_m_8_128.fits', 'topFour_9_m_8.fits', 1900, 2100, 800, 1270)	
#pipe_methods.plotMedian('KMOS_SPEC_OBS258_0009_m_8_raw.fits', 'KMOS_SPEC_OBS258_0009_m_8.fits', \
# 'segmentsSubtracted_9_m_8_128.fits', 'topFour_9_m_8.fits', 400, 600, 800, 1270)

#pipe_methods.maskFile(skyFile, badPMap)
#pipe_methods.badPixelextend(badpmap=badPMap)
#rho = pipe_methods.crossCorr(ext=1, objFile=objFile, skyFile=skyFile, badpmap=badPMap, y1=100, y2=1800, x1=100, x2=1800)
#print rho
#pipe_methods.shiftImage(ext=1, infile=objFile, skyfile=skyFile, badpmap=badPMap, \
# interp_type = 'spline3', stepsize=0.01, xmin=-0.1, xmax=0.1, ymin=-0.1, ymax=0.1)
#pipe_methods.rotateImage(ext=1, infile=objFile, skyfile=skyFile, interp_type = 'linear', stepsize=0.002, minAngle=-0.1, maxAngle=0.1)
#array = pipe_methods.imSplit(ext=1, infile=objFile, vertSegments=5, horSegments=4)
#pipe_methods.shiftAllExtensions(infile=objFile, skyfile=skyFile, badpmap=badPMap, \
#  	 vertSegments=1, horSegments=1, interp_type='poly3', \
#  	 stepsize=0.1, xmin=-0.1, xmax=0.1, ymin=-0.1, ymax=0.1)


#pipe_methods.applyShiftAllExtensions(fileList = 'NGC55_14-9-2014_fileNames_short.txt', badpmap=badPMap, \
#  	 vertSegments=1, horSegments=1, interp_type='spline3', \
#  	 stepsize=0.01, xmin=-0.1, xmax=0.1, ymin=-0.1, ymax=0.1)

#pipe_methods.applySubtraction(h_names)
#pipe_methods.extensionMedians(newFile1)
#pipe_methods.extensionMedians(newFile2)
#pipe_methods.extensionMedians(newFile3)

#Cross correlation test - why is it different? 

#pipe_methods.maskFilelcal(objFile, lcalMap14)
#pipe_methods.maskFile(objFile, badPMap14)
#pipe_methods.crossCorr(ext=1, objFile=newFile1, skyFile=newFile2, y1=500, y2=1500, x1=500, x2=1500)
#pipe_methods.crossCorrOne(ext=1, objFile=newFile3, skyFile=newFile2, y1=500, y2=1500, x1=500, x2=1500)

#pipe_methods.minimiseRho(objFile, skyFile, badPMap, interp_type='spline3')

#pipe_methods.shiftAllExtensionsMin(objFile, skyFile, badPMap, vertSegments=2, horSegments=2, interp_type='spline3')
#pipe_methods.subFrames(hobjframe, hskyframe)

#data, coords = pipe_methods.shiftImageSegmentsMin(ext=1, infile=objFile, skyfile=skyFile, badpmap=badPMap14,\
#  	 vertSegments=1, horSegments=1, interp_type='spline3')
#np.savetxt('testCoords.txt', coords, fmt='%10.5f')

#objhdu = fits.PrimaryHDU()
#objhdu.writeto('test.fits', clobber=True)
#fits.append('test.fits', data=data)

#pipe_methods.applyShiftAllExtensionsMin(fileList=names_14, badpmap=badPMap15,\
#  	 vertSegments=1, horSegments=1, interp_type='spline3')

####
##Routine to look at the pixels contaminated by OH emission
#flux = sky_cube.centralSpec()
##Check for where the flux exceeds a certain number of counts 
#indices = np.where(flux > 500)
##Find the sky values at these pixels
#values = flux[indices] 
##Now Loop round a list of combined frames, create a cube, 
##extract the 1D spectrum, find the values at the index points, 
##find the mean difference between these and the sky points 
##and add this to a dictionary labelled by IFU number 
#d = {}
#e = {}
#for fileName in namesOfFile:
#	tempCube = cubeOps(fileName)
#	tempFlux = tempCube.specPlot(3)
#	tempValues = tempFlux[indices]
#	#Array of the differences, take absolute values 
#	diff = abs(values - tempValues)
#	#print diff
#	#Find the average
#	meanDiff = np.median(diff)
#	#Either use the values themselves or the difference
#	d[tempCube.IFUNR] = np.median(tempValues)
#	#d[tempCube.IFUNR] = meanDiff
#for fileName in namesOfFile_1:
#	tempCube = cubeOps(fileName)
#	tempFlux = tempCube.specPlot(3)
#	tempValues = tempFlux[indices]
#	#Array of the differences, take absolute values 
#	diff = abs(values - tempValues)
#	#print diff
#	#Find the average
#	meanDiff = np.median(diff)
#	#Either use the values themselves or the difference
#	e[tempCube.IFUNR] = np.median(tempValues)
#	#d[tempCube.IFUNR] = meanDiff
##Now have a dictionary with the IFU NR and a badness indicator. Print this.
#xAxis = d.keys()
#yAxis = d.values()
#xAxis_1 = e.keys()
#yAxis_1 = e.values()
#fig, ax = plt.subplots(1, 1, figsize=(12.0,12.0))
#ax.plot(xAxis, yAxis)
#ax.plot(xAxis_1, yAxis_1)
#ax.set_title('Sky Tweak Performance vs. IFU ID')
#ax.set_xlabel('Detector ID')
#ax.set_xticks((np.arange(min(xAxis), max(xAxis)+1, 1.0)))
#ax.grid(b=True, which='both', linestyle='--')
#fig.savefig('/Users/owenturner/Documents/PhD/KMOS/KMOS_DATA/Pipeline_Execution/10-4-15_Pairs_14/all_but_9_ws/median_comparison.png')
#plt.show()
#plt.close('all')
#print d.keys()
##print d.values()	










