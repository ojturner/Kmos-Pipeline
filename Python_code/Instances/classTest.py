import os
import sys
# add the class file to the PYTHONPATH
sys.path.append('/disk1/turner/PhD/KMOS/Analysis_Pipeline/Python_code/Class')

# import the relevant modules

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

# import the classes
from pipelineClass import pipelineOps
from cubeClass import cubeOps
from galPhysClass import galPhys
from vel_field_class import vel_field






# Create an instance of the class
pipe_methods = pipelineOps()
cube = cubeOps('/disk1/turner/DATA/GOODS_K_P1_comb_0.8_10/Science/combine_sci_reconstructed_bs008543.fits')
field_instance = vel_field('/disk1/turner/DATA/SSA_HK_P1_comb_0.8_10/Science/combine_sci_reconstructed_n_c49_velocity_map.fits', 14, 15)
#galaxy = galPhys('/disk1/turner/DATA/Gals2/comb/Science/comp_spectrum.fits', 0)
#sky_cube = cubeOps(kskyCube)


##pipe_methods.computeOffsetSegments(objFile, skyFile, badPMap, lcalMap)
##The following are examples of using the functions within the class
##pipe_methods.computeOffsetTopFour('KMOS_SPEC_OBS258_0001_m_2_raw.fits', objFile)
#k_names = '/disk1/turner/PhD/KMOS/KMOS_DATA/NGCLEE/K-band/Calibrations/corrected_object_names.txt'
#h_names = '/disk1/turner/PhD/KMOS/KMOS_DATA/NGCLEE/H-band/Calibrations/corrected_object_names.txt'#
#

##Permanent Names
#raw_14 = '/disk1/turner/PhD/KMOS/KMOS_DATA/NGC55/14-9-2014/NGC55_14_Names.txt'
#names_14 = '/disk1/turner/PhD/KMOS/KMOS_DATA/NGC55/14-9-2014/NGC55_14_Corrected_Names.txt'
#names_14_short = '/disk1/turner/PhD/KMOS/KMOS_DATA/NGC55/14-9-2014/NGC55_14_Corrected_Names_short.txt'
#names_14_shifted = '/disk1/turner/PhD/KMOS/KMOS_DATA/NGC55/14-9-2014/NGC55_14_Corrected_Names_shifted.txt'#

#raw_15 = '/disk1/turner/PhD/KMOS/KMOS_DATA/NGC55/15-09-2014/NGC55_15_Names.txt'
#names_15 = '/disk1/turner/PhD/KMOS/KMOS_DATA/NGC55/15-09-2014/NGC55_15_Corrected_Names.txt'
#names_15_short = '/disk1/turner/PhD/KMOS/KMOS_DATA/NGC55/15-09-2014/NGC55_15_Corrected_Names_short.txt'
#names_15_shifted = '/disk1/turner/PhD/KMOS/KMOS_DATA/NGC55/15-09-2014/NGC55_15_Corrected_Names_shifted.txt'#
#

#noTel_names_15 = '/disk1/turner/PhD/KMOS/KMOS_DATA/NGC55/15-09-2014/NGC55_15_Corrected_Names_noTel.txt'
#noTel_names_14 = '/disk1/turner/PhD/KMOS/KMOS_DATA/NGC55/14-9-2014/NGC55_14_Corrected_Names_noTel.txt'#

##Changes depending on reduction process
#frame_check_names = '/disk1/turner/PhD/KMOS/KMOS_DATA/Pipeline_Execution/frameCheck/combNames.txt'
#sci_names_14 = '/disk1/turner/PhD/KMOS/KMOS_DATA/Pipeline_Execution/S24-3-15/5sig_Science_Output'
#sci_names_15 = '/disk1/turner/PhD/KMOS/KMOS_DATA/Pipeline_Execution/S24-3-15/5sig_Science_Output'
#sci_names_14_noTel = '/disk1/turner/PhD/KMOS/KMOS_DATA/Pipeline_Execution/10-4-15_Pairs_14/all_but_9_ws/sci_names_noTel.txt'
#sci_names_14_noTel_1 = '/disk1/turner/PhD/KMOS/KMOS_DATA/Pipeline_Execution/10-4-15_Pairs_14/all_but_9/sci_names_noTel.txt'
#sci_names_15_noTel = '/disk1/turner/PhD/KMOS/KMOS_DATA/Pipeline_Execution/10-4-15_250-1500_15/Science_Output/sci_names_noTel.txt'#
#

#newFile2 = '/disk1/turner/PhD/KMOS/KMOS_DATA/NGC55/15-09-2014/KMOS_SPEC_OBS259_0014_Corrected_22_spline3_Shifted.fits'
#newFile3 = '/disk1/turner/PhD/KMOS/KMOS_DATA/NGC55/14-9-2014/KMOS_SPEC_OBS258_0009_Corrected_Subtracted.fits'#

##namesOfFile = np.genfromtxt(sci_names_14_noTel, dtype='str')
##namesOfFile_1 = np.genfromtxt(sci_names_14_noTel_1, dtype='str')#

#objCube = '/disk1/turner/PhD/KMOS/KMOS_DATA/Pipeline_Execution/10-4-15_Pairs_14/all_but_9/sci_combined_n55_19__skytweak.fits'
#objSpec2 = '/disk1/turner/PhD/KMOS/KMOS_DATA/Pipeline_Execution/10-4-15_Pairs_14/all_but_9/sci_combined_n55_19__skytweak_spectrum.fits'
#objSpec = '/disk1/turner/PhD/KMOS/KMOS_DATA/NGCLEE/H-band/Science/Best_sci_combined_P108__skytweak_spectrum.fits'
#skySpec = '/disk1/turner/PhD/KMOS/KMOS_DATA/NGCLEE/H-band/Science/cubesky_spectrum.fits'#

#hobjframe = '/disk1/turner/PhD/KMOS/KMOS_DATA/NGCLEE/H-band/raw_frames/KMOS.2014-08-03T00:05:24.218_Corrected_11_spline3_Shifted.fits'
#hskyframe = '/disk1/turner/PhD/KMOS/KMOS_DATA/NGCLEE/H-band/raw_frames/KMOS.2014-08-03T00:03:33.904.fits'

combine_input = '/disk1/turner/DATA/GOODS_K_P2_comb_calibrated_2E17/Science/all_nights.txt'
sci_dir = '/disk1/turner/DATA/GOODS_K_P1_comb_calibrated_1.5E17/Science'
combNames = '/disk1/turner/PhD/KMOS/Analysis_Pipeline/Python_code/Instances/gals_names.txt'
obj_names = '/disk1/turner/DATA/NGC55/YJ/Calibrations/shifted_object_names.txt'

cal_dir = '/disk1/turner/DATA/esther_tester/Calibrations'
gal_dir = '/disk1/turner/DATA/SSA_HK_P2_comb_0.8_15/Science'
sky_cube_gal = gal_dir + '/combine_sci_reconstructed_arm3_sky.fits'
obj_cube_gal = gal_dir + '/combine_sci_reconstructed_s_sa22b-c10.fits'
std_cube_gal = gal_dir + '/combine_sci_reconstructed_c_stars_7656.fits'
raw_file = '/disk1/turner/DATA/Gals1/K/obs_09/raw_frames/Corrected/KMOS_SPEC_OBS344_0018_Corrected.fits'
badpixel_dark_new = '/disk1/turner/DATA/NGC55/15_2014/Calibrations/badpixel_dark_Added.fits'
object_spectrum = gal_dir + '/combine_sci_reconstructed_n_c47_spectrum.fits'
vor_infile = '/disk1/turner/PhD/KMOS/Analysis_Pipeline/Python_code/voronoi/kmos_voronoi_test.txt'
vor_output = '/disk1/turner/PhD/KMOS/Analysis_Pipeline/Python_code/voronoi/voronoi_2d_binning_output.txt'


#pipe_methods.stott_postage_stamps('/disk1/turner/DATA/all_names_calibrated.txt', 'oiii', 5)
cube = cubeOps(obj_cube_gal)
cube.stott_velocity_field('oiii', 3.28733447279, 10, 15, 18, method='median')
#pipe_methods.multi_apply_voronoi_binning('/disk1/turner/DATA/all_names_calibrated.txt', 10.0)
#pipe_methods.combine_by_name(sci_dir, combine_input, 0, 0.8, 2E-17)
#pipe_methods.voronoi_binning_by_line('oiii', obj_cube_gal, 3.47328838, 2.0, '/disk1/turner/PhD/KMOS/Analysis_Pipeline/Python_code/voronoi')
#pipe_methods.apply_voronoi_binning(vor_infile, '/disk1/turner/PhD/KMOS/Analysis_Pipeline/Python_code/voronoi', 1.5)
#hb_flux, hb_vel, hb_sig = pipe_methods.vor_output_fitting(gal_dir, vor_output, obj_cube_gal, std_cube_gal, sky_cube_gal, 19, 20, 3.47328838, stack='median', line='hb')
#oiii_flux, oiii_vel, oiii_sig = pipe_methods.vor_output_fitting(gal_dir, vor_output, obj_cube_gal, std_cube_gal, sky_cube_gal, 19, 20, 3.47328838, stack='median', line='oiii')
#hb_met = pipe_methods.hb_metallicity(oiii_flux, hb_flux)
#print np.nansum(oiii_flux)
#print np.nansum(oiii_flux[15:25, 10:20])
#pipe_methods.telluric_correct('IZ', cal_dir)
#field_instance.fit_kinematic_pa(plot=True, debug=False)
#pipe_methods.multi_plot_all_maps('/disk1/turner/DATA/SSA_HK_P2_comb_calibrated_1E16/all_names.txt', binning=False, xbin=1, ybin=1, interp='sum')
#pipe_methods.multi_plot_OIII_vel_map('/disk1/turner/DATA/GOODS_K_P1_comb_0.8_10/Science/goods_k_p1_spec_data.txt')
#pipe_methods.multi_plot_OIII_vel_map('/disk1/turner/DATA/SSA_HK_P1_comb_0.8_10/Science/ssa22_p1_spec_data.txt')
#pipe_methods.multi_plot_OIII_vel_map('/disk1/turner/DATA/GOODS_K_P2_comb_0.8_10/Science/goods_k_p2_spec_data.txt')
#pipe_methods.multi_plot_OIII_vel_map('/disk1/turner/DATA/SSA_HK_P2_comb_0.8_15/Science/ssa22_p2_spec_data.txt')
# data = cube.spaxel_binning(cube.data, 2, 2)
# print 'this is the original shape: (%s, %s)' % (cube.data.shape[1], cube.data.shape[2])
# print 'this is the new shape: (%s, %s)' % (data.shape[1], data.shape[2])
# print data[400:500,12,12]
# print cube.data
#spectrum, wl = pipe_methods.galExtract(gal_dir, std_cube_gal, obj_cube_gal, sky_cube_gal, 22, 15, 1)
#oiii_values, hb_values = pipe_methods.fit_lines_K(spectrum, wl, redshift=3.4737420)
#param_dict = {'centre': 2.24035, 'sigma': 0.00045, 'amplitude': 0.001}
#cube.OIII_vel_map(redshift=3.47374201278, 
#                  binning=True, 
#                  xbin=1, 
#                  ybin=1, 
#                  interp='mean', 
#                  savefig=True,
#                  **oiii_values)
#pipe_methods.multi_plot_K_image('/disk1/turner/DATA/GOODS_K_P2_comb_0.8_10/Science/goods_k_p2_spec_data.txt')
#pipe_methods.multi_plot_HK_image('/disk1/turner/DATA/SSA_HK_P2_comb_0.8_15/Science/ssa22_p2_spec_data.txt')
#pipe_methods.seeing_better_than(combine_input, 0.65)
#pipe_methods.av_seeing(combine_input)
#pipe_methods.singlePixelExtractMulti_OIII(sci_dir + '/sn.txt', sci_dir)
#pipe_methods.singlePixelExtract_OIII5008(sci_dir, obj_cube_gal, 17, 17, 3.2943, 1)
#pipe_methods.singlePixelExtract_OIII4960(sci_dir, obj_cube_gal, 23, 22, 3.08740, 1)
#pipe_methods.singlePixelExtract_Hb(sci_dir, obj_cube_gal, 23, 22, 3.08740, 1)
#pipe_methods.stackSpectra('/disk1/turner/DATA/combined_spectra.txt', 0.00028076)
#pipe_methods.pSTNK(object_spectrum, 3.47465)
#pipe_methods.maskExtraPixels(raw_file, badpixel_dark_new)
#galaxy.printProps()
#galaxy.plotSpec()
#galaxy.fitHbandOIII()
#pipe_methods.plotHandOII(1.666)
#pipe_methods.galExtract(gal_dir, std_cube_gal, obj_cube_gal, sky_cube_gal, 9, 17, 1)
#pipe_methods.multiGalExtract('/disk1/turner/DATA/all_names.txt', 1)
#pipe_methods.frameCheck(sci_dir, obj_names, 'n55_19')
#pipe_methods.compareSky(sci_dir, combNames)
#new_Table = pipe_methods.reduce_list_seeing(combine_input, 0.5, 1.0)
#print 'length new_Table is: %s'  % len(new_Table)
#name_Table = pipe_methods.reduce_list_name(new_Table, 'P107')
#print 'length name_Table is: %s'  % len(name_Table)
#combine_Table = pipe_methods.reduce_list_sky(name_Table, 1.2)
#print 'length combine_Table is: %s'  % len(combine_Table)
#stuff, stuff1, stuff2 = pipe_methods.compareSky(sci_dir='/disk1/turner/DATA/688/H/Science', combNames='co_names.txt')

#pipe_methods.saveSpec('/disk1/turner/PhD/KMOS/KMOS_DATA/Pipeline_Execution/May_20th_tests/corr_with/sci_combined_s152__skytweak.fits')
#pipe_methods.plotSpecs(gal_dir, object_spectrum, sky_cube_gal, 1)

#pipe_methods.multiExtractSpec(skyCube=hskyCube, frameNames=h_names)
#Now try the subtraction method 
#pipe_methods.subFrames('/disk1/turner/PhD/KMOS/KMOS_DATA/NGC55/14-9-2014/KMOS_SPEC_OBS258_0009_Corrected_11_spline3_Shifted.fits', skyFile)
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
#fig.savefig('/disk1/turner/PhD/KMOS/KMOS_DATA/Pipeline_Execution/10-4-15_Pairs_14/all_but_9_ws/median_comparison.png')
#plt.show()
#plt.close('all')
#print d.keys()
##print d.values()	









