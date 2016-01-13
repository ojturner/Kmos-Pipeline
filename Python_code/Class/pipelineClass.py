# This class houses the methods which are relevant to manual additions to the
# ESO KMOS pipeline Mainly concentrating on two procedures - 1) pedestal
# readout column correction and 2) shifting and aligning sky and object
# images before data cube reconstruction


# import the relevant modules
import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pyraf
import numpy.ma as ma
from scipy.spatial import distance
from scipy import ndimage
from copy import copy
from lmfit import Model
from itertools import cycle as cycle
from lmfit.models import GaussianModel, PolynomialModel
from scipy.optimize import minimize
from astropy.io import fits, ascii
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy import poly1d

# add the class file to the PYTHONPATH
sys.path.append('/disk1/turner/PhD'
                + '/KMOS/Analysis_Pipeline/Python_code/Class')

from cubeClass import cubeOps
from galPhysClass import galPhys


class pipelineOps(object):

    # Initialiser creates an instance of the spectrumFit object
    def __init__(self):

        self.self = self

    def compute_offset_top_four(self,
                                rawSubFile,
                                objectFile):

        """
        Def:
        Take the object image after calibration,
        along with a sky frame, bad pixel frame and
        the lcal frame and homogenise the readout columns,
        so that after subtraction gives 0 median noise

        Input:
                objectFile - object image to be corrected
                rawSubFile - raw subtracted image


        Output:
                newObjData - corrected 2D data array

        """

        corrected_extensions = []

        # Read in the tables of data
        table_o = fits.open(objectFile)

        fits_header = table_o[0].header

        header_one = table_o[1].header

        header_two = table_o[2].header

        header_three = table_o[3].header

        print fits_header
        print header_one
        print header_two
        print header_three

        table_s = fits.open(rawSubFile)

        # Loop over the fits image extensions, do the same each time
        for count in range(1, 4):

            data_o = table_o[count].data

            data_s = table_s[count].data

            # Counters for the slicing
            x = 0

            y = 64

            # 1D arrays to host the data
            sub_array = []

            obj_array = []

            for i in range(32):

                # Slice each of the data files into
                # 32 columns of 64 pixels width
                newObjectArray = data_o[:, x:y]

                new_sub_array = data_s[2044:2048, x:y]

                obj_array.append(newObjectArray)

                sub_array.append(new_sub_array)

                # Add 64 to the counters each time to create the slices
                x += 64
                y += 64

            # testData = np.hstack(sub_array)
            # fileName = 'testTopFour' + str(count) + '.fits'
            # fits.writeto(fileName, data=testData, clobber=True)

            for num in range(len(obj_array)):

                correctionMedian = np.median(sub_array[num])
                obj_array[num] -= correctionMedian
                print correctionMedian
                print sub_array[num][:, 0:1]

            # Have now made the correction to all
            # 32 of the 64 pixel width columns.
            # All that is left to do is stitch
            # the obj_array arrays back together to
            # give a single 2048x2048 array.

            newObjData = np.hstack(obj_array)

            corrected_extensions.append(newObjData)

            if count == 1:

                print 'Computed First Correction'

            elif count == 2:

                print 'Computed Second Correction'

            elif count == 3:

                print 'Computed Third Correction'

        # Create the object fits file with the
        # three corrected extensions

        fileName = raw_input('Enter a name for the'
                             + ' corrected fits file: ') + '.fits'

        #  Note that the readout complained about the header not being
        #  in the correct fits format

        hdu = fits.PrimaryHDU(header=fits_header)

        hdu.writeto(fileName,
                    clobber=True)

        fits.append(fileName,
                    data=corrected_extensions[0],
                    header=header_one)

        fits.append(fileName,
                    data=corrected_extensions[1],
                    header=header_two)

        fits.append(fileName,
                    data=corrected_extensions[2],
                    header=header_three)

    def computeOffsetSegments(self,
                              objectFile,
                              skyFile,
                              badPMap,
                              lcalMap):
        """
        Def:
        Take the object image after calibration,
        along with a sky frame, bad pixel frame and
        the lcal frame and homogenise the readout columns,
        so that after subtraction gives 0 median noise

        Input:
                objectFile - object image to be corrected
                skyFile - the corresponding sky file
                badPMap - the bad pixel frame generated by the pipeline
                lcalMap - wavelength calibration frame from the pipeline

        Output:
                objectFile_Corrected.fits - corrected 2D data array

        """
        # function should be identical to compute
        # Offset until the initial pixel loop

        # Set up vector to house the corrected extensions
        correctedExtensions = []

        # Set up vector to house the segments
        # to be vStacked. This is differnet from the
        # usual compute offset method which just
        # uses each individual readout column

        # Read in the tables of data
        table_o = fits.open(objectFile)

        fits_header = table_o[0].header

        header_one = table_o[1].header

        header_two = table_o[2].header

        header_three = table_o[3].header

        temp = sys.stdout

        sys.stdout = open('log.txt', 'w')

        print fits_header
        print header_one
        print header_two
        print header_three

        sys.stdout.close()

        sys.stdout = temp

        os.system('rm log.txt')

        table_s = fits.open(skyFile)

        bad_pixel_table = fits.open(badPMap)

        # Now choose the correct rotation angle
        lcal_table = fits.open(lcalMap)

        # This is a list of all possible rotation angles
        angleList = np.array([0, 60, 120, 180, 240, 300])

        # Select the ocs.rot.naangle keyword
        obsAngle = table_o[0].header["HIERARCH ESO OCS ROT NAANGLE"]

        print obsAngle

        if obsAngle < 0:

            obsAngle = obsAngle + 360

        # Find where the difference between
        # the observed and idealised angle is minimum
        newAngleList = abs(obsAngle - angleList)

        n = newAngleList.argmin()

        obsAngleNew = angleList[n]

        print obsAngleNew

        # Find the extension to which this corresponds
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

        # Loop over the fits image extensions, do the same each time
        for count in range(1, 4):

            vStackArray = []

            print val

            data_o = table_o[count].data

            data_s = table_s[count].data

            # Create copies of the data arrays
            # so that I can mask the bad pixels
            # and lcal pixels outside of the loop,
            # instead of wasting time inside

            manObjData = copy(data_o)

            manSkyData = copy(data_s)

            # Read in the bad pixel and lcal maps
            bad_pixel_data = bad_pixel_table[count].data

            lcal_data = lcal_table[val].data

            # Find the coordinates of the bad pixels and the slitlets
            bad_pixel_coords = np.where(bad_pixel_data == 0)

            lcal_pixel_coords = np.where(lcal_data > 0)

            # Loop around the bad pixel locations and
            # mask off on the manObjData and manSkyData
            for i in range(len(bad_pixel_coords[0])):

                # Because of the way np.where works,
                # need to define the x and y coords in this way
                xcoord = bad_pixel_coords[0][i]

                ycoord = bad_pixel_coords[1][i]

                # Now set all positions where there is
                # a dead pixel to np.nan in the object and sky
                manObjData[xcoord][ycoord] = np.nan

                manSkyData[xcoord][ycoord] = np.nan

            # Loop around the slitlet positions
            for i in range(len(lcal_pixel_coords[0])):

                # Do the same, this time for the slitlet
                # positions (substantially more will have a value)
                xcoord = lcal_pixel_coords[0][i]

                ycoord = lcal_pixel_coords[1][i]

                # Set all of these locations to nan
                manObjData[xcoord][ycoord] = np.nan
                manSkyData[xcoord][ycoord] = np.nan

            # Need to slice this 2D array into a
            # 1D array of 2D arrays, each of which is 64 pixels wide
            # so that they can be examined in turn and loop over them

            # Counters for the horizontal slicing
            hor1 = 0

            hor2 = 128

            for j in range(16):

                print hor1

                print hor2

                # Counters for the slicing vertical slicing
                x = 0
                y = 64

                # 1D arrays to host the data
                skyArray = []

                objArray = []

                manObjArray = []

                manSkyArray = []

                badPArray = []

                lcalArray = []

                testArray = []

                for i in range(32):

                    # Slice each of the data files
                    # into 32 columns of 64 pixels width
                    newObjectArray = data_o[hor1:hor2, x:y]

                    newSkyArray = data_s[hor1:hor2, x:y]

                    newManObjArray = manObjData[hor1:hor2, x:y]

                    newManSkyArray = manSkyData[hor1:hor2, x:y]

                    newPArray = bad_pixel_data[hor1:hor2, x:y]

                    newCalArray = lcal_data[hor1:hor2, x:y]

                    # newTestArray = test_array[:,x:y]
                    # testArray.append(newTestArray)

                    objArray.append(newObjectArray)

                    skyArray.append(newSkyArray)

                    manObjArray.append(newManObjArray)

                    manSkyArray.append(newManSkyArray)

                    badPArray.append(newPArray)

                    lcalArray.append(newCalArray)

                    # Add 64 to the counters each time to create the slices
                    x += 64

                    y += 64

                    # Have sliced each matrix into 2048x64.
                    # All that's left to do is slice
                    # Each of these into 8 chunks of 256x64,
                    # 16 chunks of 128x64 and 32 chunks of 64x64

                print objArray[1].shape

                # Start the loop for all the columns in the Array vectors

                for num in range(len(objArray)):

                    # Now all the pixels that
                    # shouldn't be included in the median
                    # have value nan. Can then just
                    # do np.nanmean(objTemp) which will ignore nan's
                    # then repeat the process for the sky,
                    # compare the mean's, compute and apply the offset.

                    obj_mean = np.nanmedian(manObjArray[num])

                    sky_mean = np.nanmedian(manSkyArray[num])

                    # print sky_mean
                    # print obj_mean

                    # Need to compare the two medians
                    # to see how to apply the offset.
                    # If the sky is brighter, add
                    # the difference to the object image

                    if sky_mean > obj_mean:

                        objArray[num] += abs(sky_mean - obj_mean)

                    elif obj_mean > sky_mean:

                        objArray[num] -= abs(obj_mean - sky_mean)

                # Have now made the correction to all
                # 32 of the 64 pixel width columns.
                # Previously would stack the data here,
                # but the loop goes back to the beginning
                # So need to save to a different object,
                # and then vstack at the end.

                vStackArray.append(np.hstack(objArray))

                hor1 += 128

                hor2 += 128

            # Now just need to vstack all of these
            # arrays and will have a 2048x2048 corrected array

            newObjData = np.vstack(vStackArray)

            print newObjData.shape

            correctedExtensions.append(newObjData)

            if count == 1:

                print 'Computed First Correction'

            elif count == 2:

                print 'Computed Second Correction'

            elif count == 3:

                print 'Computed Third Correction'

            val += 1

        # Create the object fits file with the three corrected extensions

        fileName = objectFile[0:-5] + '_Corrected' + '.fits'

        # Note that the readout complained about the header not being
        # in the correct fits format

        hdu = fits.PrimaryHDU(header=fits_header)

        hdu.writeto(fileName,
                    clobber=True)

        fits.append(fileName,
                    data=correctedExtensions[0],
                    header=header_one)

        fits.append(fileName,
                    data=correctedExtensions[1],
                    header=header_two)

        fits.append(fileName,
                    data=correctedExtensions[2],
                    header=header_three)

    def subFrames(self,
                  objectFile,
                  skyFile):

        """
        Def: Subtract all extensions of a skyfile from
             an object file

        Input:
                objectFile - KMOS object fits file
                skyFile - KMOS object sky file

        Output:
                objectFile_Subtracted.fits
        """

        # Read in the object and sky files
        objData = fits.open(objectFile)

        skyData = fits.open(skyFile)

        # Find the header and extensions of the new fits file
        header = objData[0].header

        headerOne = objData[1].header

        headerTwo = objData[2].header

        headerThree = objData[3].header

        temp = sys.stdout

        sys.stdout = open('log.txt', 'w')

        print header
        print headerOne
        print headerTwo
        print headerThree

        sys.stdout.close()

        sys.stdout = temp

        ext1 = objData[1].data - skyData[1].data

        ext2 = objData[2].data - skyData[2].data

        ext3 = objData[3].data - skyData[3].data

        # Write out to a different fits file, with name user specified

        nameOfFile = objectFile[:-5] + '_Subtracted.fits'

        hdu = fits.PrimaryHDU(header=header)

        hdu.writeto(nameOfFile,
                    clobber=True)

        fits.append(nameOfFile,
                    data=ext1,
                    header=headerOne)

        fits.append(nameOfFile,
                    data=ext2,
                    header=headerTwo)

        fits.append(nameOfFile,
                    data=ext3,
                    header=headerThree)

        os.system('rm log.txt')

    def applySubtraction(self,
                         fileList):

        """
        Def:
        Apply the subframes method to a list of files

        Input:
                fileList - list of files, where the top line is
                            the name and type, the two columns
                             and the name of the file
                              and either O for object and S for sky.

        Output: List of subtracted files saved in the object directory
        """

        # Read in the data from the fileList
        data = np.genfromtxt(fileList, dtype='str')

        # Save the names and types as lists
        names = data[0:, 0]

        # print names
        types = data[0:, 1]

        # Loop round all names and apply the computeOffsetSegments method
        for i in range(1, len(names)):

            if types[i] == 'O':

                objFile = names[i]

                if i == 1:

                    skyFile = names[i + 1]

                elif types[i - 1] == 'S':

                    skyFile = names[i - 1]

                else:

                    skyFile = names[i + 1]

                print '[INFO]: Subbing file: %s : ' % objFile

                # Now use the method defined within this class

                self.subFrames(objFile, skyFile)

    def pixelHistogram(self,
                       subFile,
                       subCorFile,
                       x1,
                       x2):

        """
        Def: Created a histogram of pixel values
        before and after correction.

        Input:
                subFile - raw subtracted file
                subCorFile - subtracted file post correction
                x1 - starting x pixel number
                x2 - finishing x pixel number

        Output:
                Plot showing the histograms


        """

        # First read in the files

        subData = fits.open(subFile)

        subCorData = fits.open(subCorFile)

        # At the moment we'll just consider the first extension for our data

        subData = subData[1].data

        subCorData = subCorData[1].data

        # The input numbers define the left and
        # right edges of the pixel section

        subData = subData[:, x1:x2]

        subCorData = subCorData[:, x1:x2]

        # This gives an array of arrays,
        # we just care about the numbers, not spatial info
        # Use ravel() to convert these into lists for the histogram

        subData = subData.ravel()

        subCorData = subCorData.ravel()

        print len(subData)
        print subData
        print np.median(subData)
        print np.median(subCorData)

        # Create the bins array for both histograms

        bins = np.arange(-15, 15, 1)

        plt.close('all')

        fig, ax = plt.subplots(1, 1, figsize=(12, 12))

        # Plot the histograms

        n1, bins1, patches1 = ax.hist(subData,
                                      bins=bins,
                                      histtype='step',
                                      color='green',
                                      linewidth=3,
                                      label='Before Correction')

        n2, bins2, patches2 = ax.hist(subCorData,
                                      bins=bins,
                                      histtype='step',
                                      color='blue',
                                      linewidth=3,
                                      alpha=0.5,
                                      label='After Correction')

        # Now want to fit the gaussians to the
        # histograms using lmfit gaussian models
        # gaussian model number 1

        mod1 = GaussianModel()

        # Create a new bins vector for the fit

        fitBins = np.arange(-14.5, 14.5, 1)

        print len(fitBins)

        # Take an initial guess at what the model parameters are
        # In this case the gaussian model has three parameters,
        # Which are amplitude, center and sigma

        pars1 = mod1.guess(n1, x=fitBins)

        # Perform the actual fit

        out1 = mod1.fit(n1, pars1, x=fitBins)

        # Now want to add this curve to our plot
        ax.plot(fitBins,
                out1.best_fit,
                linewidth=2.0,
                label='b.c. model',
                color='green')

        # Repeat for the corrected data

        mod2 = GaussianModel()

        pars2 = mod2.guess(n2, x=fitBins)

        out2 = mod2.fit(n2, pars2, x=fitBins)

        ax.plot(fitBins,
                out2.best_fit,
                linewidth=2.0,
                label='a.c. model',
                color='blue')

        ax.set_xlabel('Counts per pixel', fontsize=24)

        ax.set_ylabel('Number per bin', fontsize=24)

        ax.set_title('Sub Frame, Improvement After Correction',
                     fontsize=28)

        ax.tick_params(axis='both',
                       which='major',
                       labelsize=15)

        ax.legend(loc='upper left', fontsize=10)

        fig.savefig(raw_input('Enter the plot name: '))

        plt.show()

        plt.close('all')

        # Print out the fit reports to look
        # at the centre at S.D. of each model

        print out1.fit_report()

        print out2.fit_report()

    def stackLcal(self,
                  lcalFile):
        """
        Def:
        Stack all the lcal files in the different rotator angles
        to get an overall average

        Input:
                lcalXXX.fits produced by the pipeline

        Output:
                lcal1.fits - rotator angle averaged for det. 1
                lcal2.fits - rotator angle averaged for det. 2
                lcal3.fits - rotator angle averaged for det. 3
        """

        # Read in the file

        lcal = fits.open(lcalFile)

        # Want to have it so that only nan values survive. Can do this
        # By substituting 'False' values for
        # the numbers and true values for nan
        # So that when multiplied only True * True
        # all the way will survive

        # Loop round all the extensions and change the np.nan values to True
        # and the pixels with values to False

        d = {}

        for i in range(1, 19):

            data = lcal[i].data

            # Define the coordinates where there is a value

            value_coords = np.where(data > 0)

            # Define where there is np.nan

            nan_coords = np.where(np.isnan(data))

            temp = np.empty(shape=[2048, 2048], dtype=float)

            # loop over the pixel values in data and
            # change to either True or False

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

        # Now create the individual lcal frames
        # by stacking the extensions together
        # This gives a more conservative estimate
        # of the pixels we should and should
        # not be using throughout the data analysis stage

        lcal1 = d[1] * d[4] * d[7] * d[10] * d[13] * d[16]

        lcal2 = d[2] * d[5] * d[8] * d[11] * d[14] * d[17]

        lcal3 = d[3] * d[6] * d[9] * d[12] * d[15] * d[18]

        fits.writeto(filename='lcal1.fits', data=lcal1, clobber=True)

        fits.writeto(filename='lcal2.fits', data=lcal2, clobber=True)

        fits.writeto(filename='lcal3.fits', data=lcal3, clobber=True)

    def applyCorrection(self,
                        fileList,
                        badPMap,
                        lcalMap):

        """
        Def:
        Apply the computeOffsetSegments method to lots of files

        Input:
                fileList - txt file containing the names of files
                            to be corrected
                badPMap - the bad_pixel_Added file produced by
                            the pipeline and addons
                lcalMap - the wavelength calibration image produced
                            by the pipeline

        Output:
                the input list as corrected files _Corrected.fits
        """

        # Read in the data from the fileList

        data = np.genfromtxt(fileList, dtype='str')

        # Save the names and types as lists

        names = data[0:, 0]

        types = data[0:, 1]

        # Loop round all names and apply the computeOffsetSegments method

        for i in range(1, len(names)):

            if types[i] == 'O':

                objFile = names[i]

                if i == 1:

                    skyFile = names[i + 1]

                elif types[i - 1] == 'S':

                    skyFile = names[i - 1]

                else:

                    skyFile = names[i + 1]

                # Now use the method defined within this class

                self.computeOffsetSegments(objFile,
                                           skyFile,
                                           badPMap,
                                           lcalMap)

    def rowMedian(self,
                  subFile,
                  y1,
                  y2,
                  x1,
                  x2):

        """
        Def:
        Create a histogram of pixel values on the
        subtracted frame before and after correction

        Input:
                subFile - any subtracted image
                y1, y2, x1, x2 - the x and y pixel ranges for median

        Output:
                return the median array of pixel values
        """
        # First read in the files

        subData = fits.open(subFile)

        # At the moment we'll just consider the first extension for our data

        subData = subData[1].data

        # The input numbers define the left
        # and right edges of the pixel section

        subData = subData[y1:y2, x1:x2]

        # What we have now is an array
        # of 500 entries, each of which contains 500 entries.
        # These are the columns of the fits file, the pixel count values
        # Now need to compute the median
        # of each and store in an array and return

        medArray = np.median(subData, axis=0)

        return medArray

    def plotMedian(self,
                   rawSubFile,
                   subFile,
                   segmentsSubFile,
                   top4SubFile,
                   y1,
                   y2,
                   x1,
                   x2):

        """
        Def:
        Plots the median values of different subtracted frames
        together on the same axes. Gives an indication for which
        column correction method produces the smoothest results

        Input:
                rawSubFile - the raw subtracted file pre-correction
                subFile - the subtracte file post-correction
                segmentssubfile - subtracted using segments method
                top4subfile - subtracted using segments method
                x1, x2, y1, y2 - the x and y pixel ranges

        Output:
                plot of the comparison

        """

        # Find the medians of the different subtracted frames and plot

        rawMed = self.rowMedian(rawSubFile,
                                y1,
                                y2,
                                x1,
                                x2)

        normMed = self.rowMedian(subFile,
                                 y1,
                                 y2,
                                 x1,
                                 x2)

        segmentsMed = self.rowMedian(segmentsSubFile,
                                     y1,
                                     y2,
                                     x1,
                                     x2)

        top4Med = self.rowMedian(top4SubFile,
                                 y1,
                                 y2,
                                 x1,
                                 x2)

        # Generate the x-axis vector, which is just the number of pixels

        xVec = range(len(rawMed))

        xVec = np.array(xVec)

        # Plot everything against this xVec in turn

        fig, ax = plt.subplots(1, 1, figsize=(18, 10))

        ax.plot(xVec, segmentsMed, label='Segments')

        ax.plot(xVec, normMed, label='Cor')

        ax.plot(xVec, top4Med, label='top')

        ax.plot(xVec, rawMed, label='raw')

        ax.set_xlabel('Pixel Number', fontsize=24)

        ax.set_ylabel('Median', fontsize=24)

        ax.set_title('Pixel Medians for Different Methods', fontsize=30)

        ax.tick_params(axis='both', which='major', labelsize=15)

        ax.legend(loc='upper left', fontsize=10)

        fig.savefig(raw_input('Enter the File name: ') + '.png')

        plt.show()

        plt.close('all')

    def badPixelextend(self,
                       badpmap):

        """
        Def:
        Take an arbitary bad pixel mask and mask
        off the pixels in a ring around
        the original bad pixel location.

        Inputs:
                badpmap - any bad pixel mask with 0 as bad and 1 as good

        Outputs:
                badpmap_Added.fits - remasked bad pixel

        """

        # Read in data

        badpTable = fits.open(badpmap)

        extArray = []

        primHeader = badpTable[0].header

        ext1Header = badpTable[1].header

        ext2Header = badpTable[2].header

        ext3Header = badpTable[3].header

        temp = sys.stdout

        sys.stdout = open('log.txt', 'w')

        print primHeader
        print ext1Header
        print ext2Header
        print ext3Header

        sys.stdout.close()

        sys.stdout = temp

        # Loop around the extensions and make ammendments

        for i in range(1, 4):

            badpData = badpTable[i].data

            # Now have the 2048x2048 data array.
            # Want to find all the 0 values
            # And make the points surrounding this zero also

            badpCoords = np.where(badpData == 0)

            # Loop around the bad pixel locations and mask off

            for i in range(len(badpCoords[0])):

                # Because of the way np.where works,
                # need to define the x and y coords in this way

                if (badpCoords[0][i] < 2047) and (badpCoords[1][i] < 2047):

                    xcoord = badpCoords[0][i]

                    ycoord = badpCoords[1][i]

                    # Now set all positions where there
                    # is a dead pixel to np.nan in the object and sky

                    badpData[xcoord][ycoord + 1] = 0

                    badpData[xcoord][ycoord - 1] = 0

                    badpData[xcoord + 1][ycoord] = 0

                    badpData[xcoord - 1][ycoord] = 0

            extArray.append(badpData)

        # Now write out to a new fits file

        badpName = badpmap[:-5] + '_Added.fits'

        hdu = fits.PrimaryHDU(header=primHeader)

        hdu.writeto(badpName,
                    clobber=True)

        fits.append(badpName,
                    data=extArray[0],
                    header=ext1Header)

        fits.append(badpName,
                    data=extArray[1],
                    header=ext2Header)

        fits.append(badpName,
                    data=extArray[2],
                    header=ext3Header)

        os.system('rm log.txt')

    def extensionMedians(self,
                         subbedFile):

        """
        Def:
        Take the masked subtracted file and compute
        the median and standard deviation of each extension
        then print out these values to the terminal

        Input:
                subbedFile - sky subtracted fits image

        Output:
                print of the SD and Median

        """

        dataTable = fits.open(subbedFile)

        for i in range(1, 4):

            data = dataTable[i].data

            med = np.nanmedian(data)

            st = np.nanstd(data)

            print 'The median of extension %s is: %s \n' % (i,
                                                            med) \
                + ' The standard Deviation of extension %s is: %s' % (i,
                                                                      st)

    def maskFile(self,
                 inFile,
                 badpFile):

        """
        Def:
        Take the badPixelMap and apply it to any inFile
        to mask off the bad pixels. Particularly useful
        for applying to skyFiles before subtraction

        Input:
                inFile - KMOS object fits file
                badpFile - bad pixel mask from pipeline

        """

        dataTable = fits.open(inFile)

        badPTable = fits.open(badpFile)

        extArray = []

        primHeader = dataTable[0].header

        ext1Header = dataTable[1].header

        ext2Header = dataTable[2].header

        ext3Header = dataTable[3].header

        temp = sys.stdout

        sys.stdout = open('log.txt', 'w')

        print primHeader
        print ext1Header
        print ext2Header
        print ext3Header

        sys.stdout.close()

        sys.stdout = temp

        for i in range(1, 4):

            data = dataTable[i].data

            badpData = badPTable[i].data

            # Now find the bad pixel locations
            # and mask off the data appropriately

            bad_pixel_coords = np.where(badpData == 0)

            # Loop around the bad pixel locations
            # and mask off on the manObjData and manSkyData

            for i in range(len(bad_pixel_coords[0])):

                # Because of the way np.where works,
                # need to define the x and y coords in this way

                xcoord = bad_pixel_coords[0][i]

                ycoord = bad_pixel_coords[1][i]

                # Now set all positions where there
                # is a dead pixel to np.nan in the object and sky

                data[xcoord][ycoord] = np.nan

            extArray.append(data)

        # Write out the new data

        fileName = inFile[:-5] + '_masked.fits'

        hdu = fits.PrimaryHDU(header=primHeader)

        hdu.writeto(fileName,
                    clobber=True)

        fits.append(fileName,
                    data=extArray[0],
                    header=ext1Header)

        fits.append(fileName,
                    data=extArray[1],
                    header=ext2Header)

        fits.append(fileName,
                    data=extArray[2],
                    header=ext3Header)

        os.system('rm log.txt')

    def maskFilelcal(self, inFile, lcalFile):

        """
        Def:
        Take the lcalFile and apply it to any inFile
        to mask off the bad pixels. Particularly useful
        for applying to skyFiles before subtraction

        Input:
                 lcalFile - wavelength cal image from pipeline 
        """

        dataTable = fits.open(inFile)

        badPTable = fits.open(lcalFile)

        extArray = []

        primHeader = dataTable[0].header

        ext1Header = dataTable[1].header

        ext2Header = dataTable[2].header

        ext3Header = dataTable[3].header

        temp = sys.stdout

        sys.stdout = open('log.txt', 'w')

        print primHeader
        print ext1Header
        print ext2Header
        print ext3Header

        sys.stdout.close()

        sys.stdout = temp

        for i in range(1, 4):

            data = dataTable[i].data

            badpData = badPTable[i].data

            # Now find the bad pixel locations
            # and mask off the data appropriately

            bad_pixel_coords = np.where(np.isnan(badpData))

            # Loop around the bad pixel locations and mask

            # off on the manObjData and manSkyData

            for i in range(len(bad_pixel_coords[0])):

                # Because of the way np.where works, need
                # to define the x and y coords in this way

                xcoord = bad_pixel_coords[0][i]

                ycoord = bad_pixel_coords[1][i]

                # Now set all positions where there is
                # a dead pixel to np.nan in the object and sky

                data[xcoord][ycoord] = np.nan

            extArray.append(data)

        # Write out the new data

        fileName = inFile[:-5] + '_masked_lcal.fits'

        hdu = fits.PrimaryHDU(header=primHeader)

        hdu.writeto(fileName,
                    clobber=True)

        fits.append(fileName,
                    data=extArray[0],
                    header=ext1Header)

        fits.append(fileName,
                    data=extArray[1],
                    header=ext2Header)

        fits.append(fileName, data=extArray[2], header=ext3Header)

        os.system('rm log.txt')

    def crossCorr(self,
                  ext,
                  objFile,
                  skyFile,
                  y1,
                  y2,
                  x1,
                  x2):

        """
        Def:
        Compute the cross-correlation coefficient for a given object
        and skyfile. Define the pixel range over which to compute the
        coefficient.

        Inputs:
                ext - detector extension, must be either 1, 2, 3
                objFile - Input object file to compute correlation
                skyFile - sky image to compare objFile with
                y1, y2, x1, x2 - the range to compute rho over

        Output: 
                rho - the cross correlation coefficient
        """

        # Trying to compute the similarity
        # between a square grid of pixels
        # from the object file and from the skyfile.
        # Do this using the correlation coeff.

        # First read in the full 2048x2048 data
        # arrays from the object and sky files
        # Will first consider just the first detector
        # and can expand on this later

        objData = fits.open(objFile)

        objData = objData[ext].data

        skyData = fits.open(skyFile)

        skyData = skyData[ext].data

        objData = np.array(objData[y1:y2, x1:x2])

        skyData = np.array(skyData[y1:y2, x1:x2])

        # first mask off the pixels with 0 value
        # from both the object and the sky

        objDataMedian = np.nanmedian(objData)

        skyDataMedian = np.nanmedian(skyData)

        skyDataStd = np.nanstd(skyData)

        objDataStd = np.nanstd(objData)

        # MASKING THE HIGH SIGMA PIXELS FOR BETTER RHO
        # Let's try masking the pixels which are
        # bigger than 1000 counts and less than 50

        objData[objData < 500] = np.nan

        objData[objData > 10000] = np.nan

        skyData[skyData < 500] = np.nan

        skyData[skyData > 10000] = np.nan

        newobjDataMedian = np.nanmedian(objData)

        newskyDataMedian = np.nanmedian(skyData)

        firstPart = np.sqrt(np.nansum((objData - newobjDataMedian) ** 2))

        secondPart = np.sqrt(np.nansum((skyData - newskyDataMedian) ** 2))

        denom = firstPart * secondPart

        # print denom

        numer = np.nansum((objData - newobjDataMedian) *
                          (skyData - newskyDataMedian))

        rho = numer / denom

        # print rho

        return rho

    def crossCorrZeroth(self,
                        objFile,
                        skyFile,
                        y1,
                        y2,
                        x1,
                        x2):

        """
        Def:
        Compute the cross-correlation coefficient for a given object
        and skyfile. Define the pixel range over which to compute the
        coefficient.

        Inputs:
                ext - detector extension, must be either 1, 2, 3
                objFile - Input object file to compute correlation
                skyFile - sky image to compare objFile with
                y1, y2, x1, x2 - the range to compute rho over

        Output:
                rho - the cross correlation coefficient
        """

        # Trying to compute the similarity
        # between a square grid of pixels
        # from the object file and from the skyfile.
        # Do this using the correlation coeff.
        # First read in the full 2048x2048
        # data arrays from the object and sky files
        # Will first consider just the first
        # detector and can expand on this later

        objData = fits.open(objFile)

        objData = objData[0].data

        skyData = fits.open(skyFile)

        skyData = skyData[1].data

        objData = np.array(objData[y1:y2, x1:x2])

        skyData = np.array(skyData[y1:y2, x1:x2])

        # FIRST MASK OFF THE PIXELS WITH 0 VALUE FROM BOTH OBJECT AND SKY

        objDataMedian = np.nanmedian(objData)

        skyDataMedian = np.nanmedian(skyData)

        skyDataStd = np.nanstd(skyData)

        objDataStd = np.nanstd(objData)

        # MASKING THE HIGH SIGMA PIXELS FOR BETTER RHO
        # Let's try masking the pixels which are bigger than 1000
        # counts and less than 50

        # NOTE FOR FUTURE - THIS IS NOT A TOTALLY SECURE WAY OF
        # COMPUTING THE CROSS CORR.

        objData[objData < 500] = np.nan

        objData[objData > 10000] = np.nan

        skyData[skyData < 500] = np.nan

        skyData[skyData > 10000] = np.nan

        newobjDataMedian = np.nanmedian(objData)

        newskyDataMedian = np.nanmedian(skyData)

        firstPart = np.sqrt(np.nansum((objData - newobjDataMedian) ** 2))

        secondPart = np.sqrt(np.nansum((skyData - newskyDataMedian) ** 2))

        denom = firstPart * secondPart

        # print denom

        numer = np.nansum((objData - newobjDataMedian) *
                          (skyData - newskyDataMedian))

        rho = numer / denom

        # print rho

        return rho

    def crossCorrFirst(self, objFile, skyFile, y1, y2, x1, x2):

        """
        Def:
        Compute the cross-correlation coefficient for a given object
        and skyfile. Define the pixel range over which to compute the
        coefficient.

        Inputs:
                ext - detector extension, must be either 1, 2, 3
                objFile - Input object file to compute correlation
                skyFile - sky image to compare objFile with
                y1, y2, x1, x2 - the range to compute rho over

        Output:
                rho - the cross correlation coefficient
        """

        # Trying to compute the similarity between a square grid of pixels
        # from the object file and from the skyfile.
        # Do this using the correlation coeff.
        # First read in the full 2048x2048 data
        # arrays from the object and sky files
        # Will first consider just the first
        # detector and can expand on this later

        objData = fits.open(objFile)

        objData = objData[1].data

        skyData = fits.open(skyFile)

        skyData = skyData[1].data

        objData = np.array(objData[y1:y2, x1:x2])

        skyData = np.array(skyData[y1:y2, x1:x2])

        objDataMedian = np.nanmedian(objData)

        skyDataMedian = np.nanmedian(skyData)

        skyDataStd = np.nanstd(skyData)

        objDataStd = np.nanstd(objData)

        x1 = np.percentile(objData, 94)

        x2 = np.percentile(objData, 98)

        y1 = np.percentile(skyData, 94)

        y2 = np.percentile(skyData, 98)

        # print x1, x2, y1, y2

        objData[objData < 500] = np.nan

        objData[objData > 10000] = np.nan

        skyData[skyData < 500] = np.nan

        skyData[skyData > 10000] = np.nan

        newobjDataMedian = np.nanmedian(objData)

        newskyDataMedian = np.nanmedian(skyData)

        firstPart = np.sqrt(np.nansum((objData - newobjDataMedian) ** 2))

        secondPart = np.sqrt(np.nansum((skyData - newskyDataMedian) ** 2))

        denom = firstPart * secondPart

        # print denom

        numer = np.nansum((objData - newobjDataMedian) *
                          (skyData - newskyDataMedian))

        rho = numer / denom

        # print rho

        return rho

    def crossCorrOne(self,
                     ext,
                     objFile,
                     skyFile,
                     y1,
                     y2,
                     x1,
                     x2):

        """
        Def:
        Compute the cross-correlation coefficient for a given object
        and skyfile. Define the pixel range over which to compute the
        coefficient.

        Inputs:
                ext - detector extension, must be either 1, 2, 3
                objFile - Input object file to compute correlation
                skyFile - sky image to compare objFile with
                y1, y2, x1, x2 - the range to compute rho over

        """

        objData = fits.open(objFile)

        objData = objData[0].data

        skyData = fits.open(skyFile)

        skyData = skyData[ext].data

        objData = np.array(objData[y1:y2, x1:x2])

        skyData = np.array(skyData[y1:y2, x1:x2])

        objDataMedian = np.nanmedian(objData)

        skyDataMedian = np.nanmedian(skyData)

        skyDataStd = np.nanstd(skyData)

        objDataStd = np.nanstd(objData)

        objData[objData < 500] = np.nan

        objData[objData > 10000] = np.nan

        skyData[skyData < 500] = np.nan

        skyData[skyData > 10000] = np.nan

        newobjDataMedian = np.nanmedian(objData)

        newskyDataMedian = np.nanmedian(skyData)

        firstPart = np.sqrt(np.nansum((objData - newobjDataMedian) ** 2))

        secondPart = np.sqrt(np.nansum((skyData - newskyDataMedian) ** 2))

        denom = firstPart * secondPart

        # print denom

        numer = np.nansum((objData - newobjDataMedian) *
                          (skyData - newskyDataMedian))

        rho = numer / denom

        print rho

        return rho

    def shiftImage(self,
                   ext,
                   infile,
                   skyfile,
                   badpmap,
                   interp_type,
                   stepsize,
                   xmin,
                   xmax,
                   ymin,
                   ymax):

        """
        Def:
        Compute the correlation coefficient for a
        grid of pixel shift values and decide which one is best
        (if better than the original) and apply this to
        the object image to align with the sky. First because we use only the
        first extension because of the way the shiftImageSegments function is
        defined. This function now applies the bad pixel map both before cross
        correlating and before interpolation - need to ignore bad pixels.

        Inputs:
                ext - detector extension, must be either 1, 2, 3
                infile - Input object file to shift
                skyFile - sky image to compare objFile with
                interp_type - type of interpolation function for the shift.

                    -'nearest': nearest neighbour
                    -'linear':bilinear x,y, interpolation
                    -'poly3':third order interior polynomial
                    -'poly5':fifth order interior polynomial
                    -'spline3':third order spline3

                stepsize - value to increment grid by each time, increasing
                this increases the time taken for the computation
                xmin, xmax, ymin, ymax - Grid extremes for brute force shift

        """

        # We first want to apply the bad pixel map

        objTable = fits.open(infile)

        objData = objTable[ext].data

        skyTable = fits.open(skyfile)

        skyData = skyTable[ext].data

        badpTable = fits.open(badpmap)

        badpData = badpTable[ext].data

        # Find the headers of the primary HDU and chosen extension

        objPrimHeader = objTable[0].header

        objExtHeader = objTable[ext].header

        skyPrimHeader = skyTable[0].header

        skyExtHeader = skyTable[ext].header

        badpPrimHeader = badpTable[0].header

        badpExtHeader = badpTable[ext].header

        print (objPrimHeader)
        print (objExtHeader)
        print (skyPrimHeader)
        print (skyExtHeader)
        print (badpPrimHeader)
        print (badpExtHeader)

        # Find the coordinates of the bad pixels and the slitlets

        bad_pixel_coords = np.where(badpData == 0)

        # Loop around the bad pixel locations and mask off
        # on the manObjData and manSkyData

        for i in range(len(bad_pixel_coords[0])):

            # Because of the way np.where works, need
            # to define the x and y coords in this way

            xcoord = bad_pixel_coords[0][i]

            ycoord = bad_pixel_coords[1][i]

            # Now set all positions where there is a
            # dead pixel to np.nan in the object and sky

            objData[xcoord][ycoord] = np.nan

            skyData[xcoord][ycoord] = np.nan

        # Define the minimum and maximum ranges for the correlation

        xMinCorr = (1 * len(objData[0])) / 4

        xMaxCorr = (3 * len(objData[0])) / 4

        yMinCorr = (1 * len(objData)) / 4

        yMaxCorr = (3 * len(objData)) / 4

        # Write out to new temporary fits files - annoyingly need to have
        # the data in fits files to be able to use pyraf functions

        # OBJECT
        objhdu = fits.PrimaryHDU(header=objPrimHeader)

        objhdu.writeto('maskedObj.fits',
                       clobber=True)

        fits.append('maskedObj.fits',
                    data=objData,
                    header=objExtHeader)

        # SKY
        skyhdu = fits.PrimaryHDU(header=skyPrimHeader)

        skyhdu.writeto('maskedSky.fits',
                       clobber=True)

        fits.append('maskedSky.fits',
                    data=skyData,
                    header=skyExtHeader)

        # First compute the correlation coefficient
        # with just the newly saved fits file

        rhoArray = []

        rhoArray.append(self.crossCorrFirst('maskedObj.fits',
                                            'maskedSky.fits',
                                            yMinCorr,
                                            yMaxCorr,
                                            xMinCorr,
                                            xMaxCorr))

        print rhoArray

        # Working. Now create grid of fractional shift values.

        xArray = np.arange(xmin, xmax, stepsize)

        xArray = np.around(xArray, decimals=4)

        yArray = np.arange(ymin, ymax, stepsize)

        yArray = np.around(yArray, decimals=4)

        # Set up mesh grid of rho values for contour plot

        rhoGrid = np.zeros(shape=(len(xArray), len(yArray)))

        # Before attempting the interpolation, we
        # want to mask the bad pixel values,
        # and save to a fresh temporary fits file.
        # Loop over all values in the grid, shift the image by this
        # amount each time and compute the correlation coefficient

        successDict = {}

        for i in range(len(xArray)):

            for j in range(len(yArray)):

                # Perform the shift

                infileName = 'maskedObj.fits[1]'

                pyraf.iraf.imshift(input=infileName,
                                   output='temp_shift.fits',
                                   xshift=xArray[i],
                                   yshift=yArray[j],
                                   interp_type=interp_type)

                # re-open the shifted file and compute rho

                rho = self.crossCorrZeroth('temp_shift.fits',
                                           'maskedSky.fits',
                                           yMinCorr,
                                           yMaxCorr,
                                           xMinCorr,
                                           xMaxCorr)

                rhoGrid[i][j] = rho

                # If the correlation coefficient improves, append to new array

                if rho > rhoArray[0]:

                    print 'SUCCESS, made improvement!'

                    entryName = str(xArray[i]) + ' and ' + str(yArray[j])

                    entryValue = [round(xArray[i], 3), round(yArray[j], 3)]

                    successDict[str(round(rho, 4))] = entryValue

                rhoArray.append(rho)

                # Clean up by deleting the created temporary fits file

                os.system('rm temp_shift.fits')

                # Go back through loop, append next value of rho

                print 'Finished shift: %s %s, rho = %s ' % (xArray[i],
                                                            yArray[j],
                                                            rho)

        os.system('rm maskedObj.fits')

        os.system('rm maskedSky.fits')

        # Now we want to choose the best shift value and actually apply this
        # Need to find the x and y shift values
        # which correspond to the maximum rho
        # Only do this if the success dictionary is not
        # empty, if it is empty return 0.0,0.0

        print rhoGrid

        plt.contour(xArray,
                    yArray,
                    rhoGrid,
                    levels=[(np.max(rhoGrid) - (0.25 * np.std(rhoGrid))),
                            (np.max(rhoGrid) - (1 * np.std(rhoGrid))),
                            (np.max(rhoGrid) - (1.5 * np.std(rhoGrid)))])

        plt.xlabel('$\Delta x$')

        plt.ylabel('$\Delta y$')

        plt.title('Correlation Coefficient Grid')

        plt.savefig('Correlation_coefficient.png')

        plt.close('all')

        print 'Made plot Successfully'

        if successDict:

            print 'Finding Best Shift Value...'

            rhoMax = str(round(max(rhoArray), 4))

            # print rhoMax

            shiftVector = successDict[rhoMax]

            return shiftVector

        else:
            print 'No Shift Value Found'

            shiftVector = [0.0, 0.0]

            return shiftVector

    def shiftImageFirst(self,
                        ext,
                        infile,
                        skyfile,
                        badpmap,
                        interp_type,
                        stepsize,
                        xmin,
                        xmax,
                        ymin,
                        ymax):

        """
        Def:
        Compute the correlation coefficient for a grid
        of pixel shift values and decide which one is best
        (if better than the original) and apply this to
        the object image to align with the sky. First because we use only the
        first extension because of the way the shiftImageSegments function is
        defined. This function now applies the bad pixel map both before cross
        correlating and before interpolation - need to ignore bad pixels.

        Inputs:
                ext - detector extension, must be either 1, 2, 3
                infile - Input object file to shift
                skyFile - sky image to compare objFile with
                interp_type - type of interpolation function for the shift.

                    -'nearest': nearest neighbour
                    -'linear':bilinear x,y, interpolation
                    -'poly3':third order interior polynomial
                    -'poly5':fifth order interior polynomial
                    - 'spline3':third order spline3

                stepsize - value to increment grid by each time, increasing
                this increases the time taken for the computation
                xmin, xmax, ymin, ymax - Grid extremes for brute force shift

        Output:
                ShiftVector - list of shifts to apply to the infile
        """
        # We first want to apply the bad pixel map

        objTable = fits.open(infile)

        objData = objTable[1].data

        skyTable = fits.open(skyfile)

        skyData = skyTable[1].data

        badpTable = fits.open(badpmap)

        badpData = badpTable[1].data

        # Find the headers of the primary HDU and chosen extension

        objPrimHeader = objTable[0].header

        objExtHeader = objTable[1].header

        skyPrimHeader = skyTable[0].header

        skyExtHeader = skyTable[1].header

        badpPrimHeader = badpTable[0].header

        badpExtHeader = badpTable[1].header

        temp = sys.stdout

        sys.stdout = open('log.txt', 'w')

        print (objPrimHeader)
        print (objExtHeader)
        print (skyPrimHeader)
        print (skyExtHeader)
        print (badpPrimHeader)
        print (badpExtHeader)

        sys.stdout.close()

        sys.stdout = temp

        badpData[badpData == 0] = np.nan

        objData = objData * badpData

        skyData = skyData * badpData

        # Define the minimum and maximum ranges for the correlation

        xMinCorr = (1 * len(objData[0])) / 4

        xMaxCorr = (3 * len(objData[0])) / 4

        yMinCorr = (1 * len(objData)) / 4

        yMaxCorr = (3 * len(objData)) / 4

        # Write out to new temporary fits files - annoyingly need to have
        # the data in fits files to be able to use pyraf functions

        # OBJECT

        objhdu = fits.PrimaryHDU(header=objPrimHeader)

        objhdu.writeto('maskedObj.fits',
                       clobber=True)

        fits.append('maskedObj.fits',
                    data=objData,
                    header=objExtHeader)

        # SKY

        skyhdu = fits.PrimaryHDU(header=skyPrimHeader)

        skyhdu.writeto('maskedSky.fits',
                       clobber=True)

        fits.append('maskedSky.fits',
                    data=skyData,
                    header=skyExtHeader)

        # First compute the correlation coefficient
        # with just the newly saved fits file

        rhoArray = []

        rhoArray.append(self.crossCorrFirst('maskedObj.fits',
                                            'maskedSky.fits',
                                            yMinCorr,
                                            yMaxCorr,
                                            xMinCorr,
                                            xMaxCorr))

        print rhoArray

        # Working. Now create grid of fractional shift values.

        xArray = np.arange(xmin, xmax, stepsize)

        xArray = np.around(xArray, decimals=4)

        yArray = np.arange(ymin, ymax, stepsize)

        yArray = np.around(yArray, decimals=4)

        thArray = np.arange(-0.05, 0.05, 0.01)

        thArray = np.around(thArray, decimals=4)

        x, y = np.meshgrid(xArray, yArray)

        # Set up mesh grid of rho values for contour plot

        rhoGrid = np.zeros(shape=(len(x), len(x[0])))

        # Before attempting the interpolation,
        # we want to mask the bad pixel values,
        # and save to a fresh temporary fits file.

        # Loop over all values in the grid, shift the image by this
        # amount each time and compute the correlation coefficient

        successDict = {}

        for i in range(len(xArray)):

            for j in range(len(yArray)):

                # Perform the shift

                infileName = 'maskedObj.fits[1]'

                pyraf.iraf.imshift(input=infileName,
                                   output='temp_shift.fits',
                                   xshift=xArray[i],
                                   yshift=yArray[j],
                                   interp_type=interp_type)

                # re-open the shifted file and compute rho

                rho = self.crossCorrZeroth('temp_shift.fits',
                                           'maskedSky.fits',
                                           yMinCorr,
                                           yMaxCorr,
                                           xMinCorr,
                                           xMaxCorr)

                rhoGrid[j][i] = rho

                # If the correlation coefficient improves, append to new array

                if rho > rhoArray[0]:

                    print 'SUCCESS, made improvement!'

                    entryName = str(xArray[i]) + ' and ' + str(yArray[j])

                    entryValue = [round(xArray[i], 3), round(yArray[j], 3)]

                    successDict[str(round(rho, 4))] = entryValue

                rhoArray.append(rho)

                # Clean up by deleting the created temporary fits file

                os.system('rm temp_shift.fits')

                # Go back through loop, append next value of rho

                print 'Finished shift: %s %s, rho = %s ' % (xArray[i],
                                                            yArray[j],
                                                            rho)

        os.system('rm maskedObj.fits')

        os.system('rm maskedSky.fits')

        # Now we want to choose the best shift value and actually apply this
        # Need to find the x and y shift values which
        # correspond to the maximum rho
        # Only do this if the success dictionary is not
        # empty, if it is empty return 0.0,0.0

        plt.contour(x,
                    y,
                    rhoGrid,
                    levels=[(np.max(rhoGrid) - (0.25 * np.std(rhoGrid))),
                            (np.max(rhoGrid) - (1 * np.std(rhoGrid))),
                            (np.max(rhoGrid) - (1.5 * np.std(rhoGrid)))])

        plt.xlabel('$\Delta x$')

        plt.ylabel('$\Delta y$')

        plt.title('Correlation Coefficient Grid')

        plotName = infile[:-5] + '_' + str(ext) + '_' + \
            interp_type + '_CorrelationGraph.png'

        plt.savefig(plotName)
        plt.close('all')

        print 'Made plot Successfully'

        if successDict:

            print 'Finding Best Shift Value...'

            rhoMax = str(round(max(rhoArray), 4))

            # print rhoMax

            shiftVector = successDict[rhoMax]

            return shiftVector

        else:

            print 'No Shift Value Found'

            shiftVector = [0.0, 0.0]

            return shiftVector

    def imSplit(self,
                ext,
                infile,
                vertSegments,
                horSegments):

        """
        Def:
        Take an input image file and split up into a series of squares
        Main purpose is for use in the shiftImageSegments functino

        Input:
                infile - file to be divided
                vertSegments - number of vertical segments
                              (2048 must be divisible by)
                horSegments - number of horizontal segments
                             (2048 must be divisible by)
                ext - extension number

        Output:
                segmentArray - 1D array containing 2D square array segments

        """

        # Read in the data file at the given extension

        data = fits.open(infile)

        data = data[ext].data

        # Initialise the empty array

        segmentArray = []

        # We can't do this if 2048 isn't divisible by the segments
        # write an error function to check that this is the case
        # And exit if it isn't

        if ((2048 % vertSegments != 0) or (2048 % horSegments != 0)):

            raise ValueError('Please ensure that 2048 is'
                             + ' divisible by your segment choice')

        # Counters for the horizontal slicing

        hor1 = 0

        hor2 = (2048 / horSegments)

        for j in range(horSegments):

            # Counters for the vertical slicing

            x = 0

            y = (2048 / vertSegments)

            for i in range(vertSegments):

                # Slice the data according to user selection

                segmentArray.append(data[hor1:hor2, x:y])

                x += (2048 / vertSegments)

                y += (2048 / vertSegments)

            hor1 += (2048 / horSegments)

            hor2 += (2048 / horSegments)

        # print segmentArray

        return segmentArray

    def shiftImageSegments(self,
                           ext,
                           infile,
                           skyfile,
                           badpmap,
                           vertSegments,
                           horSegments,
                           interp_type,
                           stepsize,
                           xmin,
                           xmax,
                           ymin,
                           ymax):

        """
        Def:
        Lots of arguments because of using lots of different functions.
        This is taking an object and a sky image,
        splitting them into a specified
        number of segments, performing shifts to each of the segments and then
        computing the cross-correlation function to see if we can improve the
        alignment at all. Should give better results than a global shift

        Inputs:
                infile - file to be divided
                skyFile - sky image to compare objFile with
                vertSegments - number of vertical segments
                horSegments - number of horizontal segments
                ext - extension number
                interp_type - type of interpolation function for the shift.

                    -'nearest': nearest neighbour
                    -'linear':bilinear x,y, interpolation
                    -'poly3':third order interior polynomial
                    -'poly5':fifth order interior polynomial
                    -'spline3':third order spline3

                stepsize - value to increment grid by each time, increasing
                        this increases the time taken for the computation
                xmin, xmax, ymin, ymax - Grid extremes for brute force shift

        Output:
                reconstructedData - array containing the shifted data

        """

        # Create arrays of the split files using the imSplit function

        objArray = self.imSplit(ext, infile, vertSegments, horSegments)

        skyArray = self.imSplit(ext, skyfile, vertSegments, horSegments)

        badpArray = self.imSplit(ext, badpmap, vertSegments, horSegments)

        # Find the headers of the primary HDU and chosen extension

        objTable = fits.open(infile)

        objPrimHeader = objTable[0].header

        objExtHeader = objTable[ext].header

        skyTable = fits.open(skyfile)

        skyPrimHeader = skyTable[0].header

        skyExtHeader = skyTable[ext].header

        badpTable = fits.open(badpmap)

        badpPrimHeader = badpTable[0].header

        badpExtHeader = badpTable[ext].header

        temp = sys.stdout

        sys.stdout = open('log.txt', 'w')

        print (objPrimHeader)
        print (objExtHeader)
        print (skyPrimHeader)
        print (skyExtHeader)
        print (badpPrimHeader)
        print (badpExtHeader)

        sys.stdout.close()

        sys.stdout = temp

        shiftArray = []

        # Should now have two 1D arrays of 2D arrays of equal size

        for i in range(len(objArray)):

            tempObjName = infile[:-5] + str(vertSegments) + \
                str(horSegments) + str(i) + '_temp.fits'

            # Write out to new temporary fits files - annoyingly need to have
            # the data in fits files to be able to use pyraf functions

            # OBJECT

            objhdu = fits.PrimaryHDU(header=objPrimHeader)

            objhdu.writeto(tempObjName,
                           clobber=True)

            fits.append(tempObjName,
                        data=objArray[i],
                        header=objExtHeader)

            # SKY

            skyhdu = fits.PrimaryHDU(header=skyPrimHeader)

            skyhdu.writeto('tempSky.fits',
                           clobber=True)

            fits.append('tempSky.fits',
                        data=skyArray[i],
                        header=skyExtHeader)

            # BADPIXEL

            badphdu = fits.PrimaryHDU(header=badpPrimHeader)

            badphdu.writeto('tempbadp.fits',
                            clobber=True)

            fits.append('tempbadp.fits',
                        data=badpArray[i],
                        header=badpExtHeader)

            # Now need to apply the shiftImageFirst
            # function, which compares the chosen
            # extension shifted object and sky files.
            # Create an array to hold the shift coordinates
            # for each segment. This is defined
            # Outside the for loop so that I am not initialising it every time.

            print 'This is shift: %s' % i

            shiftArray.append(self.shiftImageFirst(ext,
                                                   tempObjName,
                                                   'tempSky.fits',
                                                   'tempbadp.fits',
                                                   interp_type,
                                                   stepsize,
                                                   xmin,
                                                   xmax,
                                                   ymin,
                                                   ymax))

            # Clean up the temporary fits files during each part of the loop

            os.system('rm %s' % tempObjName)

            os.system('rm tempSky.fits')

            os.system('rm tempbadp.fits')

            # Now just need to vstack all of these arrays and
            # will have a 2048x2048 corrected array

        print shiftArray

        # Now the clever part - to actually apply
        # the shifts to the unmasked infile
        # imshift can be used with a list of infile names,
        # outfile names and shift coordinates
        # If I create these lists I can imshift all at once,
        # read in the data and then recombine
        # First get the x and y vectors for the shift coordinates

        xArray = []

        yArray = []

        for item in shiftArray:

            xArray.append(item[0])

            yArray.append(item[1])

        xArray = np.array(np.around(xArray, 3))

        yArray = np.array(np.around(yArray, 3))

        # Create the .txt file with the two columns
        # specifying the shift coordinates
        np.savetxt('coords.txt',
                   np.c_[xArray, yArray],
                   fmt=('%5.3f', '%5.3f'))

        # Coordinates list sorted. Now need list
        # of input files and input file names
        # To do this need to go back to the
        # imSplit method with the unmasked file
        # Write to a list of temporary fits files, which will become temporary
        # Output fits files, which will be read
        # back in as data before recombining
        # Must clean everything up at the end by removing with os.system()
        # Create arrays of the split files using the imSplit function

        # Need to apply the shifts to a masked
        # object file. Open the bad pixel map
        # and the object file and mask the pixels and save to temporary file
        # Find the coordinates of the bad pixels and the slitlets

        objData = objTable[ext].data

        badpData = badpTable[ext].data

        badpData[badpData == 0] = np.nan

        objData = objData * badpData

        # Write out to new file which will then be read in to split up the data
        objhdu = fits.PrimaryHDU(header=objPrimHeader)

        objhdu.writeto('temp_masked.fits',
                       clobber=True)

        fits.append('temp_masked.fits',
                    data=objData,
                    header=objExtHeader)

        objArray = self.imSplit(1,
                                'temp_masked.fits',
                                vertSegments,
                                horSegments)

        inFileArray = []

        outFileArray = []

        shiftedDataArray = []

        vstackArray = []

        hstackArray = []

        # Should now have two 1D arrays of 2D arrays of equal size

        for i in range(len(objArray)):

            inFileName = 'tempObjin'+str(i)+'.fits'

            outFileName = 'tempObjout'+str(i)+'.fits'

            inFileArray.append(inFileName)

            outFileArray.append(outFileName)

            # Write out to new temporary fits files - annoyingly need to have
            # the data in fits files to be able to use pyraf functions

            # OBJECT

            objhdu = fits.PrimaryHDU(header=objPrimHeader)

            objhdu.writeto(inFileName,
                           clobber=True)

            fits.append(inFileName,
                        data=objArray[i],
                        header=objExtHeader)

            inFileName = inFileName + '[1]'

            # Now apply imshift with all the parameters
            pyraf.iraf.imshift(input=inFileName,
                               output=outFileName,
                               xshift=xArray[i],
                               yshift=yArray[i],
                               interp_type=interp_type)

            # We want a 1D array of 2D arrays again,
            # read the data files back in

            data = fits.open(outFileName)

            data = data[0].data

            shiftedDataArray.append(data)

            # Go back to the top of the loop and grab the next file

        # The final problem is that we have a 1D arrays of 2D arrays that needs
        # to be recombined into the original 2048 x 2048 which created it
        # Ordering depends on the number of vertical and horizontal segments

        x = len(shiftedDataArray) / horSegments

        a = 0

        while x <= len(shiftedDataArray):

            for i in range(a, x):

                print i

                hstackArray.append(shiftedDataArray[i])

            vstackArray.append(np.hstack(hstackArray))

            hstackArray = []

            x += len(shiftedDataArray) / horSegments

            a += len(shiftedDataArray) / horSegments

        # Reconstruct by vstacking the final array

        reconstructedData = np.vstack(vstackArray)

        # Clean up by getting rid of uneeded files

        for item in inFileArray:

            os.system('rm %s' % item)

        for item in outFileArray:

            os.system('rm %s' % item)

        os.system('rm coords.txt')

        os.system('rm temp_masked.fits')

        os.system('rm log.txt')

        return reconstructedData

    def shiftAllExtensions(self,
                           infile,
                           skyfile,
                           badpmap,
                           vertSegments,
                           horSegments,
                           interp_type,
                           stepsize,
                           xmin,
                           xmax,
                           ymin,
                           ymax):

        """
        Def:
        Uses the shiftImageSegments method for
        each extensions and then combines
        all of these together into a single shifted fits file

        Inputs:
                infile - file to be divided
                skyFile - sky image to compare objFile with
                vertSegments - number of vertical segments
                horSegments - number of horizontal segments
                ext - extension number
                interp_type - type of interpolation function for the shift.

                    -'nearest': nearest neighbour
                    -'linear':bilinear x,y, interpolation
                    -'poly3':third order interior polynomial
                    -'poly5':fifth order interior polynomial
                    -'spline3':third order spline3

                stepsize - value to increment grid by each time, increasing
                        this increases the time taken for the computation
                xmin, xmax, ymin, ymax - Grid extremes for brute force shift

        Output:
                infile_Shifted.fits - containing all three shifted extensions

        """

        # Prepare the headers for writing out the fits file

        objTable = fits.open(infile)

        objPrimHeader = objTable[0].header

        objExtHeader1 = objTable[1].header

        objExtHeader2 = objTable[2].header

        objExtHeader3 = objTable[3].header

        temp = sys.stdout

        sys.stdout = open('log.txt', 'w')

        print (objPrimHeader)
        print (objExtHeader1)
        print (objExtHeader2)
        print (objExtHeader3)

        sys.stdout.close()

        sys.stdout = temp

        # Set up the array

        reconstructedDataArray = []

        # Use the shifted image segment function

        for i in range(1, 4):

            print 'Shifting Extension: %s' % i

            reconstructedDataArray.append(
                self.shiftImageSegments(i,
                                        infile,
                                        skyfile,
                                        badpmap,
                                        vertSegments,
                                        horSegments,
                                        interp_type,
                                        stepsize,
                                        xmin,
                                        xmax,
                                        ymin,
                                        ymax))

        # Name the shifted data file

        shiftedName = infile[:-5] + '_' + str(vertSegments) + \
            str(horSegments) + '_' + interp_type + '_Shifted.fits'

        print 'Saving %s' % shiftedName

        objhdu = fits.PrimaryHDU(header=objPrimHeader)

        objhdu.writeto(shiftedName,
                       clobber=True)

        fits.append(shiftedName,
                    data=reconstructedDataArray[0],
                    header=objExtHeader1)

        fits.append(shiftedName,
                    data=reconstructedDataArray[1],
                    header=objExtHeader2)

        fits.append(shiftedName,
                    data=reconstructedDataArray[2],
                    header=objExtHeader3)

    def applyShiftAllExtensions(self,
                                fileList,
                                badpmap,
                                vertSegments,
                                horSegments,
                                interp_type,
                                stepsize,
                                xmin,
                                xmax,
                                ymin,
                                ymax):

        # Read in the data from the fileList

        data = np.genfromtxt(fileList, dtype='str')

        # Save the names and types as lists

        names = data[0:, 0]

        types = data[0:, 1]

        # Loop round all names and apply the computeOffsetSegments method

        for i in range(1, len(names)):

            if types[i] == 'O':

                objFile = names[i]

                if i == 1:

                    skyFile = names[i + 1]

                elif types[i - 1] == 'S':

                    skyFile = names[i - 1]

                else:

                    skyFile = names[i + 1]

                print 'Shifting file: %s : %s' % (i, objFile)

                # Now use the method defined within this class

                self.shiftAllExtensions(objFile,
                                        skyFile,
                                        badpmap,
                                        vertSegments,
                                        horSegments,
                                        interp_type,
                                        stepsize,
                                        xmin,
                                        xmax,
                                        ymin,
                                        ymax)

                # Which will loop through all and save
                # the corrected object file
                # as objectFile_Corrected.fits.
                # These are then fed through the pipeline.

    def shiftPlot(self,
                  coords_file):

        """
        Def:
        Plotting function, takes the output from applyShiftAllExtensionsMin and
        plots the x and y shift vectors against the detector ID

        Input: 
                coords_file - directly from shiftAllExtensionsMin
        Output: 
                shift_plot.png
        """

        # Load the coordinates

        coords = np.loadtxt(coords_file)

        # The coordinates are 2D arrays. Need to set up the vectors

        d_x = []

        d_y = []

        # x_vector

        for entry in coords:

            for i in (0, 2, 4):

                d_x.append(entry[i])

                i + 1

        # y_vector

        for entry in coords:

            for i in (1, 3, 5):

                d_y.append(entry[i])

                i + 1

        # Set the frame ID vectors

        f_ID = np.arange(0,
                         len(d_x),
                         1)

        # Now make the plots for both nights,
        # want the same x-axis for all three layers

        f, (ax2, ax3) = plt.subplots(2, sharex=True, figsize=(18.0, 10.0))

        ax2.plot(f_ID, d_x, color='g')
        ax2.set_title('x shift (pixels) vs. ID', fontsize=24)
        ax2.grid(b=True, which='both', linestyle='--')
        ax2.tick_params(axis='both', which='major', labelsize=15)

        ax3.plot(f_ID, d_y, color='r')
        ax3.set_title('y shift (pixels) vs. ID', fontsize=24)
        ax3.set_xlabel('Detector ID', fontsize=20)
        ax3.set_xticks((np.arange(min(f_ID), max(f_ID)+1, 1.0)))
        ax3.tick_params(axis='both', which='major', labelsize=15)
        ax3.grid(b=True, which='both', linestyle='--')

        plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
        plt.savefig('shift_plot.png')
        plt.close('all')

    #  alternate set of functions using minimisation to find the
    #  best shift position

    def shiftImageFirstMin(self,
                           xArray,
                           infile,
                           skyfile,
                           badpmap,
                           interp_type):

        """
        Def:
        Compute the correlation coefficient
        for a grid of pixel shift values and
        decide which one is best (if better than the original)
        and apply this to the object image to align with the sky.
        First because we use only the
        first extension because of the way the
        shiftImageSegments function is
        defined. This function now applies the
        bad pixel map both before cross
        correlating and before interpolation - need to ignore bad pixels.

        Inputs:
                ext - detector extension, must be either 1, 2, 3
                infile - Input object file to shift
                skyFile - sky image to compare objFile with
                interp_type - type of interpolation function for the shift.

                    -'nearest': nearest neighbour
                    -'linear':bilinear x,y, interpolation
                    -'poly3':third order interior polynomial
                    -'poly5':fifth order interior polynomial
                    -'spline3':third order spline3

                stepsize - value to increment grid by each time, increasing
                        this increases the time taken for the computation
                xmin, xmax, ymin, ymax - Grid extremes for brute force shift

        """

        # We first want to apply the bad pixel map

        objTable = fits.open(infile)

        objData = objTable[1].data

        skyTable = fits.open(skyfile)

        skyData = skyTable[1].data

        badpTable = fits.open(badpmap)

        badpData = badpTable[1].data

        # Find the headers of the primary HDU and chosen extension

        objPrimHeader = objTable[0].header

        objExtHeader = objTable[1].header

        skyPrimHeader = skyTable[0].header

        skyExtHeader = skyTable[1].header

        badpPrimHeader = badpTable[0].header

        badpExtHeader = badpTable[1].header

        temp = sys.stdout

        sys.stdout = open('log.txt', 'w')

        print (objPrimHeader)
        print (objExtHeader)
        print (skyPrimHeader)
        print (skyExtHeader)
        print (badpPrimHeader)
        print (badpExtHeader)

        sys.stdout.close()

        sys.stdout = temp

        # Find the coordinates of the bad pixels and the slitlets
        # Instead of looping, do much faster multiplication
        # Want to avoid loops at all costs

        badpData[badpData == 0] = np.nan

        objData = objData * badpData

        skyData = skyData * badpData

        # Define the minimum and maximum ranges for the correlation

        xMinCorr = (1 * len(objData[0])) * 0.0625
        xMaxCorr = (15 * len(objData[0])) * 0.0625
        yMinCorr = (1 * len(objData)) * 0.0625
        yMaxCorr = (15 * len(objData)) * 0.0625

        # Write out to new temporary fits files - annoyingly need to have
        # the data in fits files to be able to use pyraf functions

        # OBJECT

        objhdu = fits.PrimaryHDU(header=objPrimHeader)

        objhdu.writeto('maskedObj.fits',
                       clobber=True)

        fits.append('maskedObj.fits',
                    data=objData,
                    header=objExtHeader)

        # SKY

        skyhdu = fits.PrimaryHDU(header=skyPrimHeader)

        skyhdu.writeto('maskedSky.fits',
                       clobber=True)

        fits.append('maskedSky.fits',
                    data=skyData,
                    header=skyExtHeader)

        # Before attempting the interpolation,
        # we want to mask the bad pixel values,
        # and save to a fresh temporary fits file.
        # Loop over all values in the grid, shift the image by this

        infileName = 'maskedObj.fits[1]'

        print '[INFO]: Shifting: %s %s' % (xArray[0], xArray[1])

        pyraf.iraf.imshift(input=infileName,
                           output='temp_shift.fits',
                           xshift=xArray[0],
                           yshift=xArray[1],
                           interp_type=interp_type)

        # re-open the shifted file and compute rho

        rho = self.crossCorrZeroth('temp_shift.fits',
                                   'maskedSky.fits',
                                   yMinCorr,
                                   yMaxCorr,
                                   xMinCorr,
                                   xMaxCorr)

        # Tidy up

        os.system('rm temp_shift.fits')
        os.system('rm maskedObj.fits')
        os.system('rm maskedSky.fits')

        # Return the correlation coefficient

        print (1.0 / rho)

        return (1.0 / rho)

    def minimiseRho(self,
                    infile,
                    skyfile,
                    badpmap,
                    interp_type):

        # Minimise the function recipShift with respect to x and y
        # Define the shift starting points

        x0 = [1.0, 1.0]

        minimizer_kwargs = {'method': 'Nelder-Mead', 'args': (infile,
                                                              skyfile,
                                                              badpmap,
                                                              interp_type,)}

        # First Method - using simple downhill simplex
        res = minimize(self.shiftImageFirstMin,
                       x0,
                       args=(infile,
                             skyfile,
                             badpmap,
                             interp_type,),
                       method = 'Nelder-Mead',
                       tol=0.005,
                       options={'disp': True})

        # Return the shift array which minimises
        # the inverse correlation coefficient

        print res.x

        return res.x


    def shiftImageSegmentsMin(self,
                              ext,
                              infile,
                              skyfile,
                              badpmap,
                              vertSegments,
                              horSegments,
                              interp_type):

        """
        Def:
        Lots of arguments because of using lots of different functions.
        This is taking an object and a sky image,
        splitting them into a specified
        number of segments, performing shifts to each of the segments and then
        computing the cross-correlation function to see if we can improve the
        alignment at all. Should give better results than a global shift

        Inputs:
                infile - file to be divided
                skyFile - sky image to compare objFile with
                vertSegments - number of vertical segments
                horSegments - number of horizontal segments
                ext - extension number
                interp_type - type of interpolation function for the shift.

                    -'nearest': nearest neighbour
                    -'linear':bilinear x,y, interpolation
                    -'poly3':third order interior polynomial
                    -'poly5':fifth order interior polynomial
                    -'spline3':third order spline3

                stepsize - value to increment grid by each time, increasing
                this increases the time taken for the computation
                xmin, xmax, ymin, ymax - Grid extremes for brute force shift

        """

        # Create arrays of the split files using the imSplit function

        objArray = self.imSplit(ext, infile, vertSegments, horSegments)

        skyArray = self.imSplit(ext, skyfile, vertSegments, horSegments)

        badpArray = self.imSplit(ext, badpmap, vertSegments, horSegments)

        # Find the headers of the primary HDU and chosen extension

        objTable = fits.open(infile)

        objPrimHeader = objTable[0].header

        objExtHeader = objTable[ext].header

        skyTable = fits.open(skyfile)

        skyPrimHeader = skyTable[0].header

        skyExtHeader = skyTable[ext].header

        badpTable = fits.open(badpmap)

        badpPrimHeader = badpTable[0].header

        badpExtHeader = badpTable[ext].header

        temp = sys.stdout

        sys.stdout = open('log.txt', 'w')

        print (objPrimHeader)
        print (objExtHeader)
        print (skyPrimHeader)
        print (skyExtHeader)
        print (badpPrimHeader)
        print (badpExtHeader)

        sys.stdout.close()

        sys.stdout = temp

        shiftArray = []

        # Should now have two 1D arrays of 2D arrays of equal size

        for i in range(len(objArray)):

            tempObjName = infile[:-5] + str(vertSegments) + \
                str(horSegments) + str(i) + '_temp.fits'

            # Write out to new temporary fits files - annoyingly need to have
            # the data in fits files to be able to use pyraf functions

            # OBJECT

            objhdu = fits.PrimaryHDU(header=objPrimHeader)

            objhdu.writeto(tempObjName,
                           clobber=True)

            fits.append(tempObjName,
                        data=objArray[i],
                        header=objExtHeader)

            # SKY

            skyhdu = fits.PrimaryHDU(header=skyPrimHeader)

            skyhdu.writeto('tempSky.fits',
                           clobber=True)

            fits.append('tempSky.fits',
                        data=skyArray[i],
                        header=skyExtHeader)

            # BADPIXEL

            badphdu = fits.PrimaryHDU(header=badpPrimHeader)

            badphdu.writeto('tempbadp.fits',
                            clobber=True)

            fits.append('tempbadp.fits',
                        data=badpArray[i],
                        header=badpExtHeader)

            # Now need to apply the shiftImageFirst
            # function, which compares the chosen
            # extension shifted object and sky files.
            # Create an array to hold the shift
            # coordinates for each segment. This is defined
            # Outside the for loop so that I am not initialising it every time.

            xMinCorr = (1 * len(objArray[0])) * 0.0625

            xMaxCorr = (15 * len(objArray[0])) * 0.0625

            yMinCorr = (1 * len(objArray[0])) * 0.0625

            yMaxCorr = (15 * len(objArray[0])) * 0.0625

            print 'Before shifting, the correlation is: %s ' % \
                (1.0 / self.crossCorrFirst(tempObjName,
                                           'tempSky.fits',
                                           yMinCorr,
                                           yMaxCorr,
                                           xMinCorr,
                                           xMaxCorr))

            print '[INFO]: This is shift: %s' % i

            shiftArray.append(self.minimiseRho(tempObjName,
                                               'tempSky.fits',
                                               'tempbadp.fits',
                                               interp_type))

            # Clean up the temporary fits files during each part of the loop

            os.system('rm %s' % tempObjName)
            os.system('rm tempSky.fits')
            os.system('rm tempbadp.fits')

            # Now just need to vstack all of these arrays
            # and will have a 2048x2048 corrected array

        print shiftArray

        # apply the shifts to the unmasked infile
        # imshift can be used with a list of infile names,
        # outfile names and shift coordinates
        # create these lists and imshift all at once,
        # read in the data and then recombine
        # First get the x and y vectors for the shift coordinates

        xArray = []

        yArray = []

        for item in shiftArray:

            xArray.append(item[0])

            yArray.append(item[1])

        xArray = np.array(np.around(xArray, 3))

        yArray = np.array(np.around(yArray, 3))

        # Create the .txt file with the two columns
        # specifying the shift coordinates

        np.savetxt('coords.txt', np.c_[xArray, yArray], fmt=('%5.3f', '%5.3f'))

        # Need to apply the shifts to a masked
        # object file. Open the bad pixel map
        # and the object file and mask the pixels and save to temporary file
        # Find the coordinates of the bad pixels and the slitlets

        objData = objTable[ext].data

        badpData = badpTable[ext].data

        badpData[badpData == 0] = np.nan

        objData = objData * badpData

        # Write out to new file which will then
        # be read in to split up the data

        objhdu = fits.PrimaryHDU(header=objPrimHeader)

        objhdu.writeto('temp_masked.fits',
                       clobber=True)

        fits.append('temp_masked.fits',
                    data=objData,
                    header=objExtHeader)

        objArray = self.imSplit(1,
                                'temp_masked.fits',
                                vertSegments,
                                horSegments)

        inFileArray = []

        outFileArray = []

        shiftedDataArray = []

        vstackArray = []

        hstackArray = []

        # Should now have two 1D arrays of 2D arrays of equal size

        for i in range(len(objArray)):

            inFileName = 'tempObjin' + str(i) + '.fits'

            outFileName = 'tempObjout' + str(i) + '.fits'

            inFileArray.append(inFileName)

            outFileArray.append(outFileName)

            # Write out to new temporary fits files - annoyingly need to have
            # the data in fits files to be able to use pyraf functions

            # OBJECT

            objhdu = fits.PrimaryHDU(header=objPrimHeader)

            objhdu.writeto(inFileName,
                           clobber=True)

            fits.append(inFileName,
                        data=objArray[i],
                        header=objExtHeader)

            inFileName = inFileName + '[1]'

            # Now apply imshift with all the parameters
            pyraf.iraf.imshift(input=inFileName,
                               output=outFileName,
                               xshift=xArray[i],
                               yshift=yArray[i],
                               interp_type=interp_type)

            # We want a 1D array of 2D arrays again,
            # read the data files back in

            data = fits.open(outFileName)

            data = data[0].data

            shiftedDataArray.append(data)

            # Go back to the top of the loop and grab the next file

        # The final problem is that we have a 1D arrays of 2D arrays that needs
        # to be recombined into the original 2048 x 2048 which created it
        # Ordering depends on the number of vertical and horizontal segments

        x = len(shiftedDataArray) / horSegments

        a = 0

        while x <= len(shiftedDataArray):

            for i in range(a, x):

                print i

                hstackArray.append(shiftedDataArray[i])

            vstackArray.append(np.hstack(hstackArray))

            hstackArray = []

            x += len(shiftedDataArray) / horSegments

            a += len(shiftedDataArray) / horSegments

        # Reconstruct by vstacking the final array

        reconstructedData = np.vstack(vstackArray)

        # Clean up by getting rid of uneeded files

        for item in inFileArray:

            os.system('rm %s' % item)

        for item in outFileArray:

            os.system('rm %s' % item)

        os.system('rm coords.txt')
        os.system('rm temp_masked.fits')
        os.system('rm log.txt')

        return reconstructedData, shiftArray

    def shiftAllExtensionsMin(self,
                              infile,
                              skyfile,
                              badpmap,
                              vertSegments,
                              horSegments,
                              interp_type):

        """
        Def:
        Uses the shiftImageSegments method for
        each extensions and then combines
        all of these together into a single shifted fits file

        Inputs:
                infile - file to be divided
                skyFile - sky image to compare objFile with
                vertSegments - number of vertical segments
                horSegments - number of horizontal segments
                ext - extension number
                interp_type - type of interpolation function for the shift.

                    -'nearest': nearest neighbour
                    -'linear':bilinear x,y, interpolation
                    -'poly3':third order interior polynomial
                    -'poly5':fifth order interior polynomial
                    -'spline3':third order spline3

                stepsize - value to increment grid by each time, increasing
                        this increases the time taken for the computation
                xmin, xmax, ymin, ymax - Grid extremes for brute force shift

        Output:
                infile_Shifted.fits - containing all three shifted extensions

        """
        # First do the extra masking of the whole infile,
        # this saves as masked_infile.fits

        self.maskExtraPixels(infile, badpmap)

        # Prepare the headers for writing out the fits file

        objTable = fits.open('masked_infile.fits')

        objPrimHeader = objTable[0].header

        objExtHeader1 = objTable[1].header

        objExtHeader2 = objTable[2].header

        objExtHeader3 = objTable[3].header

        temp = sys.stdout

        sys.stdout = open('log.txt', 'w')

        print (objPrimHeader)
        print (objExtHeader1)
        print (objExtHeader2)
        print (objExtHeader3)

        sys.stdout.close()

        sys.stdout = temp

        # Set up the array

        reconstructedDataArray = []

        ShiftArrayList = []

        # Use the shifted image segment function

        for i in range(1, 4):

            print '[INFO]: Shifting Extension: %s' % i

            reconstructedData, shiftArray = \
                self.shiftImageSegmentsMin(i,
                                           'masked_infile.fits',
                                           skyfile,
                                           badpmap,
                                           vertSegments,
                                           horSegments,
                                           interp_type)

            reconstructedDataArray.append(reconstructedData)

            ShiftArrayList.append(shiftArray)

        # Name the shifted data file
        shiftedName = infile[:-5] + '_' + str(vertSegments) + \
            str(horSegments) + '_Shifted.fits'

        print 'Saving %s' % shiftedName

        objhdu = fits.PrimaryHDU(header=objPrimHeader)

        objhdu.writeto(shiftedName,
                       clobber=True)

        fits.append(shiftedName,
                    data=reconstructedDataArray[0],
                    header=objExtHeader1)

        fits.append(shiftedName,
                    data=reconstructedDataArray[1],
                    header=objExtHeader2)

        fits.append(shiftedName,
                    data=reconstructedDataArray[2],
                    header=objExtHeader3)

        os.system('rm masked_infile.fits')

        return ShiftArrayList

    def applyShiftAllExtensionsMin(self,
                                   fileList,
                                   badpmap,
                                   vertSegments,
                                   horSegments,
                                   interp_type):
        """
        Def:
        Apply the shift all extensions min method to a list of files

        Input:
                fileList - txt file containing the object name in the first
                            column and ocs.arm1.type in second column
                badpmap - badpixel map produced by pipeline
                vertSegments - how many vertical segments for the shift
                horSegments - how many horizontal segments for the shift
                interp_type - type of interpolation function for the shift.

                    -'nearest': nearest neighbour
                    -'linear':bilinear x,y, interpolation
                    -'poly3':third order interior polynomial
                    -'poly5':fifth order interior polynomial
                    -'spline3':third order spline3
        Output:
                infiles_Shifted.fits - with shifts applied to all extensions
        """

        # Check for existence of masked_infile

        if os.path.isfile('masked_infile.fits'):

            os.system('rm masked_infile.fits')

        # Read in the data from the fileList

        data = np.genfromtxt(fileList, dtype='str')

        # Save the names and types as lists

        names = data[0:, 0]

        types = data[0:, 1]

        shiftList = []

        # Loop round all names and apply the computeOffsetSegments method

        for i in range(1, len(names)):

            if types[i] == 'O':

                objFile = names[i]

                if i == 1:

                    skyFile = names[i + 1]

                elif types[i - 1] == 'S':

                    skyFile = names[i - 1]

                else:

                    skyFile = names[i + 1]

                print '[INFO]: Shifting file: %s : %s' % (i, objFile)

                # Now use the method defined within this class
                # Do additional masking before starting the shift

                shiftList.append(self.shiftAllExtensionsMin(objFile,
                                                            skyFile,
                                                            badpmap,
                                                            vertSegments,
                                                            horSegments,
                                                            interp_type))

                # Which will loop through all and save
                # the corrected object file
                # as objectFile_Corrected.fits. These are then
                # fed through the pipeline.

        saveName = 'Shift_Coords.txt'

        print shiftList

        g = []

        for entry in shiftList:

            g.append(np.hstack(entry))

        h = np.vstack(g)

        np.savetxt(saveName, h, fmt='%10.5f')

    def maskExtraPixels(self,
                        infile,
                        badpixel_dark):

        """
        Def:
        Take an input raw data file and mask extra
        bad pixels above a certain flux level.
        Although these aren't included in the cross
        correlation computation, they mess up the
        pyraf interpolation, appearing as funny blobs.
        Depending on which waveband is fed in,
        need to be careful about thermal signals on the detector.

        Input:
                infile - fits file to identify the extra bad pixels in

        Output:
                overwrite the file with a version of the
                    fits file with the bad pixels masked. It is
                    important that this is local to the shiftImage method
                     and that we are not overwriting any of the raw fits files.
        """

        # First read in the fits file and find which waveband is being used

        objTable = fits.open(infile)

        darkTable = fits.open(badpixel_dark)

        objHeader = objTable[0].header

        objHeader_one = objTable[1].header

        objHeader_two = objTable[2].header

        objHeader_three = objTable[3].header

        objFilter = objHeader['HIERARCH ESO INS FILT1 ID']

        new_obj_data = []

        new_obj_copy = []

        # Loop over the extension number

        for i in range(1, 4):

            # assign both the object data and object data copies

            objData = objTable[i].data

            darkData = darkTable[i].data

            # Mask the bad pixels we already know about before starting

            darkData[darkData == 0] = np.nan

            objCopy = copy(objData)

            objCopy = objCopy * darkData

            # First subtract thermal spectrum

            if objFilter == 'K' or objFilter == 'HK':

                print 'FOUND K or HK'

                # Need to subtract the thermal spectrum from these wavebands
                # Use the blackbody function defined beneath

                black_flux = []

                black_flux_2d = []

                # HARDWIRED K-BAND STARTING WAVELENGTH

                start_L = 1.92499995231628

                delta_L = 0.000280761742033064

                # define wavelength array

                wave_array = start_L + np.arange(0, 2048 * (delta_L), delta_L)

                # for each of the values in the wavelength
                # array evaluate the blackbody

                for i in range(len(wave_array)):

                    # append the evaluated blackbody to the black_flux array

                    black_flux.append(self.blackbody(250,
                                                     17.8,
                                                     wave_array[i] * 1E-6))

                black_flux = np.array(black_flux)

                black_flux = black_flux[::-1]

                # Create 2D detector array of the blackbody spectrum

                for i in range(2048):

                    black_flux_2d.append(black_flux)

                black_flux_2d = np.array(black_flux_2d)

                black_flux_2d = np.transpose(black_flux_2d)

                # Now have 2D thermal spectrum rising in
                # the K-band with the same dimensions as detector
                # Subtract this from the objCopy

                objCopy = objCopy - black_flux_2d

            # Now continue with additional masking as before
            # Point is to identify the extra bad pixels on the
            # object copy and then mask on the actual object and save
            # Initialise i_array and j_array for
            # storing the bad pixel coordinates

            i_array = []

            j_array = []

            # Can't loop over all pixels as this takes too long
            # Also only looking for the very worst pixels to mask

            print 'This is the percentile: %s %s %s' % \
                (np.nanpercentile(objCopy, 99.9),
                 np.nanpercentile(objCopy, 99.99),
                 np.nanpercentile(objCopy, 99.999))

            index = np.where(objCopy > (np.nanpercentile(objCopy, 99.9)))

            coords = np.array(index)

            print len(coords[0])

            # Now mask these pixels on the actual data

            for i, j in zip(coords[0], coords[1]):

                print i, j

                # Check that the indices will be in range
                if i >= 4 and i <= 2043 and j >= 4 and j <= 2043:

                    # Take the average of the surrounding
                    # pixels - if unusually high don't mask
                    # Taking a more horizontally extended
                    # chunk to check for the presence of sky lines

                    pixAv = np.nanmean([objCopy[i - 2][j],
                                        objCopy[i - 3][j],
                                        objCopy[i - 4][j],
                                        objCopy[i + 2][j],
                                        objCopy[i + 3][j],
                                        objCopy[i + 4][j],
                                        objCopy[i - 2][j + 1],
                                        objCopy[i - 2][j - 1],
                                        objCopy[i - 1][j + 1],
                                        objCopy[i - 1][j - 1],
                                        objCopy[i + 2][j + 1],
                                        objCopy[i + 2][j - 1],
                                        objCopy[i + 1][j + 1],
                                        objCopy[i + 1][j - 1]])

                    print 'PixAv and percentile: %s %s' % \
                        (pixAv, np.nanpercentile(objCopy, 95))

                    if pixAv < np.nanpercentile(objCopy, 95):

                        print 'YES - Adding to mask'

                        # Append these coordinates to an
                        # array and do all the masking at the end

                        i_array.append(i)

                        j_array.append(j)

            # Do the actual masking

            for i, j in zip(i_array, j_array):

                print 'Masking %s %s' % (i, j)

                # Mask off a cross around the offending pixel

                objData[i][j] = np.nan

                objData[i + 1][j] = np.nan

                objData[i - 1][j] = np.nan

                objData[i][j - 1] = np.nan

                objData[i][j + 1] = np.nan

                # For test also mask off the copy

                objCopy[i][j] = np.nan

                objCopy[i + 1][j] = np.nan

                objCopy[i - 1][j] = np.nan

                objCopy[i][j - 1] = np.nan

                objCopy[i][j + 1] = np.nan

            # Also mask all of the most extreme pixels

            worst_index = np.where(objCopy > 15000)

            coords = np.array(worst_index)

            for i, j in zip(coords[0], coords[1]):

                # Mask off the cross again for the worst pixels

                objData[i][j] = np.nan

                objData[i + 1][j] = np.nan

                objData[i - 1][j] = np.nan

                objData[i][j - 1] = np.nan

                objData[i][j + 1] = np.nan

            new_obj_data.append(objData)

            new_obj_copy.append(objCopy)

        temp = sys.stdout

        sys.stdout = open('log.txt', 'w')

        print objHeader
        print objHeader_one
        print objHeader_two
        print objHeader_three

        sys.stdout.close()

        sys.stdout = temp

        # Write out to a different fits file
        nameOfFile = 'masked_infile.fits'

        hdu = fits.PrimaryHDU(header=objHeader)

        hdu.writeto(nameOfFile,
                    clobber=True)

        fits.append(nameOfFile,
                    data=new_obj_data[0],
                    header=objHeader_one)

        fits.append(nameOfFile,
                    data=new_obj_data[1],
                    header=objHeader_two)

        fits.append(nameOfFile,
                    data=new_obj_data[2],
                    header=objHeader_three)

    def blackbody(self,
                  T,
                  scaling,
                  L):

        """
        Def:
        Return the blackbody function

        Input:
                T - temperature of KMOS detector
                scaling - some blackbody scaling value
                L

        Output:
                evaluated blackbody function using these parameters
        """

        # Define constants before writing the function

        h = 6.62606957E-34

        c = 299792458

        k = 1.3806488E-23

        # Initially try a temperature of 10K

        # T = 10

        preFactor = (2 * h * c ** 2) / L ** 5

        # print 'This is the pre-factor: %s' % preFactor

        expPower = (h * c) / (L * k * T)

        # print 'This is the exponent: %s' % expPower

        denom = np.exp(expPower) - 1

        # print 'This is the denominator: %s' % denom

        B_L = scaling * (preFactor / denom)

        return B_L

    def blackbodyMod(self):
        """
        Def:
        Create a blackbody model with lmfit

        Output:
                blackbody model
        """

        mod = Model(self.blackbody,
                    independent_vars=['L'],
                    param_names=['T', 'scaling'],
                    missing='drop')

        return mod

    def blackbodyModFit(self,
                        wavelength,
                        flux):

        """
        Def:
        Fit for the blackbody temperature given
        input wavelength array and flux data

        Input:
                wavelength - wavelength array
                flux - flux values at the wavelength points
        """

        mod = self.blackbodyMod()

        # Set the parameter hints from the initialPars method

        mod.set_param_hint('T', value=120, vary=False)

        mod.set_param_hint('scaling', value=1E9)

        # Initialise a parameters object to use in the fit

        fit_pars = mod.make_params()

        # Guess isn't implemented for this model

        mod_fit = mod.fit(flux, L=wavelength, params=fit_pars)

        print mod_fit.fit_report

        return mod_fit

    def quickSpecPlot(self,
                      flux,
                      wavelength):

        """
        Def:
        Helper function to quickly plot a spectrum given a
        flux and a wavelength
        Input: Flux, wavlength 1D arrays
        Output: Plot of 1D spectrum

        """

        fig, ax = plt.subplots(1, 1, figsize=(18, 10))
        ax.plot(wavelength, flux)
        ax.set_title('Object and Sky Comparison', fontsize=24)
        ax.set_xlabel(r'Wavelength$\AA$')
        ax.set_ylabel('Flux')
        ax.set_ylim(0, max(flux))
        ax.tick_params(axis='y', which='major', labelsize=15)
        plt.show()

    def compareSky(self,
                   sci_dir,
                   combNames):

        """
        Def:
        From an input sky cube and set of sci_combined file names,
        work out the median difference between the bright sky lines
        and the corresponding object pixels

        Input: skyCube - Any reconstructed skyCube only filename
               combNames - List of the sci_combined names

        Ouptut: medianVal - median difference between the sky and object values

        """

        # The combNames should be generated from within the Next routine
        namesOfFiles = np.genfromtxt(combNames, dtype='str')

        # Initialise an empty dictionary
        medVals = {}

        cubeNames = {}

        # Loop round each of the cubes in combNames, create a cube
        # and store the median
        print '[INFO:] Computing Sky Performance Statistic'

        for fileName in namesOfFiles:

            print 'Fitting %s' % fileName

            tempCube = cubeOps(sci_dir + '/' + fileName)

            IFUNR = tempCube.IFUNR

            # print IFUNR
            sky_name = sci_dir + '/combine_sci_reconstructed_arm' \
                + str(IFUNR) + '_sky.fits'

            # Create instance from the skycube
            sky_cube = cubeOps(sky_name)

            # Extract the sky flux
            sky_flux = sky_cube.centralSpec()

            wavelength = sky_cube.wave_array

            # print 'This is the sky wavelength
            # range: %s %s' % (min(wavelength), max(wavelength))
            # There is a sky_flux value at every wavelength value
            # Want to restrict the data to just the
            # inner regions of the filter. This has to be done conditionally
            # depending on what waveband we're in

            if np.nanmedian(wavelength) > 1.40 \
                    and np.nanmedian(wavelength) < 1.80:

                # we're in the H-band
                print 'H-band found'

                filter_id = 'H'

                filter_indices = np.where(np.logical_and(wavelength > 1.45,
                                                         wavelength < 1.75))[0]

                sky_flux = sky_flux[filter_indices]

                wavelength = wavelength[filter_indices]

            elif np.nanmedian(wavelength) > 2.0 \
                    and np.nanmedian(wavelength) < 2.3:

                # we're in the K-band
                print 'K-band found'

                filter_id = 'K'

                filter_indices = np.where(np.logical_and(wavelength > 2.10,
                                                         wavelength < 2.3))[0]

                sky_flux = sky_flux[filter_indices]

                wavelength = wavelength[filter_indices]

            elif np.nanmedian(wavelength) > 1.00 \
                    and np.nanmedian(wavelength) < 1.35:

                # we're in the YJ-band
                print 'YJ-band found'

                filter_id = 'YJ'

                filter_indices = np.where(np.logical_and(wavelength > 1.05,
                                                         wavelength < 1.35))[0]

                sky_flux = sky_flux[filter_indices]

                wavelength = wavelength[filter_indices]

            elif np.nanmedian(wavelength) > 1.8 \
                    and np.nanmedian(wavelength) < 2.1:

                # we're in the HK-band
                print 'HK-band found'

                filter_id = 'HK'

                filter_indices = np.where(np.logical_and(wavelength > 1.6,
                                                         wavelength < 2.3))[0]

                sky_flux = sky_flux[filter_indices]

                wavelength = wavelength[filter_indices]

            elif np.nanmedian(wavelength) > 0.8 \
                    and np.nanmedian(wavelength) < 1.05:

                # we're in the IZ-band
                print 'IZ-band found'

                filter_id = 'IZ'

                filter_indices = np.where(np.logical_and(wavelength > 0.8,
                                                         wavelength < 1.05))[0]

                sky_flux = sky_flux[filter_indices]

                wavelength = wavelength[filter_indices]

            # Could be a different wavelength, will just
            # leave the flux and wavelength unchanged right now
            # Check for where the flux exceeds a certain number of counts

            if filter_id == 'HK' or filter_id == 'YJ' or \
                    filter_id == 'H' or filter_id == 'K':

                emission_indices = np.where(sky_flux > 1000)[0]

            else:

                emission_indices = np.where(sky_flux > 50)[0]

            # Grow the emission_indices to
            # include values either side

            add_array = []

            for i in range(len(emission_indices)):

                value = emission_indices[i]

                plus_value = value + 1

                sub_value = value - 1

                if plus_value >= len(sky_flux):

                    plus_value = len(sky_flux) - 1

                if sub_value < 0:

                    sub_value = 0

                # insert these into the emission_indices array
                add_array.append(plus_value)

                add_array.append(sub_value)

            emission_indices = list(emission_indices)

            for item in add_array:

                emission_indices.append(item)

            # Now take the set of unique values from emission_indices
            new_emission_indices = np.sort(list(set(emission_indices)))

            # print 'The unique emission_indices
            # are: %s' % new_emission_indices
            # Find the sky values at these pixels
            sky_emission_lines = sky_flux[new_emission_indices]

            # Take the absolute value of the fluxes
            # to account for P-Cygni profiles
            sky_emission_lines = abs(sky_emission_lines)

            # Extract the object 1D spectrum,
            # this is from the whole cube now. Should help
            # eliminate spectrum curvature as much as possible
            object_spectrum = tempCube.total_spec[filter_indices]

            object_wavelength = tempCube.wave_array[filter_indices]

            # Find the object flux at the same pixels
            tempValues = object_spectrum[new_emission_indices]

            # Take the absolute value to account for P-Cygni profiles
            tempValues = abs(tempValues)

            # we want to normalise this by the mean object flux
            # at the wavelength values NOT contaminated by skylines
            # Find the list of unique emission_indices
            # without the sky_line emission_indices
            total_emission_indices = np.arange(0,
                                               len(object_wavelength), 1)
            contm_emission_indices = list(set(total_emission_indices)
                                          - set(new_emission_indices))

            # Take the mean value of these as the mean object flux
            tempObjContinuum = object_spectrum[contm_emission_indices]

            tempObjWavelength = object_wavelength[contm_emission_indices]

            object_continuum = abs(np.median(tempObjContinuum))

            # Now some of the entries in the object_spectrum
            # could be nan, in which case the
            # model fitting won't work. Need to
            # record the indices at which they appear and
            # remove the entries from both the object
            # spectrum and wavelength array at which they appear
            nan_list = []

            for i in range(len(tempObjContinuum)):

                if np.isnan(tempObjContinuum[i]):

                    nan_list.append(i)
            # If there is something in the nan_list delete
            # from wavelength array and object spectrum
            if nan_list:

                tempObjContinuum = [i for j,
                                    i in enumerate(tempObjContinuum)
                                    if j not in nan_list]

                tempObjWavelength = [i for j,
                                     i in enumerate(tempObjWavelength)
                                     if j not in nan_list]

            # Fit a polynomial to the continuum and subtract
            # Create the polynomial model from
            # lmFit (from lmfit import PolynomialModel)
            mod = PolynomialModel(4)

            # Have an initial guess at the model parameters
            pars = mod.guess(tempObjContinuum, x=tempObjWavelength)

            # Use the parameters for the full model fit
            out = mod.fit(tempObjContinuum, pars, x=tempObjWavelength)

            # The output of the model is the fitted continuum
            continuum = out.best_fit

            # extrapolate the model to the full wavelength range
            full_continuum = out.eval(x=object_wavelength)

            # print 'The length of the continuum is: %s' % len(continuum)
            fig, ax = plt.subplots(1, 1, figsize=(18, 10))

            ax.plot(tempObjWavelength, continuum)

            ax.plot(tempObjWavelength, tempObjContinuum)

            # plt.show()
            # Now subtract the object continuum from
            # the full object_spectrum and repeat analysis
            new_object_spectrum = object_spectrum - full_continuum

            divided_object_spectrum = object_spectrum / full_continuum

            fig, ax = plt.subplots(1, 1, figsize=(18, 10))

            ax.plot(tempCube.wave_array[filter_indices],
                    abs(new_object_spectrum))

            # plt.show()

            # print 'The object continuum value is: %s' % object_continuum

            # normalise the tempValues by this
            # norm_values = tempValues / object_continuum
            # print 'The normalised flux values are: %s ' % norm_values
            # Find the median of these norm_values
            # as the sky subtraction performance indicator
            pos_spec = np.abs(new_object_spectrum)

            std_dev = np.std(pos_spec)

            pos_spec = \
                np.nanmedian(np.array(pos_spec[np.where(pos_spec > std_dev)]))

            # print 'This is the positive object_spectrum: %s' % pos_spec
            medVals[int(tempCube.IFUNR)] = pos_spec

            cubeNames[int(tempCube.IFUNR)] = tempCube.IFUName

        # print 'The resultant dictionary is: %s' % (medVals)
        # print 'The names of the files are: %s' % cubeNames

        medVector = np.nanmean(medVals.values())

        # print 'This is the median values Dictionary: %s' % medVals
        # print medVector

        plt.close('all')

        return np.array(medVals.values()), \
            list(cubeNames.values()), np.array(medVals.keys())

    def gaussFit(self,
                 combNames):

        """
        Def:
        Uses the psfMask function from cubeClass to
        loop round the list of combNames and return a
        list of FWHM values and a list of mask profiles

        Input:
                combNames - list of sci_combined files produced
                            by the pipeline
        """

        namesOfFiles = np.genfromtxt(combNames, dtype='str')

        fwhmArray = []

        psfArray = []

        paramsArray = []

        try:

            # Loop round and create an instance of the class each time

            for fileName in namesOfFiles:

                print 'Gaussian fitting science frame: %s' % fileName

                tempCube = cubeOps(fileName)

                shift_list = [tempCube.xDit, tempCube.yDit]

                params, psfProfile, fwhm, offList = tempCube.psfMask()

                fwhmArray.append(fwhm)

                psfArray.append(psfProfile)

            return fwhmArray, psfArray, offList, params, shift_list

        except TypeError:

            # Only one file in the list of files

            print namesOfFiles

            cubeName = str(namesOfFiles)

            print 'Gaussian fitting sky frame: %s' % namesOfFiles

            tempCube = cubeOps(cubeName)

            shift_list = [tempCube.xDit, tempCube.yDit]

            params, psfProfile, fwhm, offList = tempCube.psfMask()

            return fwhm, psfProfile, offList, params, shift_list

    def combFrames(self,
                   sci_dir,
                   frame_array):

        """
        Def:
        Helper Function to combine a list of object files
        categorised by the frameCheck method. First appends
        sci_reconstructed_ to each of the names and then
        combines these together using the KMOS pipeline
        recipe kmo_combine

        Input: 
                frame_array - list of object names as specified
                                by the frameCheck method

        """

        # Remove .sof file if this exists

        if os.path.isfile('%s/sci_combine.sof' % sci_dir):

            os.system('rm %s/sci_combine.sof' % sci_dir)

        if len(frame_array) == 0:

            print 'Empty Array'

        else:

            # Create new list of arrays with prepended names

            new_names = []

            with open(sci_dir + '/sci_combine.sof', 'a') as f:

                for entry in frame_array:

                    # If the entry doesn't contain a backslash, the entry
                    # is the object name and can prepend directly

                    if entry.find("/") == -1:

                        name = sci_dir + '/' + 'sci_reconstructed_' + entry

                        f.write('%s\n' % name)

                    # Otherwise the directory structure is included and have to
                    # search for the backslash and omit up to the last one

                    else:

                        objName = entry[len(entry) - entry[::-1].find("/"):]

                        name = sci_dir + '/' + 'sci_reconstructed_' + objName

                        f.write('%s\n' % name)

            # Now execute the recipe
            os.system('esorex --output-dir=%s kmos_combine --method="header"'
                      + ' --edge_nan=TRUE %s/sci_combine.sof' % (sci_dir,
                                                                 sci_dir))

    def frameCheck(self,
                   sci_dir,
                   frameNames,
                   tracked_list):

        """
        Def:
        Loop round a list of frameNames of a given type
        and apply the science reduction to each pair. Return the vectors
        which contain the sky tweak performance for each pair and
        for each IFU in each pair, and the vector
        containing the fwhm of the tracked star. The tracked star is defined in
        the multiExtract method.
        Also bin objects based on their fwhm and return
        a dictionary containing the different bins.
        Input: skyCube - a given reconstructed skyCube.
        Note this assumed that the
        sky values will not vary significantly
        in wavelength from frame to frame.
        Need to check this.
        Input:
                sci_dir - path to the current science directory
                frameNames - file containing a list of object and sky frames
                            in the standard format, with column 1
                            the file name and column 2 the file type
                tracked_list - list of the names of three stars tracked
                                by KMOS, 1 star per detector. If this isn't
                                available write the name of one star three
                                times
        Output:
                IFUValVec - mean for each IFU of skytweak performance
                frameValVec - mean for each frame of skytweak performance
                fwhmValVec - tracked star fwhm for each frame
                fwhmDict - keys of 'Best', 'Good', 'Okay', 'Bad' corresponding
                to the fwhm of the tracked star, and values being the name of
                the files which fall into each of these bins

        Uses: self.compareSky
              self.gaussFit
              cubeOps

        """
        # Initially look for the sci_reduc.sof file

        if not os.path.isfile('sci_reduc.sof'):

            raise ValueError('Cannot find sci_reduc.sof file')

        # First read in the names and the types from frameNames

        data = np.genfromtxt(frameNames, dtype='str')

        # Save the names and types as lists

        names = data[0:, 0]

        types = data[0:, 1]

        # Set the sci_comb names defined in the cubeclass

        if types[1] == 'O':

            combNames = cubeOps(names[1]).combNames

            rec_combNames = cubeOps(names[1]).rec_combNames

        elif types[2] == 'O':

            combNames = cubeOps(names[2]).combNames

            rec_combNames = cubeOps(names[2]).rec_combNames

        elif types[3] == 'O':

            combNames = cubeOps(names[3]).combNames

            rec_combNames = cubeOps(names[3]).rec_combNames

        else:

            print 'Having difficulty setting sci_comb names'

        print 'This is the order of the combined names: %s' % combNames

        # Define the tracked star ID
        # As input to the method is a list of three names
        # these are for the tracked stars on each detector
        # if there is only a single tracked star then the names
        # should all be the same and everything will work

        track_name_one = tracked_list[0]

        track_name_two = tracked_list[1]

        track_name_three = tracked_list[2]

        # Loop round the list of combNames until each tracked_name appears
        # these are used later on in the method to find the pixel
        # scale of that standard star

        for entry in combNames:

            if entry.find(track_name_one) != -1:

                tracked_star_one = sci_dir + '/' + entry

        for entry in combNames:

            if entry.find(track_name_two) != -1:

                tracked_star_two = sci_dir + '/' + entry

        for entry in combNames:

            if entry.find(track_name_three) != -1:

                tracked_star_three = sci_dir + '/' + entry

        # Initialise loop counter and variables

        counter = 0

        frameValVec = []

        IFUValVec = []

        namesVec = []

        fwhmValVec_one = []
        fwhmValVec_two = []
        fwhmValVec_three = []

        # Set up the different bins for fwhm of tracked star

        a_fwhm_names = []

        b_fwhm_names = []

        c_fwhm_names = []

        d_fwhm_names = []

        # Remove the combine input file if it exists

        if os.path.isfile('combine_input.txt'):

            os.system('rm combine_input.txt')

        # Loop around each name, assign sky pair and populate
        # the sky tweak performance and fwhm variables

        for i in range(1, len(names)):

            if types[i] == 'O':

                counter += 1

                objFile = names[i]

                if i == 1:

                    skyFile = names[i + 1]

                elif types[i - 1] == 'S':

                    skyFile = names[i - 1]

                else:

                    skyFile = names[i + 1]

                print '[INFO]: reducing file: %s : ' % objFile

                namesVec.append(objFile)

                # Create the new .sof file at each stage
                # for just the object sky pair
                # the .sof file must be in the directory,
                # and must have all the other
                # files required already specified.

                # Create a copy of the sci_reduc.sof in a new temporary file

                with open('sci_reduc.sof') as f:

                    with open('sci_reduc_temp.sof', 'w') as f1:

                        for line in f:

                                f1.write(line)

                # Append the current object and skyfile
                # names to the newly created .sof file

                with open('sci_reduc_temp.sof', 'a') as f:

                    f.write('\n%s\tSCIENCE' % objFile)

                    f.write('\n%s\tSCIENCE' % skyFile)

                # Now just execute the esorex recipe for this new file
                os.system('esorex --output-dir=%s kmos_sci_red' % sci_dir
                          + '  --pix_scale=0.1 --oscan=FALSE --sky_tweak=TRUE'
                          + ' --edge_nan=TRUE sci_reduc_temp.sof')

                # We have all the science products now
                # execute the above method for each
                # of the pairs. Should think of a better
                # way to create the combNames file

                print 'Checking IFU sky tweak performance'

                # This is the array of IFU values for each frame

                medVals, cube_name_list, IFU_number_array = \
                    self.compareSky(sci_dir, combNames)

                # print the value of this vector

                print 'The sky statistic values as a' \
                    + 'function of IFU are: %s' % medVals

                print 'The median value of this array' \
                    + ' is: %s' % np.nanmedian(medVals)

                # Append the full vector for a more detailed plot

                IFUValVec.append(medVals)

                # Append the mean value for the average plot

                frameValVec.append(np.nanmedian(medVals))

                # Move onto FWHM analysis of chosen tracked stars
                # The gaussian fitting will be executed three times
                # once for each of the tracked stars

                print 'Checking PSF of tracked stars'

                fwhm_one, psfProfile_one, offList_one, params_one, \
                    shift_list_one = self.gaussFit('tracked_one.txt')

                fwhmValVec_one.append(fwhm_one)

                fwhm_two, psfProfile_two, offList_two, params_two, \
                    shift_list_two = self.gaussFit('tracked_two.txt')

                fwhmValVec_two.append(fwhm_two)

                fwhm_three, psfProfile_three, offList_three, params_three, \
                    shift_list_three = self.gaussFit('tracked_three.txt')

                fwhmValVec_three.append(fwhm_three)

                # remove the temporary .sof file and
                # go back to the start of the loop

                os.system('rm sci_reduc_temp.sof')

                # Change FWHM to arcsecond scale, by using the pixel scale
                # recover the pixel scale by creating a cubeClass instance

                pixel_scale = float(cubeOps(tracked_star_one).pix_scale)

                # define the three fwhm values from the stars

                arc_fwhm_one = fwhm_one * pixel_scale

                arc_fwhm_two = fwhm_two * pixel_scale

                arc_fwhm_three = fwhm_three * pixel_scale

                arc_fwhm = np.mean([arc_fwhm_one,
                                    arc_fwhm_two,
                                    arc_fwhm_three])

                # Now since we are tracking three stars need to create arrays
                # for writing to the file. Tricky as some of the arms aren't
                # operational and this must be taken into account. Constructing
                # arrays for the arc_fwhm, shift_list and params.

                first = np.arange(1, 9, 1)

                second = np.arange(9, 17, 1)

                third = np.arange(17, 25, 1)

                # loop round the offList and reduce the size of
                # first second and third if they contain IFU numbers
                # which are in the offList

                for item in offList_one:
                    if item in first:
                        first = first[1:]

                for item in offList_one:
                    if item in second:
                        second = second[1:]

                for item in offList_one:
                    if item in third:
                        third = third[1:]

                # now that the arrays are of the correct size
                # create the overall arrays for each parameter

                # the fwhm of the stars

                arc_fwhm_one_array = np.repeat(arc_fwhm_one,
                                               len(first))

                arc_fwhm_two_array = np.repeat(arc_fwhm_two,
                                               len(second))

                arc_fwhm_three_array = np.repeat(arc_fwhm_three,
                                                 len(third))

                arc_fwhm_array = np.hstack([arc_fwhm_one_array,
                                           arc_fwhm_two_array,
                                           arc_fwhm_three_array])

                print 'CHECKING FOR EQUAL LENGTHS'
                print 'FWHM array has length %s' % len(arc_fwhm_array)
                print 'cubelist has length %s' % len(cube_name_list)

                # the x center of the gaussian fit

                x_center_one_array = np.repeat(params_one['center_x'],
                                               len(first))

                x_center_two_array = np.repeat(params_two['center_x'],
                                               len(second))

                x_center_three_array = np.repeat(params_three['center_x'],
                                                 len(third))

                x_center_array = np.hstack([x_center_one_array,
                                           x_center_two_array,
                                           x_center_three_array])

                # the y center of the gaussian fit

                y_center_one_array = np.repeat(params_one['center_y'],
                                               len(first))

                y_center_two_array = np.repeat(params_two['center_y'],
                                               len(second))

                y_center_three_array = np.repeat(params_three['center_y'],
                                                 len(third))

                y_center_array = np.hstack([y_center_one_array,
                                           y_center_two_array,
                                           y_center_three_array])

                # the x dither value

                x_dither_one_array = np.repeat(shift_list_one[0],
                                               len(first))

                x_dither_two_array = np.repeat(shift_list_two[0],
                                               len(second))

                x_dither_three_array = np.repeat(shift_list_three[0],
                                                 len(third))

                x_dither_array = np.hstack([x_dither_one_array,
                                           x_dither_two_array,
                                           x_dither_three_array])

                # the y dither value

                y_dither_one_array = np.repeat(shift_list_one[1],
                                               len(first))

                y_dither_two_array = np.repeat(shift_list_two[1],
                                               len(second))

                y_dither_three_array = np.repeat(shift_list_three[1],
                                                 len(third))

                y_dither_array = np.hstack([y_dither_one_array,
                                           y_dither_two_array,
                                           y_dither_three_array])

                # Define the reconstructed objFile name
                # If the entry doesn't contain a backslash, the entry
                # is the object name and can prepend directly

                if objFile.find("/") == -1:
                    rec_objFile = sci_dir + '/' + \
                        'sci_reconstructed_' + objFile

                # Otherwise the directory structure is included and have to
                # search for the backslash and omit up to the last one

                else:

                    objName = objFile[len(objFile) - objFile[::-1].find("/"):]

                    rec_objFile = sci_dir + '/' + \
                        'sci_reconstructed_' + objName

                # Make the output file for combining
                # - have all the information now
                # For each object

                with open('combine_input.txt', 'a') as f:

                    for i in range(len(cube_name_list)):

                        f.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %
                                (rec_objFile,
                                 skyFile,
                                 IFU_number_array[i],
                                 cube_name_list[i],
                                 medVals[i],
                                 arc_fwhm_array[i],
                                 x_center_array[i],
                                 y_center_array[i],
                                 x_dither_array[i] / pixel_scale,
                                 y_dither_array[i] / pixel_scale))

                    if (arc_fwhm > 0.0 and arc_fwhm < 0.6):

                        print '[INFO]: Placing object in best bin'

                        a_fwhm_names.append(objFile)

                    elif (arc_fwhm > 0.6 and arc_fwhm < 1.0):

                        print '[INFO]: Placing object in good bin'

                        b_fwhm_names.append(objFile)

                    elif (arc_fwhm > 1.0 and arc_fwhm < 1.5):

                        print '[INFO]: Placing object in okay bin'

                        c_fwhm_names.append(objFile)

                    else:

                        print '[INFO]: Placing object in bad bin'

                        d_fwhm_names.append(objFile)

        # Should now have populated the frameValVec,
        # IFUValVec and incremented counter
        # Only for the objects that have passed the skytweak performance test

        IFUValVec = np.array(IFUValVec)

        fwhmValVec_one = np.array(fwhmValVec_one)
        fwhmValVec_two = np.array(fwhmValVec_two)
        fwhmValVec_three = np.array(fwhmValVec_three)

        offList = np.array(offList_one)

        # Convert the FWHM to arcseconds instead of pixels

        fwhmValVec_one = pixel_scale * np.array(fwhmValVec_one)
        fwhmValVec_two = pixel_scale * np.array(fwhmValVec_two)
        fwhmValVec_three = pixel_scale * np.array(fwhmValVec_three)

        ID = np.arange(0.0, counter, 1.0)

        # Make the dictionary of fwhm values
        fwhmDict = {'Best': a_fwhm_names,
                    'Good': b_fwhm_names,
                    'Okay': c_fwhm_names,
                    'Bad': d_fwhm_names}

        # Shift all plots to a plots folder within the science directory
        plot_dir_name = sci_dir + '/Plots'

        if os.path.isdir(plot_dir_name):

            os.system('rm -rf %s' % plot_dir_name)

        os.system('mkdir %s' % plot_dir_name)

        # Move all of the newly created Shifted files into this directory

        os.system('mv %s %s' % (sci_dir + '/*.png', plot_dir_name))

        # Return all of these values

        return ID, offList, namesVec, IFUValVec, \
            frameValVec, fwhmValVec_one, \
            fwhmValVec_two, fwhmValVec_three, fwhmDict

    def meanIFUPlot(self,
                    offList,
                    namesVec,
                    IFUValVec):

        """
        Def:
        Takes the output from frameCheck and plots a line graph
        of skytweak performance against IFUID for each input frame.
        Each of the IFUs which are not operational are plotted as
        np.nan

        Input:
                ID - Vector from 1 - len(number of frames)
                offList - List of the IFUs which are not operational
                namesVec - The names of the input object files
                IFUValVec - Vector of means for each IFU

        Output - Plot of performance against IFU, recognising the
                    IFUs which are not illuminated
        """

        # Construct ID array of length 24

        IFUID = np.arange(1.0, 25, 1.0)

        # Insert np.nan at the locations where the IFU is off
        # Initialise the counter for the frame naming

        val = 0

        colors_plot = cycle(cm.rainbow(np.linspace(0, 1, len(IFUValVec))))

        colors_scatter = cycle(cm.rainbow(np.linspace(0, 1, len(IFUValVec))))

        # Collape all the information at the IFU level onto a single plot

        fig, ax = plt.subplots(1, 1, figsize=(14.0, 14.0))

        for entry in IFUValVec:

            # Extend the value array to match the IFUID array

            for value in offList:

                entry = np.insert(entry, value, np.nan)

            ax.plot(IFUID,
                    entry,
                    label=namesVec[val],
                    color=next(colors_plot))

            ax.scatter(IFUID,
                       entry,
                       color=next(colors_scatter))

            val += 1

        ax.set_title('Sky Tweak Performance vs. IFU ID')
        ax.set_xlabel('IFU ID')
        ax.set_xlim(1, 24)
        ax.set_xticks((np.arange(min(IFUID), max(IFUID) + 1, 1.0)))
        ax.grid(b=True, which='both', linestyle='--')

        # plt.legend(prop={'size':10})

        fig.savefig('IFU_by_Frame.png')

        # plt.show()

        plt.close('all')

    def indIFUPlot(self,
                   offList,
                   ID,
                   IFUValVec):

        """
        Def: Uses the output of frameCheck to
        Plot for each individual IFU the performance of skytweak
        against frame ID. More detail as to how well skytweak is
        performing.

        Input:
                - Offlist: List of IFU numbers which aren't illuminated
                - ID: np.arange between 0 and total number of operational IFUs
                - IFUValVec: 2D array of median sky tweak performance values

        Output:
                - subplot array showing how well the sky has been subtracted
                    in each individual IFU
        """

        # Make a plot for each IFU in a subplot array

        fig, axArray = plt.subplots(3, 8, figsize=(20.0, 15.0))

        IFUCount = 0

        dataCount = 0

        # Have the data now - populate the subplots

        for col in range(3):

            for row in range(8):

                # Only plot if the IFU is functioning

                if IFUCount not in offList:

                    frameVec = IFUValVec[:, dataCount]

                    axArray[col][row].plot(ID, frameVec)
                    axArray[col][row].scatter(ID, frameVec)
                    axArray[col][row].set_title('IFU %s' % (IFUCount + 1))
                    axArray[col][row].set_xlabel('Frame ID')
                    axArray[col][row].set_xticks((np.arange(min(ID),
                                                            max(ID) + 1, 1.0)))
                    axArray[col][row].grid(b=True,
                                           which='both',
                                           linestyle='--')
                    axArray[col][row].set_xlim(0, len(ID))

                    dataCount += 1

                # Increment the IFUCount number

                IFUCount += 1

        # Subplots populated, save the overall figure

        fig.savefig('IFU_subplots.png')

        # plt.show()

        plt.close('all')

    def multiIndIFUPlot(self,
                        offList,
                        ID,
                        IFUValVec1,
                        IFUValVec2):

        """
        Def: Uses the output of frameCheck to
        Plot for each individual IFU the performance of skytweak
        against frame ID. More detail as to how well skytweak is
        performing.

        Input:
                - Offlist: List of IFU numbers which aren't illuminated
                - ID: np.arange between 0 and total number of operational IFUs
                - IFUValVec: 2D array of median sky tweak performance values

        Output:
                - subplot array showing how well the sky has been subtracted
                    in each individual IFU
        """

        # Make a plot for each IFU in a subplot array

        fig, axArray = plt.subplots(3, 8, figsize=(20.0, 15.0))

        IFUCount = 0

        dataCount = 0

        # Have the data now - populate the subplots

        for col in range(3):

            for row in range(8):

                # Only plot if the IFU is functioning

                if IFUCount not in offList:

                    frameVec1 = IFUValVec1[:, dataCount]

                    frameVec2 = IFUValVec2[:, dataCount]

                    axArray[col][row].plot(ID, frameVec1, color='blue')
                    axArray[col][row].scatter(ID, frameVec1, color='blue')
                    axArray[col][row].plot(ID, frameVec2, color='red')
                    axArray[col][row].scatter(ID, frameVec2, color='red')
                    axArray[col][row].set_title('IFU %s' % (IFUCount + 1))
                    axArray[col][row].set_xlabel('Frame ID')
                    axArray[col][row].set_xticks((np.arange(min(ID),
                                                            max(ID) + 1, 1.0)))
                    axArray[col][row].grid(b=True,
                                           which='both',
                                           linestyle='--')
                    axArray[col][row].set_xlim(0, len(ID))

                    dataCount += 1

                # Increment the IFUCount number

                IFUCount += 1

        # Subplots populated, save the overall figure

        fig.savefig('IFU_subplots_double.png')

        # plt.show()

        plt.close('all')

    def meanFramePlot(self,
                      ID,
                      frameValVec):

        """
        Def:
        Uses the output from frameCheck to make a simple
        plot of the mean sky tweak performance against frame

        Input:
                - ID np.arange between 0 and count of number of frames
                - frameValVec: 1D array of mean sky tweak performance
        """

        # Make the overall mean plot of performance for the frames

        # Create a figure and plot the results

        fig, ax = plt.subplots(1, 1, figsize=(14.0, 14.0))

        ax.plot(ID, frameValVec)

        ax.scatter(ID, frameValVec)

        ax.set_title('Sky Tweak Performance vs. Frame ID')

        ax.set_xlabel('Frame ID')

        ax.set_xticks((np.arange(min(ID), max(ID) + 1, 1.0)))

        ax.set_xlim(0, len(ID))

        ax.grid(b=True, which='both', linestyle='--')

        fig.savefig('frame_performance.png')

        # plt.show()

        plt.close('all')

    def multiMeanFramePlot(self,
                           ID,
                           frameValVec1,
                           frameValVec2):

        """
        Def:
        Uses the output from frameCheck to make a simple
        plot of the mean sky tweak performance against frame
        Input:
                - ID: np.arange between 0 and count of number of frames
                - frameValVec: 1D array of mean sky tweak performance
        """

        # Make the overall mean plot of performance for the frames
        # Create a figure and plot the results

        fig, ax = plt.subplots(1, 1, figsize=(14.0, 14.0))

        ax.plot(ID, frameValVec1, color='blue')
        ax.scatter(ID, frameValVec1, color='blue')
        ax.plot(ID, frameValVec2, color='red')
        ax.scatter(ID, frameValVec2, color='red')
        ax.set_title('Sky Tweak Performance vs. Frame ID')
        ax.set_xlabel('Frame ID')
        ax.set_xlim(0, len(ID))
        ax.set_xticks((np.arange(min(ID), max(ID)+ 1, 1.0)))
        ax.grid(b=True, which='both', linestyle='--')

        fig.savefig('frame_performance_double.png')

        # plt.show()

        plt.close('all')

    def meanFWHMPlot(self,
                     ID,
                     fwhmValVec,
                     name):

        """
        Def:
        Uses the output from frameCheck to make a simple
        plot of the tracked star FWHM against frame ID

        Input:
                - ID: np.arange between 0 and count of number of frames
                - fwhmValVec: 1D array of tracked star fwhm in each frame
        """

        fig, ax = plt.subplots(1, 1, figsize=(14.0, 14.0))

        ax.plot(ID, fwhmValVec)
        ax.scatter(ID, fwhmValVec)
        ax.set_title('Average fwhm vs. Frame ID')
        ax.set_xlabel('Frame ID')
        ax.set_xlim(0, len(ID))
        ax.set_xticks((np.arange(min(ID), max(ID) + 1, 1.0)))
        ax.grid(b=True, which='both', linestyle='--')

        save_name = 'frame_fwhm' + name + '.png'

        fig.savefig(save_name)

        # plt.show()

        plt.close('all')

    def multiMeanFWHMPlot(self,
                          ID,
                          fwhmValVec1,
                          fwhmValVec2,
                          name):

        """
        Def:
        Uses the output from frameCheck to make a simple
        plot of the tracked star FWHM against frame ID
        Input - ID: np.arange between 0 and count of number of frames
              - fwhmValVec: 1D array of tracked star fwhm in each frame
        """

        fig, ax = plt.subplots(1, 1, figsize=(14.0, 14.0))

        ax.plot(ID, fwhmValVec1, color='blue')
        ax.scatter(ID, fwhmValVec1, color='blue')
        ax.plot(ID, fwhmValVec2, color='red')
        ax.scatter(ID, fwhmValVec2, color='red')
        ax.set_xlim(0, len(ID))
        ax.set_title('Average fwhm vs. Frame ID')
        ax.set_xlabel('Frame ID')
        ax.set_xticks((np.arange(min(ID), max(ID) + 1, 1.0)))
        ax.grid(b=True, which='both', linestyle='--')

        save_name = 'fwhm_double_' + name + '.png' 
        fig.savefig(save_name)

        # plt.show()

        plt.close('all')

    def extractSpec(self,
                    sci_dir,
                    fwhmDict,
                    combNames,
                    rec_combNames,
                    tracked_name):

        """
        Def:
        Takes the grouped tracked star fwhm dictionary, combines
        the objects in each bin using the ESO pipeline and then
        extracts the optimal spectrum from each IFU in each bin
        with appended dictionary name, i.e. 'Good_sci_combined*'

        Input:
                sci_dir - current science directory
                fwhmDict - dictionary of fwhm values for the objects
                combNames - combined cube names for the objects
                rec_combNames - reconstructed combined names
                tracked_name - name of the tracked star

        Output:

        """
        # Writing out a temporary file containing the tracked star name
        # Can probably in the future get the name of
        # the IFU tracking a standard
        # Star straight from the fits header

        track_name = tracked_name

        # Loop round the list of combNames until the track_name appears

        for entry in combNames:

            if entry.find(track_name) != -1:

                tracked_star = sci_dir + '/' + entry

        # First check the directory for the
        # rec_combNames and delete if they exist

        for entry in rec_combNames:

            if os.path.isfile('%s/%s' % (sci_dir, entry)):

                os.system('rm %s/%s' % (sci_dir, entry))

            # Set the rec_combName of the
            # standard star in the same way as above

            if entry.find(track_name) != -1:

                rec_tracked_star = sci_dir + '/' + entry

        # retrieve the dictionary combining the science names and IFU numbers

        combDict = cubeOps(tracked_star).combDict

        # Remove the current sci_combined*.fits prior to this analysis

        os.system('rm %s/sci_combined*.fits' % sci_dir)

        # Initialise dictionary for final fwhm values

        fwhm_values = {}

        # Loop around each of the keys in the FWHM dictionary

        for group in fwhmDict.keys():

            print '[INFO]: Extracting Spectra for the %s group' % group

            if fwhmDict[group]:

                # The case with only 1 entry in the group (complex)

                if len(fwhmDict[group]) == 1:

                    print '[INFO]: Only 1 '
                    + '%s PSF frame: Selecting %s PSF Cubes' % (group,
                                                                group)

                    # Construct the reconstructed file name
                    # If the entry doesn't contain a backslash, the entry
                    # is the object name and can prepend directly

                    if fwhmDict[group][0].find("/") == -1:

                        rec_frame = sci_dir + '/sci_reconstructed_' \
                            + fwhmDict[group][0]

                    # Otherwise the directory structure
                    # is included and have to
                    # search for the backslash and
                    # omit up to the last character

                    else:

                        objName = fwhmDict[group][0][len(fwhmDict[group][0])
                                                     - fwhmDict[group][0][::-1]
                                                     .find("/"):]

                        rec_frame = sci_dir + '/sci_reconstructed_' + objName

                    # Now have the correct name of the reconstructed file
                    # Need to select the correct extension for the tracked_star

                    ext = combDict[track_name]

                    # Open the reconstructed frame and
                    # assign the data for the correct IFU:

                    Table = fits.open(rec_frame)

                    primHeader = Table[0].header

                    dataHeader = Table[ext].header

                    cube_data = Table[ext].data

                    # Fix the issue with the fits header

                    temp = sys.stdout

                    sys.stdout = open('log.txt', 'w')

                    print (primHeader)
                    print (dataHeader)

                    sys.stdout.close()

                    sys.stdout = temp

                    os.system('rm log.txt')

                    # Write out to a new fits file

                    objhdu = fits.PrimaryHDU(header=primHeader)

                    objhdu.writeto(tracked_star,
                                   clobber=True)

                    fits.append(tracked_star,
                                data=cube_data,
                                header=dataHeader)

                    # Now on the same footing as below

                    tracked_cube = cubeOps(tracked_star)

                    # Extract the PSFProfile from the stacked, tracked star

                    params, psfProfile, FWHM, offList = tracked_cube.psfMask()

                    # Add best stacked FWHM to the dictionary
                    # just to check afterwards

                    fwhm_values[group] = FWHM * float(tracked_cube.pix_scale)

                    # Set the rec_combNames as an iterable
                    # for specifying the file name

                    iterCombNames_one = cycle(combNames)

                    iterCombNames_two = cycle(combNames)

                    # Now get all the operational IFUs and do the same

                    for name in combDict.keys():

                        ext = combDict[name]

                        # Open the reconstructed frame:

                        Table = fits.open(rec_frame)

                        primHeader = Table[0].header

                        dataHeader = Table[ext].header

                        cube_data = Table[ext].data

                        temp = sys.stdout

                        sys.stdout = open('log.txt', 'w')

                        print (primHeader)
                        print (dataHeader)

                        sys.stdout.close()

                        sys.stdout = temp

                        os.system('rm log.txt')

                        # Write out to a new fits file
                        objhdu = fits.PrimaryHDU(header=primHeader)

                        objhdu.writeto(sci_dir + '/' + next(iterCombNames_one),
                                       clobber=True)
                        fits.append(sci_dir + '/' + next(iterCombNames_two),
                                    data=cube_data,
                                    header=dataHeader)

                    new_name_vec = []

                    # Append best to all the rec_combNames

                    for entry in combNames:

                        group_name = sci_dir + '/' + group + '_' + entry

                        new_name_vec.append(group_name)

                        os.system('mv %s %s' % (sci_dir + '/' + entry,
                                                group_name))

                # The Case where there is more
                # than one entry in the group (normal)

                elif len(fwhmDict[group]) > 1:

                    print 'Combining %s PSF Cubes' % group

                    self.combFrames(sci_dir, fwhmDict[group])

                    # The output from this is the
                    # combined_sci file names contained in rec_combNames
                    # These are all stacked data cubes in bins of seeing

                    tracked_cube = cubeOps(rec_tracked_star)

                    # Extract the PSFProfile from the stacked, tracked star

                    params, psfProfile, FWHM, offList = tracked_cube.psfMask()

                    # Add best stacked FWHM to the
                    # dictionary just to check afterwards

                    fwhm_values[group] = FWHM * float(tracked_cube.pix_scale)

                    new_name_vec = []

                    # Append best to all the rec_combNames

                    for entry in zip(rec_combNames, combNames):

                        current_name = sci_dir + '/' + entry[0]

                        group_name = sci_dir + '/' + group + '_' + entry[1]

                        new_name_vec.append(group_name)

                        os.system('mv %s %s' % (current_name, group_name))

                # Should now have the working directory filled with objects
                # That have the same name, regardless
                # of whether they passed through
                # the single object in a bin route or
                # the multi-route, with these names
                # stored in the new_name_vec array

                print 'Fitting Gaussian to each of' \
                      + 'the stacked %s objects' % group

                # Record the centre of the tracked star

                tracked_centre = [params['center_y'], params['center_x']]

                tracked_profile = copy(psfProfile)

                tracked_fwhm = copy(FWHM)

                # From here onwards should definitely
                # be fine with the new sci_dir convention
                # First Fit a gaussian to each of
                # the objects to determine the center!

                spectra_names = []

                for name in new_name_vec:

                    spec_name = name[:-5] + '_spectrum.fits'

                    spectra_names.append(spec_name)

                    cube = cubeOps(name)

                    # Find the central value of the object flux by
                    # fitting a gaussian to the image

                    params, objProfile, FWHM, offList = cube.psfMask()

                    obj_centre = [params['center_y'], params['center_x']]

                    print '[INFO]: The standard star' \
                          + ' centre is: %s' % tracked_centre

                    print '[INFO]: The Object centre is: %s' % obj_centre

                    # Find the difference between
                    # the tracked centre and obj centre

                    x_shift = obj_centre[0] - tracked_centre[0]

                    y_shift = obj_centre[1] - tracked_centre[1]

                    x_shift = int(np.round(x_shift))

                    y_shift = int(np.round(y_shift))

                    print '[INFO]: Shifting Profile' \
                          + ' by: %s %s' % (x_shift, y_shift)

                    # Use numpy.roll to shift the psfMask
                    # to the location of the object

                    new_mask = np.roll(tracked_profile, y_shift, axis=0)

                    # For the x_shift need to loop
                    # round the elements of the new_mask

                    final_new_mask = []

                    for arr in new_mask:

                        final_new_mask.append(np.roll(arr, x_shift))

                    final_new_mask = np.array(final_new_mask)

                    # Check to see that the gaussian and shifted profile align

                    colFig, colAx = plt.subplots(1, 1, figsize=(14.0, 14.0))

                    colCax = colAx.imshow(final_new_mask,
                                          interpolation='bicubic')

                    colFig.colorbar(colCax)

                    # Extract each optimal spectrum

                    optimal_spec = \
                        cube.optimalSpecFromProfile(final_new_mask,
                                                    tracked_fwhm,
                                                    params['center_y'],
                                                    params['center_x'])

                    # Save the optimal spectrum for each object
                    # Need to create a new fits table for this

                    tbhdu = \
                        fits.new_table(
                            fits.ColDefs([fits.Column(name='Wavelength',
                                                      format='E',
                                                      array=cube.wave_array),
                                          fits.Column(name='Flux',
                                                      format='E',
                                                      array=optimal_spec)]))

                    prihdu = fits.PrimaryHDU(header=cube.primHeader)

                    thdulist = fits.HDUList([prihdu, tbhdu])

                    thdulist.writeto(spec_name, clobber=True)

                    plot_sky_name = sci_dir + '/combine_sci_reconstructed_arm' \
                        + str(cube.IFUNR) + '_sky.fits'

                    print 'The skycube name for the' \
                          + ' plotting routine is: %s' % plot_sky_name

                    self.plotSpecs(sci_dir, spec_name, plot_sky_name, 1)

                # Create a new sub-directory in the Science
                # directory to house the spectra for this grouP

                new_dir_name = sci_dir + '/' + group

                if os.path.isdir(new_dir_name):

                    os.system('rm -rf %s' % new_dir_name)

                os.system('mkdir %s' % new_dir_name)

                # Move all of the newly created group
                # objects into this directory

                os.system('mv %s %s' % (sci_dir + '/' + group + '*',
                                        new_dir_name))

        print '[INFO]: %s' % fwhm_values

    def reduce_list_seeing(self,
                           combine_file,
                           seeing_lower,
                           seeing_upper):

        """
        Def:
        Helper method for reducing the combine_input.txt list of files
        down to the chosen seeing limits. Takes the list and returns a
        shorter list containing only objects appearing within the defined
        seeing limits.

        Input:
                combine_file - one of the outputs from frameCheck
                seeing_lower - lower limit for the seeing
                seeing_upper - upper limit for the seeing

        Output:
                    new_Table - shortened version of the list containing
                                the names and properties of objects to combine
        """

        # Read in the combine_file, which contains 4 columns

        Table = np.loadtxt(combine_file, dtype='str')

        # zip these together into a combined object, which can be looped

        zipped_entries = zip(Table[:, 0],
                             Table[:, 1],
                             Table[:, 2],
                             Table[:, 3],
                             Table[:, 4],
                             Table[:, 5],
                             Table[:, 6],
                             Table[:, 7],
                             Table[:, 8],
                             Table[:, 9])

        # Loop over each of the entries and
        # decide if it is in this seeing range
        # The seeing is the last entry of the table

        new_Table = []

        for row in zipped_entries:

            if float(row[5]) > seeing_lower and float(row[5]) < seeing_upper:

                new_Table.append(row)

        # Return the newly generated and reduced table

        return new_Table

    def reduce_list_name(self,
                         combine_list,
                         ifu_name):

        """
        Def:
        Helper method for reducing the output from reduce_list_seeing
        down to the chosen object name. Takes the list and returns a
        shorter list containing only objects appearing within the name defined
        within a loop.

        Input:
                 combine_list - output from reduce_list_seeing
                 ifu_name - IFU name of the object to be combined

        Output: 
                new_Table - shortened version of the list containing
                            the names and properties of objects to combine
        """

        new_Table = []

        # loop round the entries in the combine_list

        for row in combine_list:

            # The column containing the names is index number 2

            if row[3] == ifu_name:

                new_Table.append(row)

        return new_Table

    def reduce_list_sky(self,
                        combine_list,
                        performance_limit):

        """
        Def:
        Helper method for reducing the output from reduce_list_name
        down to only those which pass the sky subtraction test.
        Takes the list and returns a shorter list containing only
        objects with fourth index < performance_limit.

        Input:
                combine_list - output from reduce_list_seeing
                performance_limit - The sky subtraction statistic value

        Output:
                new_Table - shortened version of the list containing
                            the names and properties of objects to combine
        """

        new_Table = []

        # loop round the entries in the combine_list

        for row in combine_list:

            print row[4]

            # The column containing the names is index number 2

            if float(row[4]) < performance_limit:

                new_Table.append(row)

        return new_Table

    def compute_shifts(self,
                       sci_dir,
                       combine_list,
                       name):

        """
        Def:
        Helper method for taking the output from
        reduce_list_sky and calculating the list of
        shift values to pass to combine_by_name.
        This function will be called each time for
        the different object names.

        Input:
                sci_dir - path to the science directory
                combine_list - output from reduce_list_seeing
                name - the name of the object being shifted
        Output:
                *_shift_file.txt - List of shifts relative to first in the list
                plot of the increase in the difference
                between actual and theoretical
                shift values - shift_plot.png

        """

        # Record the names of the deltaFile and shiftFile

        deltaFile = sci_dir + '/' + name + '_delta_file.txt'

        shiftFile = sci_dir + '/' + name + '_shift_file.txt'

        # If the shift file exists in the current directory, delete this

        if os.path.isfile(deltaFile):

            os.system('rm %s' % deltaFile)

        if os.path.isfile(shiftFile):

            os.system('rm %s' % shiftFile)

        # print combine_list[0]

        # Find the initial central values

        center_x = float(combine_list[0][7])

        center_y = float(combine_list[0][6])

        shift_x = float(combine_list[0][8])

        shift_y = float(combine_list[0][9])

        # Now write to the shift file the actual shifts

        # write to the delta file the difference between expected and actual

        with open(shiftFile, 'a') as s:

            with open(deltaFile, 'a') as d:

                for i in range(1, len(combine_list)):

                    row = combine_list[i]

                    real_shift_x = float(row[7]) - center_x

                    real_shift_y = center_y - float(row[6])

                    dshift_x = shift_x - float(row[8])

                    dshift_y = shift_y - float(row[9])

                    dX = real_shift_x + dshift_x

                    dY = real_shift_y + dshift_y

                    s.write('%s\t%s\n' % (real_shift_x, real_shift_y))

                    d.write('%s\t%s\n' % (dX, dY))

        # shift file has now been created for each of the objects

    def combine_by_name(self,
                        sci_dir,
                        combine_file,
                        seeing_lower,
                        seeing_upper,
                        performance_limit):

        """
        Def:
        Main method for combining the frames into
        the science cubes. Right now will not do anything
        special with the produced output files,
        will just leave all of those in the science directory.
        Selects all of the frames which pass the
        sky subtraction test and are within a particular seeing
        range for each of the objects. More powerful way of combining.

        Input:
                sci_dir - science directory def in the environment variables
                combine_file - output file from frameCheck
                seeing_lower - lower limit for the seeing
                seeing_upper - upper limit for the seeing
                performance_limit - sky performance pass limit

        Output:
                sci_combined cubes for all of the objects imaged
                in the frames, provided at least a single frame
                passes the tests
        """

        # First step is to get a unique list
        # of the object names from the combine_file table
        # So that these can be looped over
        # to generate the .sof file in each case

        Table = np.loadtxt(combine_file, dtype='str')

        # The names are the third column

        ifu_names = Table[:, 3]

        # Create a unique list by taking a set

        ifu_names = list(set(ifu_names))

        # For each name execute the three helper reduce methods

        for name in ifu_names:

            # Initialise the names of the files

            deltaFile = sci_dir + '/' + name + '_delta_file.txt'

            shiftFile = sci_dir + '/' + name + '_shift_file.txt'

            print '[INFO]: Combining Object: %s ' % name

            new_Table = self.reduce_list_seeing(combine_file,
                                                seeing_lower,
                                                seeing_upper)

            name_Table = self.reduce_list_name(new_Table,
                                               name)

            combine_Table = self.reduce_list_sky(name_Table,
                                                 performance_limit)

            # Create the shift list for this object

            self.compute_shifts(sci_dir, combine_Table, name)

            # Plot the drift results

            self.plotDrift(deltaFile, name)

            # This combine_Table contains as first column
            # the reconstructed names to combine for that object
            # Want to write these out to a combine.sof file
            # - checking to see whether it exists already

            combine_name = name + '_combine.sof'

            if os.path.isfile(combine_name):

                os.system('rm %s' % combine_name)

            # Conditional execution of the recipes
            # depending on how much frames
            # Easiest if more than a single object

            if len(combine_Table) > 1:

                with open(combine_name, 'a') as f:

                    for row in combine_Table:

                        f.write('%s\n' % row[0])

                # Now execute the combine recipe
                # for this name given the sof file has been created

                os.system('esorex --output-dir=%s' % sci_dir
                          + ' kmos_combine --name=%s' % name
                          + ' --cmethod=median'
                          + ' --cpos_rej=2.0'
                          + ' --cneg_rej=2.0'
                          + ' --method="user"'
                          + ' --filename=%s' % shiftFile
                          + ' --edge_nan=TRUE %s' % combine_name)

            # If there is only a single object in this
            # seeing bin, isolate the core part of the name
            # and execute the kmo_sci_red recipe after
            # appending the object name to the sci_reduc.sof file
            # Since this is being executed in the calibrations
            # directory the sci_reduc.sof file is already there

            elif len(combine_Table) == 1:

                print 'This is the combined object' \
                      + ' name and type: %s %s' % (combine_Table[0][0],
                                                   type(combine_Table[0][0]))

                objName = combine_Table[0][0][
                    combine_Table[0][0].find('sci_reconstructed_') + 18:]

                real_name = \
                    combine_Table[0][1][: len(combine_Table[0][1])
                                        - combine_Table[0][1][::-1].find('/')] \
                    + objName

                # The sci_reduc.sof file is in the directory
                # - write out to this
                # Create a copy of the sci_reduc.sof in a new temporary file

                if os.path.isfile('sci_reduc_temp.sof'):

                    os.system('rm sci_reduc_temp.sof')
                with open(sci_dir + '/sci_reduc.sof') as f:
                    with open('sci_reduc_temp.sof', 'w') as f1:
                        for line in f:
                                f1.write(line)

                # Append the current object and skyfile
                # names to the newly created .sof file

                with open('sci_reduc_temp.sof', 'a') as f:

                    f.write('\n%s\tSCIENCE' % real_name)

                    f.write('\n%s\tSCIENCE' % combine_Table[0][1])

                # Now just execute the esorex recipe for this new file

                os.system('esorex --output-dir=%s' % sci_dir
                          + ' kmos_sci_red --pix_scale=0.1 --oscan=FALSE'
                          + ' --name=%s --sky_tweak=TRUE' % name
                          + ' --edge_nan=TRUE'
                          + ' sci_reduc_temp.sof')

                os.system('rm sci_reduc_temp.sof')

            # Final case - if the list is empty, do nothing

            else:

                print 'Nothing to combine for object %s' % name


    def plotDrift(self,
                  deltaFile,
                  name):

        """
        Def:
        Simple plotting method to visualise the drift
        away from given shift values.
        Takes the output delta shift file
        from compute_shifts and plots seperate
        graphs for the x_y_shift from each object.

        Input:
                deltaFile - file produced by compute_shifts,
                            x in column 1, y in column 2

                name - name of the object being plotted

        Output:
                x plot and y plot for each object
        """

        # Read in the x and y columns from the input file

        Table = np.loadtxt(deltaFile, dtype='str')

        # The names are the third column

        x_drift = Table[:, 0]

        y_drift = Table[:, 1]

        fileCount = np.arange(1, len(x_drift) + 1, 1)

        # Plot the results

        # Now make the plots for both nights,
        # want the same x-axis for all three layers

        f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(18.0, 10.0))

        ax1.plot(fileCount, x_drift, color='b')
        ax1.set_title('X-Shift ' + name, fontsize=24)
        ax1.tick_params(axis='y', which='major', labelsize=15)
        ax1.grid(b=True, which='both', linestyle='--')

        nbins = len(ax1.get_xticklabels())

        ax2.plot(fileCount, y_drift, color='g')
        ax2.set_title('Y-Shift ' + name, fontsize=24)
        ax2.set_xlabel(r'Frame Number', fontsize=24)
        ax2.tick_params(axis='both', which='major', labelsize=15)
        ax2.grid(b=True, which='both', linestyle='--')
        ax2.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
        ax2.set_xlim(min(fileCount) + 0.1, max(fileCount) - 0.1)

        f.subplots_adjust(hspace=0.001)
        f.tight_layout()
        # plt.show()
        f.savefig(name + '_drift.png')

    def multiExtractSpec(self,
                         sci_dir,
                         frameNames,
                         tracked_list,
                         **kwargs):

        """
        Def:
        Executes the kmo_sci_red recipe with skytweak
        for each of the individual object sky pairs and then
        applies compareSky to the science products. Plots a graph
        for the identification of bad science frames. Also fits a gaussian
        function to the collapsed data for each IFU, for each object sky pair.

        Input:
                reconstructed object cube from the current set
                skyCube - any reconstructed sky cube
                frameNames - list of object/sky pairs with
                            the names in the first column and type in second
                tracked_list - list of the names of three stars tracked
                                by KMOS, 1 star per detector. If this isn't
                                available write the name of one star three
                                times

        Outpt:
                Plot of frame performance against ID
        """

        # Need the sci_reduc.sof file in the directory

        if ((not (os.path.isfile('sci_reduc.sof')))):

            raise ValueError("Missing reduction .sof file")

        # remove the temporary sof file if it exists
        # also remove the tracked star names

        if os.path.isfile('sci_reduc_temp.sof'):

            os.system('rm sci_reduc_temp.sof')

        if os.path.isfile('tracked.txt'):

            os.system('rm tracked.txt')

        if os.path.isfile('tracked_one.txt'):

            os.system('rm tracked_one.txt')

        if os.path.isfile('tracked_two.txt'):

            os.system('rm tracked_two.txt')

        if os.path.isfile('tracked_three.txt'):

            os.system('rm tracked_three.txt')

        # Read in the data from the frameNames

        data = np.genfromtxt(frameNames, dtype='str')

        # Save the names and types as lists
        names = data[0:, 0]

        types = data[0:, 1]

        # Set the sci_comb names defined in the cubeclass

        if types[1] == 'O':

            combNames = cubeOps(names[1]).combNames

            rec_combNames = cubeOps(names[1]).rec_combNames

        elif types[2] == 'O':

            combNames = cubeOps(names[2]).combNames

            rec_combNames = cubeOps(names[2]).rec_combNames

        elif types[3] == 'O':

            combNames = cubeOps(names[3]).combNames

            rec_combNames = cubeOps(names[3]).rec_combNames

        else:

            print 'Having difficulty setting sci_comb names'

        track_name_one = tracked_list[0]

        track_name_two = tracked_list[1]

        track_name_three = tracked_list[2]

        # Loop round the list of combNames until each tracked_name appears

        for entry in combNames:

            if entry.find(track_name_one) != -1:

                tracked_star_one = sci_dir + '/' + entry

        for entry in combNames:

            if entry.find(track_name_two) != -1:

                tracked_star_two = sci_dir + '/' + entry

        for entry in combNames:

            if entry.find(track_name_three) != -1:

                tracked_star_three = sci_dir + '/' + entry

        # write to file each of these tracked names

        with open('tracked_one.txt', 'a') as f:

            f.write(tracked_star_one)

        with open('tracked_two.txt', 'a') as f:

            f.write(tracked_star_two)

        with open('tracked_three.txt', 'a') as f:

            f.write(tracked_star_three)

        # execute the framecheck method to complete the science reduction

        ID, offList, namesVec, IFUValVec, \
            frameValVec, fwhmValVec_one, \
            fwhmValVec_two, fwhmValVec_three, \
            fwhmDict = self.frameCheck(sci_dir,
                                       frameNames,
                                       tracked_list)

        # There are the IFU sky tweak performance plots
        # The mean sky tweak performance across each IFU

        print '[INFO]: Plotting frame performance against IFU number'

        self.meanIFUPlot(offList, namesVec, IFUValVec)

        # The performance of skytweak in each IFU

        print '[INFO]: Plotting Individual IFU performance with frame'

        self.indIFUPlot(offList, ID, IFUValVec)

        # Mean skytweak as a function of frame

        print '[INFO]: Plotting mean sky subtraction performance'

        self.meanFramePlot(ID, frameValVec)

        # FWHM PLOT - monitoring seeing across the frames

        print '[INFO]: Plotting evolution of tracked star FWHM'

        self.meanFWHMPlot(ID, fwhmValVec_one, 'star_one')
        self.meanFWHMPlot(ID, fwhmValVec_two, 'star_two')
        self.meanFWHMPlot(ID, fwhmValVec_three, 'star_three')

        # If supplying an additional list of files, construct the double plots
        if kwargs:

            print '[INFO]: Additional keyword arguments supplied' \
                  + ' - Checking additional frames and double plotting'

            additional_frameNames = kwargs.values()[0]

            ID1, offList1, namesVec1, IFUValVec1, frameValVec1, \
                fwhmValVec1_one, fwhmValVec1_two, fwhmValVec1_three, \
                fwhmDict1 = self.frameCheck(sci_dir,
                                            additional_frameNames,
                                            tracked_list)

            # Now have two sets of all the parameters
            # and can make the double plots using the multiplot methods

            print '[INFO]: Plotting multi mean IFU performance'

            self.multiMeanFramePlot(ID, frameValVec, frameValVec1)

            print '[INFO]: Plotting multi IFU subplots'

            self.multiIndIFUPlot(offList, ID, IFUValVec, IFUValVec1)

            print '[INFO]: Plotting multi FWHM plot'

            self.multiMeanFWHMPlot(ID,
                                   fwhmValVec_one,
                                   fwhmValVec1_one,
                                   'star_one')

            self.multiMeanFWHMPlot(ID,
                                   fwhmValVec_two,
                                   fwhmValVec1_two,
                                   'star_two')

            self.multiMeanFWHMPlot(ID,
                                   fwhmValVec_three,
                                   fwhmValVec1_three,
                                   'star_three')

        # list the objects in the different seeing bins

        print '[INFO]: These are the best names: %s ' % fwhmDict['Best']

        print '[INFO]: These are the Good names: %s ' % fwhmDict['Good']

        print '[INFO]: These are the Okay names: %s ' % fwhmDict['Okay']

        print '[INFO]: These are the Bad names: %s ' % fwhmDict['Bad']

    def saveSpec(self,
                 sci_dir,
                 cubeName):

        """
        Def:
        Extract a spectrum from a given cube optimally
        and save the spectrum to a fits file in the same
        format as the frameCheck method
        Input - cube: Any reconstructed cube, object or sky
        """

        # Create an instance of the cube class
        cube = cubeOps(cubeName)

        # extract the properties
        wave_arr = cube.wave_array

        spec = cube.centralSpec()

        # Save to fits file
        tbhdu = fits.new_table(fits.ColDefs(
            [fits.Column(name='Wavelength',
                         format='E',
                         array=wave_arr),
             fits.Column(name='Flux',
                         format='E',
                         array=spec)]))

        prihdu = fits.PrimaryHDU(header=cube.primHeader)

        thdulist = fits.HDUList([prihdu, tbhdu])

        # Find the name of the cube

        if cubeName.find("/") == -1:

            sky_name = sci_dir + '/' + cubeName

        # Otherwise the directory structure is included and have to
        # search for the backslash and omit up to the last one

        else:

            objName = cubeName[len(cubeName) - cubeName[::-1].find("/"):]
            sky_name = sci_dir + '/' + objName

        thdulist.writeto(sky_name[:-5] + '_spectrum.fits', clobber=True)

    def plotSpecs(self,
                  sci_dir,
                  objSpec,
                  skyCube,
                  n):

        """
        Def:
        Takes the object and sky spectra, bins according to n
        which must be a factor of the wavelength array and plots both
        on the same axes
        Input
                - objSpec: Input spectrum, must be in the fits format specified
                            in the frameCheck recipe
                            i.e. Table data with one column that has header
                            FLUX and one with header WAVELENGTH
              - skySpec: Input sky spectrum in the same format
              - n: Binning order, must be an integer
        """

        if type(n) != int:

            raise TypeError('n must be an integer')

        if n <= 0:

            raise ValueError('You have specified a value'
                             + ' less than or equal 0 for n')

        # Read in the files

        objTable = fits.open(objSpec)

        obj_spec = objTable[1].data['Flux']

        obj_wave = objTable[1].data['Wavelength']

        # skyCube

        self.saveSpec(sci_dir, skyCube)

        # Find the name of the skyCube in the science directory

        if skyCube.find("/") == -1:

            sky_name = sci_dir + '/' + skyCube

        # Otherwise the directory structure is included and have to
        # search for the backslash and omit up to the last one

        else:

            objName = skyCube[len(skyCube) - skyCube[::-1].find("/"):]
            sky_name = sci_dir + '/' + objName

        skyTable = fits.open(sky_name[:-5] + '_spectrum.fits')

        sky_spec = skyTable[1].data['Flux']

        sky_wave = skyTable[1].data['Wavelength']

        if n == 1:

            new_obj_spec = copy(obj_spec)
            new_obj_wave = copy(obj_wave)
            new_sky_spec = copy(sky_spec)
            new_sky_wave = copy(sky_wave)

        elif n > 1:

            # Variables to house the new binned spectra

            new_obj_spec = []
            new_obj_wave = []
            new_sky_spec = []
            new_sky_wave = []

            # Counters for the binning
            lower = 0
            upper = copy(n)

            # Bin the data

            for i in range(len(obj_wave) / n):

                # The binned spectra are the sum over the ranges

                new_obj_spec.append(sum(obj_spec[lower:upper]))

                new_sky_spec.append(sum(sky_spec[lower:upper]))

                # The binned wavelengths are the median over the ranges

                new_obj_wave.append(np.median(obj_wave[lower:upper]))

                new_sky_wave.append(np.median(sky_wave[lower:upper]))

                lower += n
                upper += n

            # Print to make sure they are the same length

            print len(new_obj_spec), len(new_sky_spec)

            new_obj_spec = np.array(new_obj_spec)
            new_obj_wave = np.array(new_obj_wave)
            new_sky_spec = np.array(new_sky_spec)
            new_sky_wave = np.array(new_sky_wave)

        # Plot the results
        # Now make the plots for both nights,
        # want the same x-axis for all three layers

        f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(18.0, 10.0))

        ax1.plot(new_obj_wave, new_obj_spec, color='b')
        ax1.set_title('Object and Sky Comparison', fontsize=24)
        ax1.set_ylim(-0.25E-17, 3E-17)
        ax1.tick_params(axis='y', which='major', labelsize=15)

        nbins = len(ax1.get_xticklabels())

        ax2.plot(new_sky_wave, new_sky_spec, color='g')
        ax2.set_xlabel(r'Wavelength ($\mu m$)', fontsize=24)
        ax2.tick_params(axis='both', which='major', labelsize=15)
        ax2.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
        ax2.set_xlim(1.95, 2.35)

        f.subplots_adjust(hspace=0.001)
        f.tight_layout()

        # plt.show()
        f.savefig(objSpec[:-5] + '.png')

    def telluric_correct(self,
                         grating,
                         cal_dir):

        """
        Def: Currently the standard pipeline telluric correction isn't
        working properly - it doesn't do a good job for the absorption
        features. This method fits a polynomial to the three telluric
        spectra stored in IFUs 3, 12 and 18 and replaces the more detailed
        telluric correction with this fit.

        Input:
                grating - the grating ID in the fits header

        Output:
                telluric_XXX - where X is the grating value
        """

        # start by checking the value of the grating
        if grating == 'K':

            table = fits.open(cal_dir + '/telluric_KKK.fits', mode='update')

        elif grating == 'HK':

            table = fits.open(cal_dir + '/telluric_HKHKHK.fits', mode='update')

        elif grating == 'H':

            table = fits.open(cal_dir + '/telluric_HHH.fits', mode='update')

        elif grating == 'YJ':

            table = fits.open(cal_dir + '/telluric_YJYJYJ.fits', mode='update')

        elif grating == 'IZ':

            table = fits.open(cal_dir + '/telluric_IZIZIZ.fits', mode='update')

        else:

            raise ValueError('The supplied grating ID'
                             + ' is not recognised')

        # construct the wavelength array

        header = table[5].header

        start_l = header['CRVAL1']

        delta = header['CDELT1']

        wave_array = start_l + np.arange(0, 2048 * delta, delta)

        # assign the data for IFUs 3, 12, 18

        data_1 = table[5].data

        data_2 = table[23].data

        data_3 = table[35].data

        # create masked versions of these for the polynomial fit

        data_1_masked = ma.masked_where(np.isnan(data_1), data_1, copy=True)

        data_2_masked = ma.masked_where(np.isnan(data_2), data_2, copy=True)

        data_3_masked = ma.masked_where(np.isnan(data_3), data_3, copy=True)

        # construct polynomial model from lmfit

        mod = PolynomialModel(6)

        pars = mod.guess(data_1_masked, x=wave_array)

        # for the masked array to work need to assign the parameters

        pars['c0'].set(value=0)
        pars['c1'].set(value=0)
        pars['c2'].set(value=0)
        pars['c3'].set(value=0)
        pars['c4'].set(value=0)
        pars['c5'].set(value=0)
        pars['c6'].set(value=0)
        # pars['c7'].set(value=0)

        # these will be the initial parameters for all the fits

        out_1 = mod.fit(data_1_masked, pars, x=wave_array)
        out_1_spec = out_1.best_fit

        out_2 = mod.fit(data_2_masked, pars, x=wave_array)
        out_2_spec = out_2.best_fit

        out_3 = mod.fit(data_3_masked, pars, x=wave_array)
        out_3_spec = out_3.best_fit

        # quickly plot these

        fig, ax = plt.subplots(3, 1, figsize=(18, 18))

        ax[0].plot(wave_array, data_1_masked)
        ax[0].plot(wave_array, out_1_spec)

        ax[1].plot(wave_array, data_2_masked)
        ax[1].plot(wave_array, out_2_spec)

        ax[2].plot(wave_array, data_3_masked)
        ax[2].plot(wave_array, out_3_spec)

        # plt.show()

        # now re-insert these data into the telluric file
        # and overwrite the values that are there currently

        table[5].data = out_1_spec

        table[23].data = out_2_spec

        table[35].data = out_3_spec

        table.flush()

        table.close()

        # original file has now been updated and closed



    def arrayCompare(self,
                     profile,
                     imData):

        """
        Def: Function to compare the shapes of the optimal extraction profile
        and the object image data which may not necessarily be the same. But in
        order for the optimal extraction from
        profile to work they have to be. This
        method spits out an adjusted profile so that
        the shape matches that of the
        object data, allowing them to be multiplied together.

        Input:
                profile - the optimal extraction profile
                imData - the object cube data

        Output:
                new_profile - adjusted profile with
                              shape matches the last two indices
                              of the object data
        """

        if imData.shape[1] != profile.shape[0]:

            print 'Must adjust rows'

            # decide whether to add rows or clip rows

            if imData.shape[1] - profile.shape[0] > 0:

                # need to add rows to the profile

                print 'Adding Rows'

                row_addition = np.zeros([abs(imData.shape[1]
                                        - profile.shape[0]), profile.shape[1]])

                # vstack this together with the original profile

                profile = np.vstack((profile, row_addition))

            else:

                # need to clip rows from the profile,
                # assuming we're at the edge so it won't affect the process

                print 'Deleting Rows'

                profile = np.delete(profile,
                                    np.arange(imData.shape[1],
                                              profile.shape[0],
                                              1),
                                    axis=0)

            # Now deal with the columns

        if imData.shape[2] != profile.shape[1]:

            print 'Must adjust columns'

            # Decide whether to add or clip columns

            if imData.shape[2] - profile.shape[1] > 0:

                # need to add columns

                print 'Adding Columns'

                column_addition = np.zeros([profile.shape[0],
                                           abs(imData.shape[2] -
                                           profile.shape[1])])

                # hstack this togehter with the original profile

                profile = np.hstack((profile, column_addition))

            else:

                # need to clip the columns

                print 'Deleting columns'

                profile = np.delete(profile,
                                    np.arange(imData.shape[2],
                                              profile.shape[1],
                                              1),
                                    axis=1)

        # return the profile

        return profile

    # Starting to write some functions for extracting
    # spectra from galaxies and measuring their properties
    def galExtract(self,
                   sci_dir,
                   std_cube,
                   obj_cube,
                   sky_cube,
                   center_x,
                   center_y,
                   n):

        """
        Def: Function for extracting the optimal spectrum from a galaxy at the
        location specified by the user, after
        examining the object cube in qfits.
        Recovers the object spectrum and makes a plot of this in the usual way
        with the sky spectrum beneath.

        Input:
                std_cube - data cube for the IFU dedicated to observing a star
                obj_cube - the cube containing the
                             galaxy to extract a spectrum from
                sky_cube - one of the sky cubes extracted
                             in the runPipeline.py analysis
                             note that this is only used in the plotting part.
                center_x - x coordinate of the object
                           center determined from qfits
                center_y - y coordinate of the object
                            center determined from qfits
                n - binning for spectrum plot

        Output: Plot of the object and sky spectra using plotSpecs
                If required, a data file containing the optimally
                extracted object spectrum
                returns - 1D optimally extracted spectrum
        """

        # First extract the profile from the standard star cube

        std_star_cube = cubeOps(std_cube)

        gal_obj_cube = cubeOps(obj_cube)

        wl = gal_obj_cube.wave_array

        # Find the cube name

        if obj_cube.find("/") == -1:

            gal_name = sci_dir + '/' + obj_cube

        # Otherwise the directory structure is included and have to
        # search for the backslash and omit up to the last one

        else:

            objName = obj_cube[len(obj_cube) - obj_cube[::-1].find("/"):]

            gal_name = sci_dir + '/' + objName

        # Extract the PSF profile using the method in cubeClass

        params, std_profile, FWHM, offList = std_star_cube.psfMask()

        # Now roll this profile over to the central location of the object
        # Use numpy.roll to shift the psfMask to the location of the object
        # Find the shift values

        x_shift = int(center_x - params['center_x'])

        y_shift = int(center_y - params['center_y'])

        new_mask = np.roll(std_profile, x_shift, axis=0)

        # For the x_shift need to loop round the elements of the new_mask

        final_new_mask = []

        for arr in new_mask:

            final_new_mask.append(np.roll(arr, y_shift))

        final_new_mask = np.array(final_new_mask)

        # Check the shapes of the profile and object data, adjust accordingly

        final_new_mask = self.arrayCompare(final_new_mask, gal_obj_cube.data)

        # Check to see that the gaussian and shifted profile align

        colFig, colAx = plt.subplots(1, 1, figsize=(14.0, 14.0))

        colCax = colAx.imshow(final_new_mask, interpolation='bicubic')

        colFig.colorbar(colCax)

        # plt.show()

        plt.close('all')

        # Extract each optimal spectrum

        optimal_spec = gal_obj_cube.optimalSpecFromProfile(final_new_mask,
                                                           FWHM,
                                                           center_y,
                                                           center_x)
        # Save the optimal spectrum for each object
        # Need to create a new fits table for this

        tbhdu = fits.new_table(fits.ColDefs(
            [fits.Column(name='Wavelength',
                         format='E',
                         array=gal_obj_cube.wave_array),
             fits.Column(name='Flux',
                         format='E',
                         array=optimal_spec)]))

        prihdu = fits.PrimaryHDU(header=gal_obj_cube.primHeader)

        thdulist = fits.HDUList([prihdu, tbhdu])

        thdulist.writeto(gal_name[:-5] + '_spectrum.fits', clobber=True)

        # Create a plot of both the sky and the object next to one another
        self.plotSpecs(sci_dir, gal_name[:-5] + '_spectrum.fits', sky_cube, n)

        return optimal_spec, wl

    def multiGalExtract(self,
                        inFile,
                        res):

        """
        Def: Function to perform galExtract for a selection of object supplied
        within an infile.

        Input:
                inFile - with top line containing the gal_dir, std_cube_gal and
                         sky_gal. All other lines contain the
                         object cube and x y central positions
                res - resolution with which to extract the spectra
        """

        # first read in the top line containing the standard quantities

        Table = np.genfromtxt(inFile, dtype='str')

        gal_dir = Table[0, :][0]

        std_gal = Table[0, :][1]

        sky_gal = Table[0, :][2]

        # Now assign the object names and x, y coordinates

        obj_names = Table[1:, 0]

        x_cen = Table[1:, 3]

        y_cen = Table[1:, 2]

        # Loop around these objects and use galExtract

        for item in zip(obj_names, x_cen, y_cen):

            optimal_spec, wl = self.galExtract(gal_dir,
                                               std_gal,
                                               item[0],
                                               sky_gal,
                                               int(item[1]),
                                               int(item[2]),
                                               res)

    def singlePixelExtract_OIII5008(self,
                                    sci_dir,
                                    obj_cube,
                                    centre_x,
                                    centre_y,
                                    z,
                                    n):

        """
        Def: Function for extracting the spectrum
        from a stacked galaxy image pixel by pixel at
        locations specified by the user, after examining
        the object cube in qfits.
        Recovers the object spectrum and makes a plot of this in the usual way
        with the sky spectrum beneath.

        Input:
                obj_cube - galaxy cube to extract a spectrum from
                sky_cube - sky cubes extracted in the runPipeline.py analysis
                centre_x - galaxy central pixel location in x-direction
                centre_y - galaxy central pixel location in y-direction
                n - binning for spectrum plot

        Output:
                Plot of the object and sky spectra using plotSpecs
                If required, a data file containing the
                optimally extracted object spectrum
                Data file containing the optimally extracted sky spectrum
        """

        z = float(z)

        centre_x = int(centre_x)

        centre_y = int(centre_y)

        # Find the cube name

        if obj_cube.find("/") == -1:

            gal_name = sci_dir + '/' + obj_cube

        # Otherwise the directory structure is included and have to
        # search for the backslash and omit up to the last one

        else:

            objName = obj_cube[len(obj_cube) - obj_cube[::-1].find("/"):]

            gal_name = sci_dir + '/' + objName

        # First find the noise by extracting a
        # spectrum from the skycube
        # Define the lower and upper wavelength
        # ranges for the gaussian fit to the OIII line

        OIII5008_lower = (0.5008 * (1 + z)) - 0.01

        OIII5008_upper = (0.5008 * (1 + z)) + 0.01

        # Now compute the noise value. This will be done by
        # extracting a spectrum from a single spaxel
        # where there is no object flux, restricting
        # the pixel range to the region around the OIII line
        # plotting a histogram of the flux values, fitting
        # with a gaussian and taking the sigma of the
        # gaussian as the noise level

        objCube = cubeOps(obj_cube)

        objWavelength = objCube.wave_array

        fluxArray = objCube.singlePixelExtract(5, 5)

        indices = \
            np.where(np.logical_and(objWavelength > OIII5008_lower - 0.05,
                                    objWavelength < OIII5008_upper + 0.05))[0]

        fit_flux = fluxArray[indices]

        # Have now the flux array, the bins
        # will span from -2 to 2 in 0.2 intervals

        bins = np.arange(-2, 2, 0.2)

        hist, edges = np.histogram(fit_flux, bins=bins)

        # Now fit a gaussian to this data

        gmod = GaussianModel()

        # Guess and set the initial parameter values

        pars = gmod.make_params()
        pars['center'].set(0)
        pars['sigma'].set(0.5)
        pars['amplitude'].set(1000)

        # perform the fit

        out = gmod.fit(hist, pars, x=bins[0:-1])

        # Quickly plot the results to make sure all good

        fig, ax = plt.subplots(1, figsize=(18, 10))

        ax.plot(bins[0:-1], hist, color='black')
        ax.plot(bins[0:-1], out.best_fit, color='green')

        plt.show()

        # Use the width of the fitted gaussian to find the noise value

        noise_value = out.best_values['sigma']

        # Now have our estimate of the noise level of the measurements
        # Extract a spectrum using the cubeOps class

        central_spec = objCube.singlePixelExtract(centre_x, centre_y)

        # Now split the flux array into the
        # region around the OIII emission line
        # Will construct x and y arrays housing the S/N value.
        # Need to loop around both of those
        # and plug in the pixel values to extract the spectrum

        x_loop_array = np.arange(1, len(objCube.data[0][0]) - 1)

        y_loop_array = np.arange(1, len(objCube.data[0]) - 1)

        x_SN_array = []

        y_SN_array = []

        print 'Checking the S/N of object: %s' % (objName)

        # for each of the values in the x_loop_array
        # and y_loop_array compute the signal to noise of OIII

        for item in x_loop_array:

            flux = objCube.singlePixelExtract(item, centre_y)

            indices = \
                np.where(np.logical_and(objWavelength > OIII5008_lower,
                                        objWavelength < OIII5008_upper))[0]

            fit_wavelength = objWavelength[indices]

            fit_flux = flux[indices]

            # Now have both the fit flux and wavelength
            # - fit a gaussian to measure the amplitude of emission line

            best_values, covar = self.fitSingleGauss(fit_wavelength,
                                                     fit_flux,
                                                     0.5008 * (1 + z))
            # print covar

            # define some maths helper variables

            amp = abs(best_values['amplitude'])

            num = 100 * np.sqrt(covar[1][1])

            if covar is None:

                x_SN_array.append(0.0)

            elif num / amp > 30 or best_values['sigma'] > 0.007:

                x_SN_array.append(0.0)

            else:

                # Note the 0.00028 is the distance between
                # adjacent wavelength pixels for normalisation

                x_SN_array.append(abs(best_values['amplitude']) /
                                  (noise_value * 0.00028))

        # Now for the y_array

        print 'Moving onto y array'

        for item in y_loop_array:

            flux = objCube.singlePixelExtract(item, centre_x)

            indices = \
                np.where(np.logical_and(objWavelength > OIII5008_lower,
                                        objWavelength < OIII5008_upper))[0]

            fit_wavelength = objWavelength[indices]

            fit_flux = flux[indices]

            # Now have both the fit flux and wavelength
            # - fit a gaussian to measure the amplitude of emission line

            best_values, covar = self.fitSingleGauss(fit_wavelength,
                                                     fit_flux,
                                                     0.5008 * (1 + z))
            # print covar

            # define some maths helper variables

            amp = abs(best_values['amplitude'])

            num = 100 * np.sqrt(covar[1][1])

            if covar is None:

                y_SN_array.append(0.0)

            elif num / amp > 30 or best_values['sigma'] > 0.007:

                y_SN_array.append(0.0)

            else:

                # Note the 0.00028 is the distance between
                # adjacent wavelength pixels for normalisation

                y_SN_array.append(abs(best_values['amplitude']) /
                                  (noise_value * 0.00028))

        f, (ax1, ax2) = plt.subplots(2, 1, sharex=False, figsize=(18.0, 10.0))

        ax1.plot(x_loop_array, x_SN_array, color='b')
        ax1.set_title('%s x-direction S/N OIII5008' % objName[26:-5],
                      fontsize=24)
        ax1.set_xlabel(r'x-pixel Position at y = %s' % centre_y, fontsize=24)
        ax1.tick_params(axis='y', which='major', labelsize=15)

        nbins = len(ax1.get_xticklabels())

        ax2.plot(y_loop_array, y_SN_array, color='g')
        ax2.set_xlabel(r'y-pixel Position at x = %s' % centre_x, fontsize=24)
        ax2.tick_params(axis='both', which='major', labelsize=15)
        ax2.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))

        f.subplots_adjust(hspace=0.001)
        f.tight_layout()

        plt.show()

        f.savefig(gal_name[:-5] + '_OIII5008SN.png')

        # Take the median of this array where the flux values don't exceed 50
        # plot an initial evaluation of the spectrum

        f, ax1 = plt.subplots(1, 1, sharex=True, figsize=(18.0, 10.0))

        ax1.plot(objWavelength, fluxArray, color='b')
        ax1.plot(objWavelength, central_spec, color='red')
        ax1.set_title('Sky Spectrum', fontsize=30)
        ax1.set_ylim(-1, 10)
        ax1.tick_params(axis='y', which='major', labelsize=15)
        ax1.set_xlabel(r'Wavelength ($\mu m$)', fontsize=24)
        ax1.set_ylabel(r'Flux', fontsize=24)

        f.tight_layout()

        plt.show()

    def singlePixelExtractMulti_OIII5008(self,
                                         inFile,
                                         sci_dir):

        """
        Def:
        Applies the singlePixelExtract function to a list of files,
        in the infile we have the name of the object to read in,
        the redshift and the x,y positions
        of the center of the galaxy from Qfits
        Input:
                inFile - Containing the object attributes listed above
        """

        # Want to apply the above function to everything within a file
        # First read in the whole file as a Table

        Table = np.loadtxt(inFile, dtype='str')

        for row in Table:
            self.singlePixelExtract_OIII5008(sci_dir,
                                             row[0],
                                             row[2],
                                             row[3],
                                             row[1],
                                             1)

    def singlePixelExtract_OIII4960(self,
                                    sci_dir,
                                    obj_cube,
                                    centre_x,
                                    centre_y,
                                    z,
                                    n):
        """
        Def: Function for extracting the spectrum from a
        stacked galaxy image pixel by pixel at
        locations specified by the user,
        after examining the object cube in qfits.
        Recovers the object spectrum and makes a plot of this in the usual way
        with the sky spectrum beneath.

        Input:
                obj_cube - the cube containing the galaxy
                            to extract a spectrum from
                sky_cube - one of the sky cubes extracted
                            in the runPipeline.py analysis
                centre_x - galaxy central pixel location in x-direction
                centre_y - galaxy central pixel location in y-direction
                n - binning for spectrum plot

        Output:
                Plot of the object and sky spectra using plotSpecs
                If required, a data file containing the
                optimally extracted object spectrum.
                Data file containing the optimally extracted sky spectrum

        """

        z = float(z)

        centre_x = int(centre_x)

        centre_y = int(centre_y)

        # Find the cube name

        if obj_cube.find("/") == -1:

            gal_name = sci_dir + '/' + obj_cube

        # Otherwise the directory structure is included and have to
        # search for the backslash and omit up to the last one

        else:

            objName = obj_cube[len(obj_cube) - obj_cube[::-1].find("/"):]

            gal_name = sci_dir + '/' + objName

        # First find the noise by extracting a spectrum from the skycube
        # Define the lower and upper wavelength ranges
        # for the gaussian fit to the OIII line

        OIII4960_lower = (0.4960 * (1 + z)) - 0.01

        OIII4960_upper = (0.4960 * (1 + z)) + 0.01

        # Now compute the noise value. This will be done by
        # extracting a spectrum from a single spaxel
        # where there is no object flux, restricting
        # the pixel range to the region around the OIII line
        # plotting a histogram of the flux values, fitting
        # with a gaussian and taking the sigma of the
        # gaussian as the noise level

        objCube = cubeOps(obj_cube)

        objWavelength = objCube.wave_array

        fluxArray = objCube.singlePixelExtract(5, 5)

        indices = \
            np.where(np.logical_and(objWavelength > OIII4960_lower - 0.05,
                                    objWavelength < OIII4960_upper + 0.05))[0]

        fit_flux = fluxArray[indices]

        # Have now the flux array, the bins will
        # span from -2 to 2 in 0.2 intervals

        bins = np.arange(-2, 2, 0.2)

        hist, edges = np.histogram(fit_flux, bins=bins)

        # Now fit a gaussian to this data

        gmod = GaussianModel()

        # Guess and set the initial parameter values

        pars = gmod.make_params()
        pars['center'].set(0)
        pars['sigma'].set(0.5)
        pars['amplitude'].set(1000)

        # perform the fit

        out = gmod.fit(hist, pars, x=bins[0:-1])

        # Quickly plot the results to make sure all good

        fig, ax = plt.subplots(1, figsize=(18, 10))

        ax.plot(bins[0:-1], hist, color='black')

        ax.plot(bins[0:-1], out.best_fit, color='green')

        plt.show()

        # Use the width of the fitted gaussian to find the noise value

        noise_value = out.best_values['sigma']

        # Now have our estimate of the noise level of the measurements
        # Extract a spectrum using the cubeOps class

        central_spec = objCube.singlePixelExtract(centre_x, centre_y)

        # Now split the flux array into the
        # region around the OIII emission line
        # Will construct x and y arrays housing
        # the S/N value. Need to loop around
        # both of those and plug in the pixel values to extract the spectrum

        x_loop_array = np.arange(1, len(objCube.data[0][0]) - 1)

        y_loop_array = np.arange(1, len(objCube.data[0]) - 1)

        x_SN_array = []

        y_SN_array = []

        print 'Checking the S/N of object: %s' % (objName)

        # for each of the values in the x_loop_array
        # and y_loop_array compute the signal to noise of OIII

        for item in x_loop_array:

            flux = objCube.singlePixelExtract(item, centre_y)

            indices = \
                np.where(np.logical_and(objWavelength > OIII4960_lower,
                         objWavelength < OIII4960_upper))[0]

            fit_wavelength = objWavelength[indices]

            fit_flux = flux[indices]

            # Now have both the fit flux and wavelength -
            # fit a gaussian to measure the amplitude of emission line

            best_values, covar = self.fitSingleGauss(fit_wavelength,
                                                     fit_flux,
                                                     0.4960 * (1 + z))

            # define some maths helper variables

            amp = abs(best_values['amplitude'])

            num = 100 * np.sqrt(covar[1][1])

            if covar is None:

                x_SN_array.append(0.0)

            elif num / amp > 30 or best_values['sigma'] > 0.007:

                x_SN_array.append(0.0)

            else:

                # Note the 0.00028 is the distance between
                # adjacent wavelength pixels for normalisation

                x_SN_array.append(amp / (noise_value * 0.00028))

        # Now for the y_array

        print 'Moving onto y array'

        for item in y_loop_array:

            flux = objCube.singlePixelExtract(item, centre_x)

            indices = \
                np.where(np.logical_and(objWavelength > OIII4960_lower,
                         objWavelength < OIII4960_upper))[0]

            fit_wavelength = objWavelength[indices]

            fit_flux = flux[indices]

            # Now have both the fit flux and wavelength -
            # fit a gaussian to measure the amplitude of emission line

            best_values, covar = self.fitSingleGauss(fit_wavelength,
                                                     fit_flux,
                                                     0.4960 * (1 + z))

            # define some maths helper variables

            amp = abs(best_values['amplitude'])

            num = 100 * np.sqrt(covar[1][1])

            if covar is None:

                y_SN_array.append(0.0)

            elif num / amp > 30 or best_values['sigma'] > 0.007:

                y_SN_array.append(0.0)

            else:

                # Note the 0.00028 is the distance between
                # adjacent wavelength pixels for normalisation

                y_SN_array.append(amp / (noise_value * 0.00028))

        f, (ax1, ax2) = plt.subplots(2, 1, sharex=False, figsize=(18.0, 10.0))

        ax1.plot(x_loop_array, x_SN_array, color='b')
        ax1.set_title('%s x-direction S/N OIII4960' % objName[26:-5],
                      fontsize=24)
        ax1.set_xlabel(r'x-pixel Position at y = %s' % centre_y, fontsize=24)
        ax1.tick_params(axis='y', which='major', labelsize=15)

        nbins = len(ax1.get_xticklabels())

        ax2.plot(y_loop_array, y_SN_array, color='g')
        ax2.set_xlabel(r'y-pixel Position at x = %s' % centre_x, fontsize=24)
        ax2.tick_params(axis='both', which='major', labelsize=15)
        ax2.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))

        f.subplots_adjust(hspace=0.001)
        f.tight_layout()

        plt.show()

        f.savefig(gal_name[:-5] + '_OIII4960SN.png')

        # Take the median of this array where the flux values don't exceed 50
        # plot an initial evaluation of the spectrum

        f, ax1 = plt.subplots(1, 1, sharex=True, figsize=(18.0, 10.0))

        ax1.plot(objWavelength, fluxArray, color='b')
        ax1.plot(objWavelength, central_spec, color='red')
        ax1.set_title('Sky Spectrum', fontsize=30)
        ax1.set_ylim(-1, 10)
        ax1.tick_params(axis='y', which='major', labelsize=15)
        ax1.set_xlabel(r'Wavelength ($\mu m$)', fontsize=24)
        ax1.set_ylabel(r'Flux', fontsize=24)

        f.tight_layout()

        plt.show()

    def singlePixelExtractMulti_OIII4960(self,
                                         inFile,
                                         sci_dir):

        """
        Def:
        Applies the singlePixelExtract function to a list of files,
        in the infile we have the name of the object to read in,
        the redshift and the x,y positions
        of the center of the galaxy from Qfits.

        Input:
                inFile - Containing the object attributes listed above
        """
        # Want to apply the above function to everything within a file
        # First read in the whole file as a Table

        Table = np.loadtxt(inFile, dtype='str')

        for row in Table:

            self.singlePixelExtract_OIII4960(sci_dir,
                                             row[0],
                                             row[2],
                                             row[3],
                                             row[1],
                                             1)

    def singlePixelExtract_Hb(self,
                              sci_dir,
                              obj_cube,
                              centre_x,
                              centre_y,
                              z,
                              n):

        """
        Def: Function for extracting the spectrum from a
        stacked galaxy image pixel by pixel at
        locations specified by the user,
        after examining the object cube in qfits.
        Recovers the object spectrum and makes a plot of this in the usual way
        with the sky spectrum beneath.

        Input:
                obj_cube - the cube containing the galaxy
                            to extract a spectrum from
                sky_cube - one of the sky cubes extracted
                            in the runPipeline.py analysis
                centre_x - galaxy central pixel location in x-direction
                centre_y - galaxy central pixel location in y-direction
                n - binning for spectrum plot

        Output:
                Plot of the object and sky spectra using plotSpecs
                If required, a data file containing the
                optimally extracted object spectrum.
                Data file containing the optimally extracted sky spectrum

        """

        z = float(z)

        centre_x = int(centre_x)

        centre_y = int(centre_y)

        # Find the cube name

        if obj_cube.find("/") == -1:

            gal_name = sci_dir + '/' + obj_cube

        # Otherwise the directory structure is included and have to
        # search for the backslash and omit up to the last one

        else:

            objName = obj_cube[len(obj_cube) - obj_cube[::-1].find("/"):]

            gal_name = sci_dir + '/' + objName

        # First find the noise by extracting a spectrum from the skycube
        # Define the lower and upper wavelength ranges
        # for the gaussian fit to the OIII line

        Hb_lower = (0.4861 * (1 + z)) - 0.01

        Hb_upper = (0.4861 * (1 + z)) + 0.01

        # Now compute the noise value. This will be done by
        # extracting a spectrum from a single spaxel
        # where there is no object flux, restricting
        # the pixel range to the region around the OIII line
        # plotting a histogram of the flux values, fitting
        # with a gaussian and taking the sigma of the
        # gaussian as the noise level

        objCube = cubeOps(obj_cube)

        objWavelength = objCube.wave_array

        fluxArray = objCube.singlePixelExtract(5, 5)

        indices = \
            np.where(np.logical_and(objWavelength > Hb_lower - 0.05,
                                    objWavelength < Hb_upper + 0.05))[0]

        fit_flux = fluxArray[indices]

        # Have now the flux array, the bins will
        # span from -2 to 2 in 0.2 intervals

        bins = np.arange(-2, 2, 0.2)

        hist, edges = np.histogram(fit_flux, bins=bins)

        # Now fit a gaussian to this data

        gmod = GaussianModel()

        # Guess and set the initial parameter values

        pars = gmod.make_params()
        pars['center'].set(0)
        pars['sigma'].set(0.5)
        pars['amplitude'].set(1000)

        # perform the fit

        out = gmod.fit(hist, pars, x=bins[0:-1])

        # Quickly plot the results to make sure all good

        fig, ax = plt.subplots(1, figsize=(18, 10))

        ax.plot(bins[0:-1], hist, color='black')

        ax.plot(bins[0:-1], out.best_fit, color='green')

        plt.show()

        # Use the width of the fitted gaussian to find the noise value

        noise_value = out.best_values['sigma']

        # Now have our estimate of the noise level of the measurements
        # Extract a spectrum using the cubeOps class

        central_spec = objCube.singlePixelExtract(centre_x, centre_y)

        # Now split the flux array into the
        # region around the OIII emission line
        # Will construct x and y arrays housing
        # the S/N value. Need to loop around
        # both of those and plug in the pixel values to extract the spectrum

        x_loop_array = np.arange(1, len(objCube.data[0][0]) - 1)

        y_loop_array = np.arange(1, len(objCube.data[0]) - 1)

        x_SN_array = []

        y_SN_array = []

        print 'Checking the S/N of object: %s' % (objName)

        # for each of the values in the x_loop_array
        # and y_loop_array compute the signal to noise of OIII

        for item in x_loop_array:

            flux = objCube.singlePixelExtract(item, centre_y)

            indices = \
                np.where(np.logical_and(objWavelength > Hb_lower,
                         objWavelength < Hb_upper))[0]

            fit_wavelength = objWavelength[indices]

            fit_flux = flux[indices]

            # Now have both the fit flux and wavelength -
            # fit a gaussian to measure the amplitude of emission line

            best_values, covar = self.fitSingleGauss(fit_wavelength,
                                                     fit_flux,
                                                     0.4861 * (1 + z))

            # define some maths helper variables

            amp = abs(best_values['amplitude'])

            num = 100 * np.sqrt(covar[1][1])

            if covar is None:

                x_SN_array.append(0.0)

            elif num / amp > 30 or best_values['sigma'] > 0.007:

                x_SN_array.append(0.0)

            else:

                # Note the 0.00028 is the distance between
                # adjacent wavelength pixels for normalisation

                x_SN_array.append(amp / (noise_value * 0.00028))

        # Now for the y_array

        print 'Moving onto y array'

        for item in y_loop_array:

            flux = objCube.singlePixelExtract(item, centre_x)

            indices = \
                np.where(np.logical_and(objWavelength > Hb_lower,
                         objWavelength < Hb_upper))[0]

            fit_wavelength = objWavelength[indices]

            fit_flux = flux[indices]

            # Now have both the fit flux and wavelength -
            # fit a gaussian to measure the amplitude of emission line

            best_values, covar = self.fitSingleGauss(fit_wavelength,
                                                     fit_flux,
                                                     0.4861 * (1 + z))

            # define some maths helper variables

            amp = abs(best_values['amplitude'])

            num = 100 * np.sqrt(covar[1][1])

            if covar is None:

                y_SN_array.append(0.0)

            elif num / amp > 30 or best_values['sigma'] > 0.007:

                y_SN_array.append(0.0)

            else:

                # Note the 0.00028 is the distance between
                # adjacent wavelength pixels for normalisation

                y_SN_array.append(amp / (noise_value * 0.00028))

        f, (ax1, ax2) = plt.subplots(2, 1, sharex=False, figsize=(18.0, 10.0))

        ax1.plot(x_loop_array, x_SN_array, color='b')
        ax1.set_title('%s x-direction S/N Hb' % objName[26:-5],
                      fontsize=24)
        ax1.set_xlabel(r'x-pixel Position at y = %s' % centre_y, fontsize=24)
        ax1.tick_params(axis='y', which='major', labelsize=15)

        nbins = len(ax1.get_xticklabels())

        ax2.plot(y_loop_array, y_SN_array, color='g')
        ax2.set_xlabel(r'y-pixel Position at x = %s' % centre_x, fontsize=24)
        ax2.tick_params(axis='both', which='major', labelsize=15)
        ax2.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))

        f.subplots_adjust(hspace=0.001)
        f.tight_layout()

        plt.show()

        f.savefig(gal_name[:-5] + '_HbSN.png')

        # Take the median of this array where the flux values don't exceed 50
        # plot an initial evaluation of the spectrum

        f, ax1 = plt.subplots(1, 1, sharex=True, figsize=(18.0, 10.0))

        ax1.plot(objWavelength, fluxArray, color='b')
        ax1.plot(objWavelength, central_spec, color='red')
        ax1.set_title('Sky Spectrum', fontsize=30)
        ax1.set_ylim(-1, 10)
        ax1.tick_params(axis='y', which='major', labelsize=15)
        ax1.set_xlabel(r'Wavelength ($\mu m$)', fontsize=24)
        ax1.set_ylabel(r'Flux', fontsize=24)

        f.tight_layout()

        plt.show()

    def singlePixelExtractMulti_Hb(self, inFile, sci_dir):

        """
        Def:
        Applies the singlePixelExtract function to a list of files,
        in the infile we have the name of the object
        to read in, the redshift and the x,y positions
        of the center of the galaxy from Qfits
        Input: inFile - Containing the object attributes listed above
        """

        # Want to apply the above function to everything within a file
        # First read in the whole file as a Table
        Table = np.loadtxt(inFile, dtype='str')

        for row in Table:

            self.singlePixelExtract_Hb(sci_dir,
                                       row[0],
                                       row[2],
                                       row[3],
                                       row[1],
                                       1)

    def fitSingleGauss(self,
                       wavelength,
                       flux,
                       center):

        """
        Def:
        Supply wavelength and flux arrays as well as a guess
        at the central wavelength value of the gaussian to perform
        a fit and recover the best fit parameters.

        Input: 
                wavelength - wavelength array
                flux - corresponding flux array
                center - central wavelength of emission line in microns

        Output: Best fitting parameters in dictionary

        """
        # Construct the gaussian model from lmfit

        mod = GaussianModel()

        # Guess and set the initial parameter values

        pars = mod.make_params()
        pars['center'].set(center)
        pars['center'].set(vary=False)
        pars['sigma'].set(0.05)
        pars['amplitude'].set(0.8)

        # Perform the model fit

        out = mod.fit(flux, pars, x=wavelength)

        # print out.fit_report()
        # plot an initial evaluation of the model on top of the spectrum

        f, ax1 = plt.subplots(1, 1, sharex=True, figsize=(18.0, 10.0))

        ax1.plot(wavelength, flux, color='b')
        ax1.plot(wavelength, out.best_fit, 'r-')
        ax1.set_title('Object Spectrum', fontsize=30)
        ax1.set_ylim(-10, 10)
        ax1.tick_params(axis='y', which='major', labelsize=15)
        ax1.set_xlabel(r'Wavelength ($\mu m$)', fontsize=24)
        ax1.set_ylabel(r'Flux', fontsize=24)

        f.tight_layout()

        # plt.show()
        plt.close('all')

        return out.best_values, out.covar

    def pSTNK(self,
              object_spectrum,
              z):

        """
        Def:
        Uses the methods in galPhys class to fit a gaussian
        to each of the K-band emission lines and print out
        the signal to noise of each line
        """

        # Create an instance of the galPhysClass with the object spectrum
        galaxy = galPhys(object_spectrum, z)

        # Compute the signal to noise
        galaxy.sToNK()

    def plotHandOII(self,
                    l_o):

        """
        Def:
        Plots the H-band sky spectrum and the location of the OII
        emission line for different objects. Saves with the object name.
        """

        # Read in from the text file containing the wavelength and flux
        Table = np.loadtxt('H_sky.txt')

        wavelength = Table[:, 0] / 1000

        flux = Table[:, 1]

        # Set up a gaussian to simulate the OII emission line
        oii_gauss = GaussianModel()

        pars = oii_gauss.make_params()
        pars['center'].set(l_o)
        pars['sigma'].set(0.0008)
        pars['amplitude'].set(1000)

        init = oii_gauss.eval(pars, x=wavelength)

        # plot an initial evaluation of the model on top of the spectrum
        f, ax1 = plt.subplots(1, 1, sharex=True, figsize=(18.0, 10.0))

        ax1.plot(wavelength, flux, color='b')
        ax1.set_title('H-band OII Position', fontsize=30)
        ax1.plot(wavelength, init, 'k--')
        ax1.axvline(x=l_o, ymin=0, ymax = max(flux), color='r')
        ax1.set_xlim(1.4, 1.8)
        ax1.tick_params(axis='y', which='major', labelsize=15)
        ax1.set_xlabel(r'Wavelength ($\mu m$)', fontsize=24)
        ax1.set_ylabel(r'Flux', fontsize=24)

        f.tight_layout()

        plt.savefig('/disk1/turner/DATA/Gals2/comb/'
                    + 'Science/OII_overplots/gals1_lbg105_OII.png')
        plt.show()

    def stackSpectra(self,
                     fitsList,
                     dL):

        """
        Def:
        Takes a collection of spectra of objects
        in a similar redshift range,
        blueshifts back to rest wavelength and
        then stacks these to create a composite
        spectrum.

        Input: fitsList - file containing the names of
                 the .fits files housing the spectra
               dL - the wavelength separation between
                 adjacent points in the spectra
        Output: compSpec.fits - composite spectrum, spanning the full range of
                 rest frame wavelengths the observations have probed
        """

        # Read in the fitsList filenames and redshift values
        data = np.genfromtxt(fitsList, dtype='str')

        # Save the names and types as lists
        fileNames = data[:, 0]

        redshifts = data[:, 1]

        for item in zip(fileNames, redshifts):
            print item

        # initialise a table and define this as the composite
        initTable = fits.open(fileNames[0])

        # Set the Wavelength and flux, blueshifting
        # the initial wavelength values
        compWavelength = initTable[1].data['Wavelength'] \
            / (1 + float(redshifts[0]))

        compFlux = initTable[1].data['Flux']

        # Start loop over the other files and redshift values
        for i in range(1, len(fileNames)):

            # Open each fits file separately and define
            # the blueshifted wavelength / flux values
            fitsTable = fits.open(fileNames[i])

            fitsWavelength = fitsTable[1].data['Wavelength'] \
                / (1 + float(redshifts[i]))

            fitsFlux = fitsTable[1].data['Flux']

            # IS THE LOWEST WAVELENGTH LOWER
            if fitsWavelength[0] < compWavelength[0]:

                print '[INFO]: The lowest new array wavelength is lower than' \
                    + ' the lowest composite array wavelength: %s ' \
                    % (compWavelength[0])

                # The lowest wavelength is lower
                # (meaning the redshift is higher
                # than the maximum in the composite)
                # Find the highest wavelength index
                # at which this is still the case

                low_counter = 0

                while fitsWavelength[low_counter] < compWavelength[0]:
                    low_counter += 1

                # low_counter is the correct index to clip until.
                # Do so and delete from both the
                # fitsWavelength and fitsFlux arrays
                lowFitsWavelength = fitsWavelength[:low_counter]
                lowFitsFlux = fitsFlux[:low_counter]
                fitsWavelength = fitsWavelength[low_counter:]
                fitsFlux = fitsFlux[low_counter:]

                # Append the lowFitsWavelength and flux
                # arrays to the beginning of the composite spectrum
                compWavelength = np.hstack([lowFitsWavelength, compWavelength])
                compFlux = np.hstack([lowFitsFlux, compFlux])

            # IS THE HIGHEST WAVELENGTH LOWER

            if fitsWavelength[len(fitsWavelength) - 1] > \
                    compWavelength[len(compWavelength) - 1]:

                print '[INFO]: The highest new array wavelength' \
                    + ' is higher than the highest composite array' \
                    + ' wavelength: %s ' \
                    % (compWavelength[len(compWavelength) - 1])

                # The highest wavelength is higher
                # (meaning the redshift is lower
                # than the minimum in the composite)
                # Find the lowest wavelength index
                # at which this is still the case

                high_counter = len(fitsWavelength) - 1

                while fitsWavelength[high_counter] > \
                        compWavelength[len(compWavelength) - 1]:

                    high_counter -= 1

                # high_counter is the correct index to clip until
                # Do so and delete from both the
                # fitsWavelength and fitsFlux arrays

                highFitsWavelength = fitsWavelength[high_counter + 1:]
                highFitsFlux = fitsFlux[high_counter + 1:]
                fitsWavelength = fitsWavelength[:high_counter + 1]
                fitsFlux = fitsFlux[:high_counter + 1]

                # Append the highFitsWavelength and
                # flux arrays to the end of the composite spectrum
                compWavelength = np.hstack([compWavelength,
                                            highFitsWavelength])

                compFlux = np.hstack([compFlux, highFitsFlux])

            # Stage 1 done - maximum and minimum sections added.
            # Remaining stage is to evaluate
            # the points remaining in fitsWavelength and fitsFlux 1 by 1
            # and stack the fluxes together in bins

            # Check all the remaining wavelength values
            print 'Stacking remaining wavelengths'

            for i in range(len(fitsWavelength)):

                counter = -1

                for j in range(len(compWavelength)):

                    counter += 1

                    if compWavelength[j] - dL <= \
                       fitsWavelength[i] <= \
                       compWavelength[j]:

                        break

                compFlux[j] = compFlux[j] + fitsFlux[i]

        # plot an initial evaluation of the spectrum
        f, ax1 = plt.subplots(1, 1, sharex=True, figsize=(18.0, 10.0))
        ax1.plot(compWavelength, compFlux, color='b')
        ax1.set_title('Composite Spectrum', fontsize=30)
        ax1.set_ylim(-1, 50)
        ax1.tick_params(axis='y', which='major', labelsize=15)
        ax1.set_xlabel(r'Wavelength ($\mu m$)', fontsize=24)
        ax1.set_ylabel(r'Flux', fontsize=24)
        f.tight_layout()
        plt.savefig('/disk1/turner/DATA/comp_spectrum.png')
        plt.show()

        # Isolate the Hbeta and OIII lines and save to a new spectrum
        index = np.where(np.logical_and(compWavelength > 0.485,
                                        compWavelength < 0.510))[0]

        saveWavelength = compWavelength[index]

        saveFlux = compFlux[index]

        tbhdu = fits.new_table(fits.ColDefs(
            [fits.Column(name='Wavelength', format='E', array=saveWavelength),
             fits.Column(name='Flux', format='E', array=saveFlux)]))

        prihdu = fits.PrimaryHDU(header=initTable[0].header)

        thdulist = fits.HDUList([prihdu, tbhdu])

        thdulist.writeto('/disk1/turner/DATA/comp_spectrum.fits', clobber=True)

    def av_seeing(self,
                  inFile):

        """
        Def:
        Compute the average seeing from the combine_input.txt file.
        Input: inFile - the combine_input.txt file
        """

        # Read in the infile
        data = np.genfromtxt(inFile, dtype='str')

        # Save the names and types as lists
        av_seeing = data[:, 5]

        ar = []

        for item in av_seeing:

            ar.append(float(item))

        print np.nanmedian(ar)

    def seeing_better_than(self,
                           inFile,
                           seeing_value):

        """
        Def:
        Print how many objects in a sample have seeing
        better than a chosen value
        Input: inFile - the combine_input.txt file
                seeing_value - the value the seeing must exceed

        """

        # Read in the infile
        data = np.genfromtxt(inFile, dtype='str')

        # assign the seeing column
        av_seeing = data[:, 5]

        print 'This is the length of the seeing array: %s' \
              % len(np.unique(av_seeing))

        # create an array to store the results
        ar = []

        # loop through and if better than seeing value, store it
        for item in av_seeing:

            if float(item) < seeing_value and float(item) > 0:

                ar.append(float(item))

        # now search for unique values in the list
        print 'These are the unique seeing values better than: %s %s' \
              % (seeing_value, np.unique(ar))

        print 'We have lost %f objects' \
              % (len(np.unique(av_seeing)) - len(np.unique(ar)))

        print 'Therefore the total on source time is: %.2f hours' \
              % ((5 * (len(np.unique(av_seeing)) - (len(np.unique(av_seeing))
                 - len(np.unique(ar))))) / 60.0)

    def multi_plot_HK_sn_map(self,
                             infile):

        """
        Def:
        Plot the sn maps for lots of cubes together.

        Input:
                infile - file listing the redshifts and cubenames
        """

        # read in the table of cube names
        Table = ascii.read(infile)

        for entry in Table:

            obj_name = entry[0]

            cube = cubeOps(obj_name)

            redshift = entry[1]

            print "\nDoing %s (redshift = %.3f) ..." % (obj_name, redshift)

            cube.plot_HK_sn_map(redshift, savefig=True)

    def multi_plot_K_sn_map(self,
                            infile):

        """
        Def:
        Plot the sn maps for lots of cubes together.

        Input:
                infile - file listing the redshifts and cubenames
        """

        # read in the table of cube names
        Table = ascii.read(infile)

        for entry in Table:

            obj_name = entry[0]

            cube = cubeOps(obj_name)

            redshift = entry[1]

            print "\nDoing %s (redshift = %.3f) ..." % (obj_name, redshift)

            cube.plot_K_sn_map(redshift, savefig=True)

    def multi_plot_HK_image(self,
                            infile):

        """
        Def:
        Plot the sn maps for lots of cubes together.

        Input:
                infile - file listing the redshifts and cubenames
        """

        # read in the table of cube names
        Table = ascii.read(infile)

        for entry in Table:

            obj_name = entry[0]

            cube = cubeOps(obj_name)

            redshift = entry[1]

            print "\nDoing %s (redshift = %.3f) ..." % (obj_name, redshift)

            cube.plot_HK_image(redshift, savefig=True)

    def multi_plot_K_image(self,
                           infile):

        """
        Def:
        Plot the sn maps for lots of cubes together.

        Input:
                infile - file listing the redshifts and cubenames
        """

        # read in the table of cube names
        Table = ascii.read(infile)

        for entry in Table:

            obj_name = entry[0]

            cube = cubeOps(obj_name)

            redshift = entry[1]

            print "\nDoing %s (redshift = %.3f) ..." % (obj_name, redshift)

            cube.plot_K_image(redshift, savefig=True)

    def multi_plot_OIII_vel_map(self,
                                infile):

        """
        Def:

        Plot the velocity maps for lots of cubes together.

        Input:
                infile - file listing the redshifts and cubenames

        """

        # read in the table of cube names

        Table = ascii.read(infile)

        for entry in Table:

            obj_name = entry[0]

            cube = cubeOps(obj_name)

            redshift = entry[1]

            print "\nDoing %s (redshift = %.3f) ..." % (obj_name, redshift)

            cube.OIII_vel_map(redshift, savefig=True)

    def multi_plot_all_maps(self,
                            infile,
                            binning,
                            xbin,
                            ybin,
                            interp):
        """
        Def:
        Plot the velocity maps for lots of cubes together.

        Input:
                infile - file listing the redshifts and cubenames
                      and also the x,y central position of the galaxy
                     also needs to contain the 'standard star cube'
                    and an associated 'skycube'
               binning - whether to bin velocity field data or not
               xbin - what the x-direction bin size should be
               ybin - what the y-direction bin size should be
               interp - what the bin interpolation should be
        """

        # read in the table of cube names
        Table = ascii.read(infile)

        # assign variables to the different items in the infile
        for entry in Table:

            obj_name = entry[0]

            cube = cubeOps(obj_name)

            redshift = entry[1]

            centre_x = entry[2]

            centre_y = entry[3]

            std_cube = entry[4]

            sky_cube = entry[5]

            # define the science directory for each cube
            sci_dir = obj_name[:len(obj_name) - obj_name[::-1].find("/") - 1]

            print "\nDoing %s (redshift = %.3f) ..." % (obj_name, redshift)

            # check whether we're looking at HK or K band

            if cube.filter == 'HK':

                # we're in the HK band, use the methods which include OII
                # and construct a 3 x 2 grid of plots
                fig, axes = plt.subplots(figsize=(14, 14), nrows=3, ncols=3)
                # fig.subplots_adjust(right=0.83)
                # cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])

                # now this is set up, get the data to populate the plots
                sn_dict = cube.plot_HK_sn_map(redshift, savefig=True)

                # these are the sn images of the lines
                # add these to the first
                for grid, ax in zip(sn_dict.values(), axes[0]):

                    # add the plot

                    im = ax.imshow(grid, aspect='auto', vmin=0.,
                                   vmax=3.,
                                   cmap=plt.get_cmap('hot'))

                    # add colourbar to each plot
                    divider = make_axes_locatable(ax)

                    cax_new = divider.append_axes('right',
                                                  size='10%',
                                                  pad=0.05)

                    plt.colorbar(im, cax=cax_new)

                # name the plots
                for name, ax in zip(sn_dict.keys(), axes[0]):

                    ax.set_title('%s_image' % name)

                # top row populated, now add the OII & Hb
                # metallicity maps and the velocity map

                # metallicity maps
                Hb_met_array, OII_met_array \
                    = cube.plot_HK_image(redshift, savefig=True)

                # plot both of these
                # Hb
                im = axes[1][1].imshow(Hb_met_array,
                                       aspect='auto',
                                       vmin=7.5,
                                       vmax=9.0,
                                       cmap=plt.get_cmap('jet'))

                # add colourbar to each plot
                divider = make_axes_locatable(axes[1][1])
                cax_new = divider.append_axes('right', size='10%', pad=0.05)
                plt.colorbar(im, cax=cax_new)

                # set the name
                axes[1][1].set_title('OIII / Hb Metallicity')

                # OII
                im = axes[1][0].imshow(OII_met_array,
                                       aspect='auto',
                                       vmin=7.5,
                                       vmax=9.0,
                                       cmap=plt.get_cmap('jet'))

                # add colourbar to each plot
                divider = make_axes_locatable(axes[1][0])
                cax_new = divider.append_axes('right', size='10%', pad=0.05)
                plt.colorbar(im, cax=cax_new)

                # set the title
                axes[1][0].set_title('OIII / OII Metallicity')

                # now want to get the initial gaussian parameters to feed
                # into the OIII, Hb and OII velocity maps

                # first optimally extract the spectrum using
                # no binning for the wavelength axis (could do this)
                spectrum, wl = self.galExtract(sci_dir,
                                               std_cube,
                                               obj_name,
                                               sky_cube,
                                               centre_y,
                                               centre_x,
                                               1)

                # now feed the output of this into fit_lines_HK
                oiii_values, hb_values, oii_values \
                    = self.fit_lines_HK(spectrum,
                                        wl,
                                        redshift=redshift)

                # OIII velocity map - using the initial gaussian
                # parameters determined in the above analysis
                OIII_vel, OIII_disp = cube.OIII_vel_map(redshift,
                                                        binning=binning,
                                                        savefig=True,
                                                        xbin=xbin,
                                                        ybin=ybin,
                                                        interp=interp,
                                                        **oiii_values)

                # use the velocity map values to get the limits of the
                # colourbar using np.nanpercentile
                try:

                    mn_vel, mx_vel = np.nanpercentile(OIII_vel,
                                                      [10.0, 90.0])

                    mn_disp, mx_disp = np.nanpercentile(OIII_disp,
                                                        [10.0, 90.0])

                except TypeError:
                    mn_vel, mx_vel = [-100, 100]
                    mn_disp, mx_disp = [0, 100]

                # mx = np.nanmin([mx, -mn])

                im = axes[1][2].imshow(OIII_vel, aspect='auto', vmin=mn_vel,
                                       vmax=mx_vel,
                                       interpolation='nearest',
                                       cmap=plt.get_cmap('jet'))

                # add colourbar to each plot
                divider = make_axes_locatable(axes[1][2])
                cax_new = divider.append_axes('right', size='10%', pad=0.05)
                plt.colorbar(im, cax=cax_new)

                # set the title
                axes[1][2].set_title('OIII Velocity')

                im = axes[2][2].imshow(OIII_disp, aspect='auto', vmin=mn_disp,
                                       vmax=mx_disp,
                                       interpolation='nearest',
                                       cmap=plt.get_cmap('jet'))

                # add colourbar to each plot
                divider = make_axes_locatable(axes[2][2])
                cax_new = divider.append_axes('right', size='10%', pad=0.05)
                plt.colorbar(im, cax=cax_new)

                # set the title
                axes[2][2].set_title('OIII Dispersion')

                # also include OII and Hb velocity for fun

                # OII velocity map
                OII_vel, OII_disp = cube.OII_vel_map(redshift,
                                                     binning=binning,
                                                     savefig=True,
                                                     xbin=xbin,
                                                     ybin=ybin,
                                                     interp=interp,
                                                     **oii_values)

                # use the velocity map values to get the limits of the
                # colourbar using np.nanpercentile
                try:
                    mn_vel, mx_vel = np.nanpercentile(OII_vel, [10.0, 90.0])
                    mn_disp, mx_disp = np.nanpercentile(OII_disp, [10.0, 90.0])

                except TypeError:
                    mn_vel, mx_vel = [-100, 100]
                    mn_disp, mx_disp = [0, 100]
                # mx = np.nanmin([mx, -mn])

                im = axes[2][0].imshow(OII_vel, aspect='auto', vmin=mn_vel,
                                       vmax=mx_vel,
                                       interpolation='nearest',
                                       cmap=plt.get_cmap('jet'))

                # add colourbar to each plot
                divider = make_axes_locatable(axes[2][0])
                cax_new = divider.append_axes('right', size='10%', pad=0.05)
                plt.colorbar(im, cax=cax_new)

                # set the title
                axes[2][0].set_title('OII Velocity')

                # Hb velocity map
                Hb_vel, Hb_disp = cube.Hb_vel_map(redshift,
                                                  binning=binning,
                                                  savefig=True,
                                                  xbin=xbin,
                                                  ybin=ybin,
                                                  interp=interp,
                                                  **hb_values)

                # use the velocity map values to get the limits of the
                # colourbar using np.nanpercentile
                try:

                    mn_vel, mx_vel = np.nanpercentile(Hb_vel, [10.0, 90.0])
                    mn_disp, mx_disp = np.nanpercentile(Hb_disp, [10.0, 90.0])

                except TypeError:

                    mn_vel, mx_vel = [-100, 100]
                    mn_disp, mx_disp = [0, 100]
                # mx = np.nanmin([mx, -mn])

                im = axes[2][1].imshow(Hb_vel, aspect='auto', vmin=mn_vel,
                                       vmax=mx_vel,
                                       interpolation='nearest',
                                       cmap=plt.get_cmap('jet'))

                # add colourbar to each plot
                divider = make_axes_locatable(axes[2][1])
                cax_new = divider.append_axes('right', size='10%', pad=0.05)
                plt.colorbar(im, cax=cax_new)

                # set the title
                axes[2][1].set_title('Hb Velocity')

                # save the big figure
                # fig.show()
                fig.savefig('%s_all_maps.pdf' % obj_name[:-5])
                plt.close('all')

            elif cube.filter == 'K':

                # we're in the K band, use the methods which dont include OII
                # and construct a 2 x 2 grid of plots
                fig, axes = plt.subplots(figsize=(10, 12), nrows=3, ncols=2)

                # now this is set up, get the data to populate the plots
                sn_dict = cube.plot_K_sn_map(redshift, savefig=True)

                # these are the sn images of the lines
                # add these to the first
                for grid, ax in zip(sn_dict.values(), axes[0]):

                    # add the plot

                    im = ax.imshow(grid, aspect='auto', vmin=0.,
                                   vmax=3.,
                                   cmap=plt.get_cmap('hot'))

                    # add colourbar to each plot
                    divider = make_axes_locatable(ax)
                    cax_new = divider.append_axes('right',
                                                  size='10%',
                                                  pad=0.05)
                    plt.colorbar(im, cax=cax_new)

                # name the plots
                for name, ax in zip(sn_dict.keys(), axes[0]):
                    ax.set_title('%s_image' % name)

                # top row populated, now add the OII & Hb
                # metallicity maps and the velocity map

                # metallicity maps
                Hb_met_array \
                    = cube.plot_K_image(redshift, savefig=True)

                # plot this
                # Hb
                im = axes[1][0].imshow(Hb_met_array, aspect='auto', vmin=7.5,
                                       vmax=9.0,
                                       cmap=plt.get_cmap('jet'))

                # add colourbar to each plot
                divider = make_axes_locatable(axes[1][0])
                cax_new = divider.append_axes('right', size='10%', pad=0.05)
                plt.colorbar(im, cax=cax_new)

                axes[1][0].set_title('OIII / Hb Metallicity')

                # now want to get the initial gaussian parameters to feed
                # into the OIII, Hb and OII velocity maps

                # first optimally extract the spectrum using
                # no binning for the wavelength axis (could do this)
                spectrum, wl = self.galExtract(sci_dir,
                                               std_cube,
                                               obj_name,
                                               sky_cube,
                                               centre_y,
                                               centre_x,
                                               1)

                # now feed the output of this into fit_lines_HK
                oiii_values, hb_values \
                    = self.fit_lines_K(spectrum,
                                       wl,
                                       redshift=redshift)

                # OIII velocity map - using the initial gaussian
                # parameters determined in the above analysis
                OIII_vel, OIII_disp = cube.OIII_vel_map(redshift,
                                                        binning=binning,
                                                        savefig=True,
                                                        xbin=xbin,
                                                        ybin=ybin,
                                                        interp=interp,
                                                        **oiii_values)

                # use the velocity map values to get the limits of the
                # colourbar using np.nanpercentile
                try:

                    mn_vel, mx_vel = np.nanpercentile(OIII_vel,
                                                      [10.0, 90.0])

                    mn_disp, mx_disp = np.nanpercentile(OIII_disp,
                                                        [10.0, 90.0])

                except TypeError:
                    mn_vel, mx_vel = [-100, 100]
                    mn_disp, mx_disp = [0, 100]
                # mx = np.nanmin([mx, -mn])

                im = axes[1][1].imshow(OIII_vel,
                                       aspect='auto',
                                       vmin=mn_vel,
                                       vmax=mx_vel,
                                       interpolation='nearest',
                                       cmap=plt.get_cmap('jet'))

                # add colourbar to each plot
                divider = make_axes_locatable(axes[1][1])
                cax_new = divider.append_axes('right', size='10%', pad=0.05)
                plt.colorbar(im, cax=cax_new)

                # set the title
                axes[1][1].set_title('OIII Velocity')

                im = axes[2][1].imshow(OIII_disp,
                                       aspect='auto',
                                       vmin=mn_disp,
                                       vmax=mx_disp,
                                       interpolation='nearest',
                                       cmap=plt.get_cmap('jet'))

                # add colourbar to each plot
                divider = make_axes_locatable(axes[2][1])
                cax_new = divider.append_axes('right', size='10%', pad=0.05)
                plt.colorbar(im, cax=cax_new)

                # set the title
                axes[2][1].set_title('OIII Dispersion')

                # Now add the Hb velocity for fun
                Hb_vel, Hb_disp = cube.Hb_vel_map(redshift,
                                                  binning=binning,
                                                  savefig=True,
                                                  xbin=xbin,
                                                  ybin=ybin,
                                                  interp=interp,
                                                  **hb_values)

                # use the velocity map values to get the limits of the
                # colourbar using np.nanpercentile

                try:
                    mn_vel, mx_vel = np.nanpercentile(Hb_vel, [10.0, 90.0])
                    mn_disp, mx_disp = np.nanpercentile(Hb_disp, [10.0, 90.0])

                except TypeError:
                    mn_vel, mx_vel = [-100, 100]
                    mn_disp, mx_disp = [0, 100]
                # mx = np.nanmin([mx, -mn])

                im = axes[2][0].imshow(Hb_vel, aspect='auto', vmin=mn_vel,
                                       vmax=mx_vel,
                                       interpolation='nearest',
                                       cmap=plt.get_cmap('jet'))

                # add colourbar to each plot
                divider = make_axes_locatable(axes[2][0])
                cax_new = divider.append_axes('right', size='10%', pad=0.05)
                plt.colorbar(im, cax=cax_new)

                # set the title
                axes[2][0].set_title('Hb Velocity')

                # save the big figure
                # plt.show()
                fig.savefig('%s_all_maps.pdf' % obj_name[:-5])
                plt.close('all')

    def fit_lines_K(self,
                    spectrum,
                    wl,
                    redshift):
        """
        Def:
        Takes the 1D spectrum output from galExtract and fits the two
        K-band emission lines, returning their parameters in a dict
        Input: spectrum - from galExtract
               redshift - the redshift of the galaxy being fitted
               wl - corresponding 1D wavelength array
        output: Dictionaries containing the gaussian parameters
        """
        # the spectrum is a 1D python array
        # first determine the OIII and Hb central wavelengths

        oiii_central = 0.500824 * (1. + redshift)
        hb_central = 0.486268 * (1. + redshift)

        # first dealing with OIII
        # trying to isolate the wavelength region around the emission line

        line_idx = np.argmin(np.abs(wl - oiii_central))

        fit_wl = wl[line_idx - 8: line_idx + 8]
        fit_flux = spectrum[line_idx - 8: line_idx + 8]

        # construct gaussian model using lmfit
        gmod = GaussianModel()
        # set the initial parameter values
        pars = gmod.make_params()

        pars['center'].set(value=oiii_central)
        pars['sigma'].set(0.0008)
        pars['amplitude'].set(1E-20)

        # perform the fit
        out = gmod.fit(fit_flux, pars, x=fit_wl)

        # plot the results
#        f, ax1 = plt.subplots(1, 1, sharex=True, figsize=(18.0, 10.0))
#        ax1.plot(fit_wl, fit_flux, color='b')
#        ax1.plot(fit_wl, out.best_fit, 'r-')
#        ax1.set_title('Object Spectrum', fontsize=30)
#        #ax1.set_ylim(0, 4)
#        ax1.tick_params(axis='y', which='major', labelsize=15)
#        ax1.set_xlabel(r'Wavelength ($\mu m$)', fontsize=24)
#        ax1.set_ylabel(r'Flux', fontsize=24)
#        f.tight_layout()
#        plt.show()

        # print out.fit_report()

        # create an oiii dictionary to house the best parameters
        oiii_values = {'centre_oiii': out.best_values['center'],
                       'sigma_oiii': out.best_values['sigma'],
                       'amplitude_oiii': out.best_values['amplitude']}

        # now dealing with hb
        # trying to isolate the wavelength region around the emission line

        line_idx = np.argmin(np.abs(wl - hb_central))

        fit_wl = wl[line_idx - 8: line_idx + 8]
        fit_flux = spectrum[line_idx - 8: line_idx + 8]

        # construct gaussian model using lmfit
        gmod = GaussianModel()
        # set the initial parameter values
        pars = gmod.make_params()

        pars['center'].set(value=hb_central)
        pars['sigma'].set(0.0008)
        pars['amplitude'].set(1E-20)

        # perform the fit
        out = gmod.fit(fit_flux, pars, x=fit_wl)

#        # plot the results
#        f, ax1 = plt.subplots(1, 1, sharex=True, figsize=(18.0, 10.0))
#        ax1.plot(fit_wl, fit_flux, color='b')
#        ax1.plot(fit_wl, out.best_fit, 'r-')
#        ax1.set_title('Object Spectrum', fontsize=30)
#        #ax1.set_ylim(0, 4)
#        ax1.tick_params(axis='y', which='major', labelsize=15)
#        ax1.set_xlabel(r'Wavelength ($\mu m$)', fontsize=24)
#        ax1.set_ylabel(r'Flux', fontsize=24)
#        f.tight_layout()
#        plt.show()

        # create an oiii dictionary to house the best parameters
        hb_values = {'centre_hb': out.best_values['center'],
                     'sigma_hb': out.best_values['sigma'],
                     'amplitude_hb': out.best_values['amplitude']}

        # simply return both of these dictionaries as a good
        # indication for what the initial gaussian
        # parameters / wavelength values should be

        return oiii_values, hb_values

    def fit_lines_HK(self,
                     spectrum,
                     wl,
                     redshift):

        """
        Def:
        Takes the 1D spectrum output from galExtract and fits the three
        HK-band emission lines, returning their parameters in a dict
        Input: spectrum - from galExtract
               redshift - the redshift of the galaxy being fitted
               wl - corresponding 1D wavelength array
        output: Dictionaries containing the gaussian parameters
        """
        # the spectrum is a 1D python array
        # first determine the OIII and Hb central wavelengths

        oiii_central = 0.500824 * (1. + redshift)
        hb_central = 0.486268 * (1. + redshift)
        oii_central = 0.3729875 * (1. + redshift)

        # first dealing with OIII
        # trying to isolate the wavelength region around the emission line

        line_idx = np.argmin(np.abs(wl - oiii_central))

        fit_wl = wl[line_idx - 8: line_idx + 8]
        fit_flux = spectrum[line_idx - 8: line_idx + 8]

        # construct gaussian model using lmfit
        gmod = GaussianModel()
        # set the initial parameter values
        pars = gmod.make_params()

        pars['center'].set(value=oiii_central)
        pars['sigma'].set(0.0008)
        pars['amplitude'].set(1E-20)

        # perform the fit
        out = gmod.fit(fit_flux, pars, x=fit_wl)

#        # plot the results
#        f, ax1 = plt.subplots(1, 1, sharex=True, figsize=(18.0, 10.0))
#        ax1.plot(fit_wl, fit_flux, color='b')
#        ax1.plot(fit_wl, out.best_fit, 'r-')
#        ax1.set_title('Object Spectrum', fontsize=30)
#        #ax1.set_ylim(0, 4)
#        ax1.tick_params(axis='y', which='major', labelsize=15)
#        ax1.set_xlabel(r'Wavelength ($\mu m$)', fontsize=24)
#        ax1.set_ylabel(r'Flux', fontsize=24)
#        f.tight_layout()
#        plt.show()

        # create an oiii dictionary to house the best parameters
        oiii_values = {'centre_oiii': out.best_values['center'],
                       'sigma_oiii': out.best_values['sigma'],
                       'amplitude_oiii': out.best_values['amplitude']}

        # now dealing with hb
        # trying to isolate the wavelength region around the emission line

        line_idx = np.argmin(np.abs(wl - hb_central))

        fit_wl = wl[line_idx - 8: line_idx + 8]
        fit_flux = spectrum[line_idx - 8: line_idx + 8]

        # construct gaussian model using lmfit
        gmod = GaussianModel()
        # set the initial parameter values
        pars = gmod.make_params()

        pars['center'].set(value=hb_central)
        pars['sigma'].set(0.0008)
        pars['amplitude'].set(1E-20)

        # perform the fit
        out = gmod.fit(fit_flux, pars, x=fit_wl)

#        # plot the results
#        f, ax1 = plt.subplots(1, 1, sharex=True, figsize=(18.0, 10.0))
#        ax1.plot(fit_wl, fit_flux, color='b')
#        ax1.plot(fit_wl, out.best_fit, 'r-')
#        ax1.set_title('Object Spectrum', fontsize=30)
#        #ax1.set_ylim(0, 4)
#        ax1.tick_params(axis='y', which='major', labelsize=15)
#        ax1.set_xlabel(r'Wavelength ($\mu m$)', fontsize=24)
#        ax1.set_ylabel(r'Flux', fontsize=24)
#        f.tight_layout()
#        plt.show()

        # create an oiii dictionary to house the best parameters
        hb_values = {'centre_hb': out.best_values['center'],
                     'sigma_hb': out.best_values['sigma'],
                     'amplitude_hb': out.best_values['amplitude']}

        # simply return both of these dictionaries as a good
        # indication for what the initial gaussian
        # parameters / wavelength values should be

        # finally dealing with OII
        # trying to isolate the wavelength region around the emission line

        line_idx = np.argmin(np.abs(wl - oii_central))

        fit_wl = wl[line_idx - 8: line_idx + 8]
        fit_flux = spectrum[line_idx - 8: line_idx + 8]

        # construct gaussian model using lmfit
        gmod = GaussianModel()
        # set the initial parameter values
        pars = gmod.make_params()

        pars['center'].set(value=oii_central)
        pars['sigma'].set(0.0008)
        pars['amplitude'].set(1E-20)

        # perform the fit
        out = gmod.fit(fit_flux, pars, x=fit_wl)

#        # plot the results
#        f, ax1 = plt.subplots(1, 1, sharex=True, figsize=(18.0, 10.0))
#        ax1.plot(fit_wl, fit_flux, color='b')
#        ax1.plot(fit_wl, out.best_fit, 'r-')
#        ax1.set_title('Object Spectrum', fontsize=30)
#        #ax1.set_ylim(0, 4)
#        ax1.tick_params(axis='y', which='major', labelsize=15)
#        ax1.set_xlabel(r'Wavelength ($\mu m$)', fontsize=24)
#        ax1.set_ylabel(r'Flux', fontsize=24)
#        f.tight_layout()
#        plt.show()

        # create an oiii dictionary to house the best parameters
        oii_values = {'centre_oii': out.best_values['center'],
                      'sigma_oii': out.best_values['sigma'],
                      'amplitude_oii': out.best_values['amplitude']}

        return oiii_values, hb_values, oii_values

    """
    #####################################################################

    Copyright (C) 2001-2014, Michele Cappellari
    E-mail: michele.cappellari_at_physics.ox.ac.uk

    Updated versions of the software are available from my web page
    http://purl.org/cappellari/software

    If you have found this software useful for your
    research, we would appreciate an acknowledgment to use of
    `the Voronoi binning method by Cappellari & Copin (2003)'.

    This software is provided as is without any warranty whatsoever.
    Permission to use, for non-commercial purposes is granted.
    Permission to modify for personal or internal use is granted,
    provided this copyright and disclaimer are included unchanged
    at the beginning of the file. All other rights are reserved.

    #####################################################################

    NAME:
        VORONOI_2D_BINNING

    AUTHOR:
          Michele Cappellari, University of Oxford
          michele.cappellari_at_physics.ox.ac.uk

    PURPOSE:
          Perform adaptive spatial binning of Integral-Field Spectroscopic
          (IFS) data to reach a chosen constant signal-to-noise ratio per bin.
          This method is required for the proper analysis of IFS
          observations, but can also be used for standard photometric
          imagery or any other two-dimensional data.
          This program precisely implements the algorithm described in
          section 5.1 of the reference below.

    EXPLANATION:
          Further information on VORONOI_2D_BINNING algorithm can be found in
          Cappellari M., Copin Y., 2003, MNRAS, 342, 345
          http://adsabs.harvard.edu/abs/2003MNRAS.342..345C

    CALLING SEQUENCE:

        binNum, xBin, yBin, xBar, yBar, sn, nPixels, scale = voronoi_2d_binning(
                    x, y, signal, noise, targetSN, plot=True, quiet=False,
                    wvt=False, cvt=True, pixelsize=None)

        The function _sn_func() below returns the S/N of a bin and it can be
        changed by the user if needed.

    INPUTS:
               X: Vector containing the X coordinate of the pixels to bin.
                  Arbitrary units can be used (e.g. arcsec or pixels).
                  In what follows the term "pixel" refers to a given
                  spatial element of the dataset (sometimes called "spaxel" in
                  the IFS community): it can be an actual pixel of a CCD
                  image, or a spectrum position along the slit of a long-slit
                  spectrograph or in the field of view of an IFS
                  (e.g. a lenslet or a fiber).
                  It is assumed here that pixels are arranged in a regular
                  grid, so that the pixel size is a well defined quantity.
                  The pixel grid however can contain holes (some pixels can be
                  excluded from the binning) and can have an irregular boundary.
                  See the above reference for an example and details.
               Y: Vector (same size as X) containing the Y coordinate
                  of the pixels to bin.
          SIGNAL: Vector (same size as X) containing the signal
                  associated with each pixel, having coordinates (X,Y).
                  If the `pixels' are actually the apertures of an
                  integral-field spectrograph, then the signal can be
                  defined as the average flux in the spectral range under
                  study, for each aperture.
                  If pixels are the actual pixels of the CCD in a galaxy
                  image, the signal will be simply the counts in each pixel.
           NOISE: Vector (same size as X) containing the corresponding
                  noise (1 sigma error) associated with each pixel.
        TARGETSN: The desired signal-to-noise ratio in the final
                  2D-binned data. E.g. a S/N~50 per pixel may be a
                  reasonable value to extract stellar kinematics
                  information from galaxy spectra.

    KEYWORDS:
             CVT: Set this keyword to skip the Centroidal Voronoi Tessellation
                  (CVT) step (vii) of the algorithm in Section 5.1 of
                  Cappellari & Copin (2003).
                  This may be useful if the noise is strongly non Poissonian,
                  the pixels are not optimally weighted, and the CVT step
                  appears to introduces significant gradients in the S/N.
                  A similar alternative consists of using the /WVT keyword below.
            PLOT: Set this keyword to produce a plot of the two-dimensional
                  bins and of the corresponding S/N at the end of the
                  computation.
         PIXSIZE: Optional pixel scale of the input data.
                  This can be the size of a pixel of an image or the size
                  of a spaxel or lenslet in an integral-field spectrograph.
                - The value is computed automatically by the program, but
                  this can take a long times when (X, Y) have many elements.
                  In those cases the PIXSIZE keyword should be given.
           QUIET: by default the program shows the progress while accreting
                  pixels and then while iterating the CVT. Set this keyword
                  to avoid printing progress results.
             WVT: When this keyword is set, the routine bin2d_cvt_equal_mass is
                  modified as proposed by Diehl & Statler (2006, MNRAS, 368, 497).
                  In this case the final step of the algorithm, after the bin-accretion
                  stage, is not a modified Centroidal Voronoi Tessellation, but it uses
                  a Weighted Voronoi Tessellation.
                  This may be useful if the noise is strongly non Poissonian,
                  the pixels are not optimally weighted, and the CVT step
                  appears to introduces significant gradients in the S/N.
                  A similar alternative consists of using the /NO_CVT keyword above.
                  If you use the /WVT keyword you should also include a reference to
                  `the WVT modification proposed by Diehl & Statler (2006).'

    OUTPUTS:
       BINNUMBER: Vector (same size as X) containing the bin number assigned
                  to each input pixel. The index goes from zero to Nbins-1.
                  This vector alone is enough to make *any* subsequent
                  computation on the binned data. Everything else is optional!
            XBIN: Vector (size Nbins) of the X coordinates of the bin generators.
                  These generators uniquely define the Voronoi tessellation.
            YBIN: Vector (size Nbins) of Y coordinates of the bin generators.
            XBAR: Vector (size Nbins) of X coordinates of the bins luminosity
                  weighted centroids. Useful for plotting interpolated data.
            YBAR: Vector (size Nbins) of Y coordinates of the bins luminosity
                  weighted centroids.
              SN: Vector (size Nbins) with the final SN of each bin.
         NPIXELS: Vector (size Nbins) with the number of pixels of each bin.
           SCALE: Vector (size Nbins) with the scale length of the Weighted
                  Voronoi Tessellation, when the /WVT keyword is set.
                  In that case SCALE is *needed* together with the coordinates
                  XBIN and YBIN of the generators, to compute the tessellation
                  (but one can also simply use the BINNUMBER vector).

    PROCEDURES USED:
          The following procedures are contained in the main VORONOI_2D_BINNING program.
              SN_FUNC           -- Example routine to calculate the S/N of a bin.
              WEIGHTED_CENTROID -- computes weighted centroid of one bin
              BIN_ROUNDNESS     -- equation (5) of Cappellari & Copin (2003)
              BIN_ACCRETION     -- steps (i)-(v) in section 5.1
              REASSIGN_BAD_BINS -- steps (vi)-(vii) in section 5.1
              CVT_EQUAL_MASS    -- the modified Lloyd algorithm in section 4.1
              COMPUTE_USEFUL_BIN_QUANTITIES -- self explanatory
              DISPLAY_PIXELS    -- plotting of colored pixels

    MODIFICATION HISTORY:
          V1.0.0: First implementation. Michele Cappellari, Leiden, June 2001
          V2.0.0: Major revisions. Stable version. MC, Leiden, 11 September 2001
          V2.1.0: First released version. Written documentation.
              MC, Vicenza, 13 February 2003
          V2.2.0: Added computation of useful bin quantities in output. Deleted some
              safety checks for zero size bins in CVT. Minor polishing of the code.
              MC, Leiden, 11 March 2003
          V2.3.0: Unified the three tests to stop the accretion of one bin.
              This can improve some bins at the border. MC, Leiden, 9 April 2003
          V2.3.1: Do *not* assume the first bin is made of one single pixel.
              Added computation of S/N scatter and plotting of 1-pixel bins.
              MC, Leiden, 13 April 2003
          V2.4.0: Addedd basic error checking of input S/N. Reintroduced the
              treatment for zero-size bins in CVT, which was deleted in V2.2.
              Thanks to Robert Sharp and Kambiz Fathi for reporting problems.
              MC, Leiden, 10 December 2003.
          V2.4.1: Added /QUIET keyword and verbose output during the computation.
              After suggestion by Richard McDermid. MC, Leiden, 14 December 2003
          V2.4.2: Use LONARR instead of INTARR to define the CLASS vector,
              to be able to deal with big images. Thanks to Tom Statler.
              MC, Leiden, 4 August 2004
          V2.4.3: Corrected bug introduced in version 2.3.1. It went undetected
              for a long time because it could only happen in special conditions.
              Now we recompute the index of the good bins after computing all
              centroids of the reassigned bins in reassign_bad_bins. Many thanks
              to Simona Ghizzardi for her clear analysis of the problem and
              the solution. MC, Leiden, 29 November 2004
          V2.4.4: Prevent division by zero for pixels with signal=0
              and noise=sqrt(signal)=0, as can happen from X-ray data.
              MC, Leiden, 30 November 2004
          V2.4.5: Added BIN2D prefix to internal routines to avoid possible
              naming conflicts. MC, Leiden, 3 December 2004
          V2.4.6: Added /NO_CVT keyword to optionally skip the CVT step of
              the algorithm. MC, Leiden, 27 August 2005
          V2.4.7: Verify that SIGNAL and NOISE are non negative vectors.
              MC, Leiden, 27 September 2005
          V2.4.8: Use geometric centroid of a bin during the bin-accretion stage,
              to allow the routine to deal with negative signal (e.g. in
              background-subtracted X-ray images). Thanks to Steven Diehl for
              pointing out the usefulness of dealing with negative signal.
              MC, Leiden, 23 December 2005
          V2.5.0: Added two new lines of code and the corresponding /WVT keyword
              to implement the nice modification to the algorithm proposed by
              Diehl & Statler (2006). MC, Leiden, 9 March 2006
          V2.5.1: Updated documentation. MC, Oxford, 3 November 2006
          V2.5.2: Print number of unbinned pixels. MC, Oxford, 28 March 2007
          V2.5.3: Fixed program stop, introduced in V2.5.0, with /NO_CVT keyword.
              MC, Oxford, 3 December 2007
          V2.5.4: Improved color shuffling for final plot.
              MC, Oxford, 30 November 2009
          V2.5.5: Added PIXSIZE keyword. MC, Oxford, 28 April 2010
          V2.5.6: Use IDL intrinsic function DISTANCE_MEASURE for
              automatic pixelSize, when PIXSIZE keyword is not given.
              MC, Oxford, 11 November 2011
          V2.5.7: Included safety termination criterion of Lloyd algorithm
              to prevent loops using /WVT. MC, Oxford, 24 March 2012
          V2.5.8: Update Voronoi tessellation at the exit of bin2d_cvt_equal_mass.
              This is only done when using /WVT, as DIFF may not be zero at the
              last iteration. MC, La Palma, 15 May 2012
          V2.6.0: Included new SN_FUNCTION to illustrate the fact that the user can
              define his own function to estimate the S/N of a bin if needed.
              MC, London, 19 March 2014
          V3.0.0: Translated from IDL into Python and tested against the original.
              MC, London, 19 March 2014
          V3.0.1: Support both Python 2.6/2.7 and Python 3. MC, Oxford, 25 May 2014
          V3.0.2: Avoid potential runtime warning while plotting.
              MC, Oxford, 2 October 2014

    """

    # ----------------------------------------------------------------------------


    def _sn_func(self, signal, noise, index):
        """
        Generic function to calculate the S/N of a bin with spaxels "index".
        The Voronoi binning algorithm does not require this function to have a
        specific form and this generic one can be changed by the user if needed.

        The S/N returned by this function does not need to be an analytic function
        of S and N. There is no need for this function to return the actual S/N.
        Instead this function could return any quantity the user needs to equalize.

        For example _sn_func could be a procedure which uses ppxf to measure the
        velocity dispersion from the coadded spectrum of spaxels "index" and
        returns the relative error in the dispersion.
        Of course an analytic approximation of S/N speeds up the calculation.

        """
        return np.sum(signal[index]) / np.sqrt(np.sum(noise[index] ** 2))

    # ----------------------------------------------------------------------


    def _weighted_centroid(self, x, y, density):
        """
        Computes weighted centroid of one bin.
        Equation (4) of Cappellari & Copin (2003)

        """
        mass = np.sum(density)
        xBar = np.sum(x * density) / mass
        yBar = np.sum(y * density) / mass

        return xBar, yBar

    # ----------------------------------------------------------------------

    def _roundness(self, x, y, pixelSize):
        """
        Implements equation (5) of Cappellari & Copin (2003)

        """
        n = x.size
        equivalentRadius = np.sqrt(n / np.pi) * pixelSize
        xBar, yBar = np.mean(x), np.mean(y)  # Geometric centroid here!
        maxDistance = np.sqrt(np.max((x - xBar) ** 2 + (y - yBar) ** 2))
        roundness = maxDistance / equivalentRadius - 1.

        return roundness

    # ----------------------------------------------------------------------

    def _accretion(self, x, y, signal, noise, targetSN, pixelSize, quiet):
        """
        Implements steps (i)-(v) in section 5.1 of Cappellari & Copin (2003)

        """
        n = x.size

        # will contain the bin number of each given pixel

        classe = np.zeros(n, dtype=int)

        # will contain 1 if the bin has been accepted as good

        good = np.zeros(n, dtype=bool)

        # For each point, find the distance to all other
        # points and select the minimum.
        # This is a robust but slow way of determining
        # the pixel size of unbinned data.
        
        if pixelSize is None:
            pixelSize = np.min(distance.pdist(np.column_stack([x, y])))

        # Start from the pixel with highest S/N

        currentBin = np.argmax(signal / noise)

        SN = signal[currentBin] / noise[currentBin]

        # Rough estimate of the expected final bin number.
        # This value is only used to give an idea of the expected
        # remaining computation time when binning very big dataset.
        #
        w = signal / noise < targetSN
        maxnum = int(np.sum((signal[w] / noise[w]) ** 2)
        / targetSN ** 2 + np.sum(~w))

        # The first bin will be assigned CLASS = 1
        # With N pixels there will be at most N bins
        #
        for ind in range(1, n + 1):

            if not quiet:
                print(ind, ' / ', maxnum)

            classe[currentBin] = ind  # Here currentBin is still made of one pixel
            xBar, yBar = x[currentBin], y[currentBin]    # Centroid of one pixels

            while True:

                if np.all(classe):
                    break  # Stops if all pixels are binned

                # Find the unbinned pixel closest to the centroid of the current bin
                #
                unBinned = np.where(classe == 0)[0]
                k = np.argmin((x[unBinned] - xBar)**2 + (y[unBinned] - yBar)**2)

                # (1) Find the distance from the closest pixel to the current bin
                #
                minDist = np.min((x[currentBin] - x[unBinned[k]])**2 + (y[currentBin] - y[unBinned[k]])**2)

                # (2) Estimate the `roundness' of the POSSIBLE new bin
                #
                nextBin = np.append(currentBin, unBinned[k])
                roundness = self._roundness(x[nextBin], y[nextBin], pixelSize)

                # (3) Compute the S/N one would obtain by adding
                # the CANDIDATE pixel to the current bin
                #
                SNOld = SN
                SN = self._sn_func(signal, noise, nextBin)

                # Test whether (1) the CANDIDATE pixel is connected to the
                # current bin, (2) whether the POSSIBLE new bin is round enough
                # and (3) whether the resulting S/N would get closer to targetSN
                #
                if (np.sqrt(minDist) > 1.2*pixelSize or roundness > 0.3
                    or abs(SN - targetSN) > abs(SNOld - targetSN)):
                    if SNOld > 0.8*targetSN:
                        good[currentBin] = 1
                    break

                # If all the above 3 tests are negative then accept the CANDIDATE
                # pixel, add it to the current bin, and continue accreting pixels
                #
                classe[unBinned[k]] = ind
                currentBin = nextBin

                # Update the centroid of the current bin
                #
                xBar, yBar = np.mean(x[currentBin]), np.mean(y[currentBin])

            # Get the centroid of all the binned pixels
            #
            binned = classe > 0
            if np.all(binned):
                break  # Stop if all pixels are binned
            xBar, yBar = np.mean(x[binned]), np.mean(y[binned])

            # Find the closest unbinned pixel to the centroid of all
            # the binned pixels, and start a new bin from that pixel.
            #
            unBinned = np.where(classe == 0)[0]
            k = np.argmin((x[unBinned] - xBar)**2 + (y[unBinned] - yBar)**2)
            currentBin = unBinned[k]    # The bin is initially made of one pixel
            SN = signal[currentBin]/noise[currentBin]

        classe *= good  # Set to zero all bins that did not reach the target S/N

        return classe, pixelSize

    #----------------------------------------------------------------------------

    def _reassign_bad_bins(self, classe, x, y):
        """
        Implements steps (vi)-(vii) in section 5.1 of Cappellari & Copin (2003)

        """
        # Find the centroid of all successful bins.
        # CLASS = 0 are unbinned pixels which are excluded.
        #
        good = np.unique(classe[classe > 0])
        xnode = ndimage.mean(x, labels=classe, index=good)
        ynode = ndimage.mean(y, labels=classe, index=good)

        # Reassign pixels of bins with S/N < targetSN
        # to the closest centroid of a good bin
        #
        bad = classe == 0
        index = np.argmin((x[bad, None] - xnode)**2 + (y[bad, None] - ynode)**2, axis=1)
        classe[bad] = good[index]

        # Recompute all centroids of the reassigned bins.
        # These will be used as starting points for the CVT.
        #
        good = np.unique(classe)
        xnode = ndimage.mean(x, labels=classe, index=good)
        ynode = ndimage.mean(y, labels=classe, index=good)

        return xnode, ynode

    #----------------------------------------------------------------------------

    def _cvt_equal_mass(self, x, y, signal, noise, xnode, ynode, quiet, wvt):
        """
        Implements the modified Lloyd algorithm
        in section 4.1 of Cappellari & Copin (2003).

        NB: When the keyword WVT is set this routine includes
        the modification proposed by Diehl & Statler (2006).

        """
        if wvt:
            dens = np.ones_like(signal)
        else:
            dens = (signal/noise)**2  # See beginning of section 4.1 of CC03
        scale = np.ones_like(xnode)   # Start with the same scale length for all bins

        for it in range(1, xnode.size):  # Do at most xnode.size iterations

            xnodeOld, ynodeOld = xnode.copy(), ynode.copy()

            # Computes (Weighted) Voronoi Tessellation of the pixels grid
            #
            classe = np.argmin(((x[:, None] - xnode)**2 + (y[:, None] - ynode)**2)/scale**2, axis=1)

            # Computes centroids of the bins, weighted by dens**2.
            # Exponent 2 on the density produces equal-mass Voronoi bins.
            # The geometric centroids are computed if WVT keyword is set.
            #
            good = np.unique(classe)
            for k in good:
                index = classe == k   # Find subscripts of pixels in bin k.
                xnode[k], ynode[k] = self._weighted_centroid(x[index], y[index], dens[index]**2)
                if wvt:
                    sn = self._sn_func(signal, noise, index)
                    scale[k] = np.sqrt(index.sum()/sn)  # Eq. (4) of Diehl & Statler (2006)

            diff = np.sum((xnode - xnodeOld)**2 + (ynode - ynodeOld)**2)

            if not quiet:
                print('Iter: %4i  Diff: %.4g' % (it, diff))

            if diff == 0:
                break

        # If coordinates have changed, re-compute (Weighted) Voronoi Tessellation of the pixels grid
        #
        if diff > 0:
            classe = np.argmin(((x[:, None] - xnode)**2 + (y[:, None] - ynode)**2)/scale**2, axis=1)
            good = np.unique(classe)  # Check for zero-size Voronoi bins

        # Only return the generators and scales of the nonzero Voronoi bins

        return xnode[good], ynode[good], scale[good], it

    #-----------------------------------------------------------------------

    def _compute_useful_bin_quantities(self, x, y, signal, noise, xnode, ynode, scale):
        """
        Recomputes (Weighted) Voronoi Tessellation of the pixels grid to make sure
        that the class number corresponds to the proper Voronoi generator.
        This is done to take into account possible zero-size Voronoi bins
        in output from the previous CVT (or WVT).

        """
        # classe will contain the bin number of each given pixel
        #
        classe = np.argmin(((x[:, None] - xnode)**2 + (y[:, None] - ynode)**2)/scale**2, axis=1)

        # At the end of the computation evaluate the bin luminosity-weighted
        # centroids (xbar, ybar) and the corresponding final S/N of each bin.
        #
        xbar = np.empty_like(xnode)
        ybar = np.empty_like(xnode)
        sn = np.empty_like(xnode)
        area = np.empty_like(xnode)
        good = np.unique(classe)
        for k in good:
            index = classe == k   # Find subscripts of pixels in bin k.
            xbar[k], ybar[k] = self._weighted_centroid(x[index], y[index], signal[index])
            sn[k] = self._sn_func(signal, noise, index)
            area[k] = index.sum()

        return classe, xbar, ybar, sn, area

    #-----------------------------------------------------------------------

    def _display_pixels(self, x, y, counts, pixelSize):
        """
        Display pixels at coordinates (x, y) coloured with "counts".
        This routine is fast but not fully general as it assumes the spaxels
        are on a regular grid. This needs not be the case for Voronoi binning.

        """
        xmin, xmax = np.min(x), np.max(x)
        ymin, ymax = np.min(y), np.max(y)
        nx = round((xmax - xmin)/pixelSize) + 1
        ny = round((ymax - ymin)/pixelSize) + 1
        img = np.full((nx, ny), np.nan)  # use nan for missing data
        j = np.round((x - xmin)/pixelSize).astype(int)
        k = np.round((y - ymin)/pixelSize).astype(int)
        img[j, k] = counts

        plt.imshow(np.rot90(img), interpolation='none', cmap='prism',
                   extent=[xmin - pixelSize/2, xmax + pixelSize/2,
                           ymin - pixelSize/2, ymax + pixelSize/2])

    #----------------------------------------------------------------------

    def voronoi_2d_binning(self, x, y, signal, noise, targetSN, cvt=True,
                             pixelsize=None, plot=True, quiet=True, wvt=True):
        """
        PURPOSE:
              Perform adaptive spatial binning of Integral-Field Spectroscopic
              (IFS) data to reach a chosen constant signal-to-noise ratio per bin.
              This method is required for the proper analysis of IFS
              observations, but can also be used for standard photometric
              imagery or any other two-dimensional data.
              This program precisely implements the algorithm described in
              section 5.1 of the reference below.

        EXPLANATION:
              Further information on VORONOI_2D_BINNING algorithm can be found in
              Cappellari M., Copin Y., 2003, MNRAS, 342, 345

        CALLING SEQUENCE:

            binNum, xBin, yBin, xBar, yBar, sn, nPixels, scale = \
                voronoi_2d_binning(x, y, signal, noise, targetSN,
                                   plot=True, quiet=False, wvt=False,
                                   cvt=True, pixelsize=None)

        """
        # This is the main program that has to be called from external programs.
        # It simply calls in sequence the different steps of the algorithms
        # and optionally plots the results at the end of the calculation.

        if not (x.size == y.size == signal.size == noise.size):
            raise ValueError('Input vectors (x, y, signal, noise) must have the same size')
        if not np.all((noise > 0) & np.isfinite(noise)):
            raise ValueError('NOISE must be a positive vector')

        # Perform basic tests to catch common input errors
        #
        if np.sum(signal)/np.sqrt(np.sum(noise**2)) < targetSN:
            raise ValueError("""Not enough S/N in the whole set of pixels.
                Many pixels may have noise but virtually no signal.
                They should not be included in the set to bin,
                or the pixels should be optimally weighted.
                See Cappellari & Copin (2003, Sec.2.1) and README file.""")
        if np.min(signal/noise) > targetSN:
            raise ValueError('All pixels have enough S/N and binning is not needed')

        # Prevent division by zero for pixels with signal=0 and
        # noise=sqrt(signal)=0 as can happen with X-ray data
        #
        noise = noise.clip(np.min(noise[noise > 0])*1e-9)

        print('Bin-accretion...')
        classe, pixelsize = self._accretion(x, y, signal, noise, targetSN, pixelsize, quiet)
        print(np.max(classe), ' initial bins.')
        print('Reassign bad bins...')
        xNode, yNode = self._reassign_bad_bins(classe, x, y)
        print(xNode.size, ' good bins.')
        if cvt:
            print('Modified Lloyd algorithm...')
            xNode, yNode, scale, it = self._cvt_equal_mass(x, y, signal, noise, xNode, yNode, quiet, wvt)
            print(it-1, ' iterations.')
        else:
            scale = 1.
        classe, xBar, yBar, sn, area = self._compute_useful_bin_quantities(x, y, signal, noise, xNode, yNode, scale)
        w = area == 1
        print('Unbinned pixels: ', np.sum(w), ' / ', x.size)
        print('Fractional S/N scatter (%):', np.std(sn[~w] - targetSN, ddof=1)/targetSN*100)

        if plot:
            plt.clf()
            plt.subplot(211)
            rnd = np.argsort(np.random.random(xNode.size))  # Randomize bin colors
            self._display_pixels(x, y, rnd[classe], pixelsize)
            plt.plot(xNode, yNode, '+w', scalex=False, scaley=False) # do not rescale after imshow()
            plt.xlabel('R (arcsec)')
            plt.ylabel('R (arcsec)')
            plt.title('Map of Voronoi bins')

            plt.subplot(212)
            rad = np.sqrt(xBar**2 + yBar**2)  # Use centroids, NOT generators
            plt.plot(rad[~w], sn[~w], 'or', label='Voronoi bins')
            plt.xlabel('R (arcsec)')
            plt.ylabel('Bin S/N')
            plt.axis([np.min(rad), np.max(rad), 0, np.max(sn)])  # x0, x1, y0, y1
            if np.sum(w) > 0:
                plt.plot(rad[w], sn[w], 'xb', label='single spaxels')
            plt.axhline(targetSN)
            plt.legend()
            plt.show()  # allow plot to appear in certain cases

        return classe, xNode, yNode, xBar, yBar, sn, area, scale

    #----------------------------------------------------------------------------
    def apply_voronoi_binning(self, infile, out_dir, target_sn):

        """
        Usage example for the procedure VORONOI_2D_BINNING.

        It is assumed below that the file voronoi_2d_binning_example.txt
        resides in the current directory. Here columns 1-4 of the text file
        contain respectively the x, y coordinates of each SAURON lens
        and the corresponding Signal and Noise.

        """

        x, y, signal, noise = np.loadtxt(infile,
                                         unpack=1,
                                         skiprows=3)

        # Perform the actual computation. The vectors
        # (binNum, xNode, yNode, xBar, yBar, sn, nPixels, scale)
        # are all generated in *output*
        #
        binNum, xNode, yNode, xBar, yBar, sn, nPixels, scale = self.voronoi_2d_binning(
            x, y, signal, noise, target_sn, plot=1, quiet=0)

        # Save to a text file the initial coordinates of each pixel together
        # with the corresponding bin number computed by this procedure.
        # binNum uniquely specifies the bins and for this reason it is the only
        # number required for any subsequent calculation on the bins.
        #
        np.savetxt(out_dir + '/voronoi_2d_binning_output.txt', np.column_stack([x, y, binNum]),
                   fmt=b'%10.6f %10.6f %8i')

    # BACK TO OWEN CODE
    def vor_output_fitting(self,
                           sci_dir,
                           vor_output,
                           incube,
                           std_cube,
                           sky_cube,
                           centre_x,
                           centre_y,
                           z,
                           stack='sum',
                           line='oiii'):

        """
        Def: Take the output from the above voronoi binning methods
        which is a txt file containing the pixel coordinates and the
        allocated bins, extract each pixel from the incube by assigned bin,
        stack these together using the chosen stack method and then
        fit a gaussian to the appropriate emission line which should now
        be stacked to the apppropriate s/n level. Finally return to the output
        file and add a new column giving the relevant velocity value for each
        of the bin numbers. Each pixel should have a velocity value associated
        with it - can then choose later to ignore some of these which clearly
        don't make any sense. Note the vor_output file shouldn't have any
        column names so that the bin_values are assigned to 'col3'

        Input:
                vor_output - txt file produced from the voronoi binning alg.
                                it contains the pixel coordinates and the
                                bin allocated to each coordinate.
                incube - The combined cube to which this analysis applies.
                z - the redshift of the object in question, used to determine
                    the position of different emission lines.
                stack - method used to stack the spaxels together in each bin
                line - the emission line under scrutiny for this test.
        Output:
                output.txt - file containing the same information as the
                            vor_output txt file but with an additional column
                            to show the velocity value associated with each
                            pixel.
        """
        # sanity checks

        if not(line != 'oiii' or line != 'oii' or line != 'hb'):

            raise ValueError('Please ensure that you have'
                             + ' chosen an appropriate emission line')

        if not(stack != 'sum' or stack != 'median' or stack != 'average'):

            raise ValueError('Please ensure that you have'
                             + ' chosen an appropriate stacking method')
        # first load the output file

        table = ascii.read(vor_output)

        # assign the bin_numbers

        bin_arr = table['col3']

        # now have the bin list in array form. Look for unique bin entries
        # and start a dictionary containing a unique bin key and a tuple of
        # spaxel coordinates which correpsond to this.

        bin_dict = dict()

        for entry in table:

            if entry[2] in bin_dict:

                bin_dict[entry[2]].append([entry[0], entry[1]])

            else:

                bin_dict[entry[2]] = [[entry[0], entry[1]]]

        # the bin dictionary is now populated with the pixel coordinates
        # and has the bin numbers as keys. Time for first external function
        # vor_pixel_stack which will create stacks of these pixels in each bin
        # for the chosen cube

        wave_array, stack_dict = self.vor_pixel_stack(incube,
                                                      bin_dict,
                                                      stack)

        # now need to fit the spectrum around the chosen emission line
        # for each of the spectra in stack_dict. use vor_gauss_fit.
        # first optimally extract the spectrum using
        # no binning for the wavelength axis (could do this)

        spectrum, wl = self.galExtract(sci_dir,
                                       std_cube,
                                       incube,
                                       sky_cube,
                                       centre_y,
                                       centre_x,
                                       1)

        # get the dictionary of initial gaussian params for chosen line

        if line == 'oiii':

            params, blank \
                = self.fit_lines_K(spectrum,
                                   wl,
                                   redshift=z)

            central_wl = 0.500824 * (1. + z)

            input_params = dict()

            input_params['centre'] = params['centre_oiii']
            input_params['sigma'] = params['sigma_oiii']
            input_params['amplitude'] = params['amplitude_oiii']

        elif line == 'hb':

            blank, params \
                = self.fit_lines_K(spectrum,
                                   wl,
                                   redshift=z)

            central_wl = 0.486268 * (1. + z)

            input_params = dict()

            input_params['centre'] = params['centre_hb']
            input_params['sigma'] = params['sigma_hb']
            input_params['amplitude'] = params['amplitude_hb']

        elif line == 'oii':

            blank_1, blank_2, params \
                = self.fit_lines_HK(spectrum,
                                    wl,
                                    redshift=z)

            central_wl = 0.3729875 * (1. + z)

            input_params = dict()

            input_params['centre'] = params['centre_oii']
            input_params['sigma'] = params['sigma_oii']
            input_params['amplitude'] = params['amplitude_oii']

        # now loop around the stack dictionary entries and prepare the
        # fit_wl and fit_flux for input into the vor_gauss_fit method
        # initialise the velocity dictionary

        vel_dict = dict()
        sig_dict = dict()
        flux_dict = dict()

        for entry in stack_dict:

            spec = stack_dict[entry]

            # print spec

            # find the appropriate line index and perform the fit
            line_idx = np.argmin(np.abs(wave_array - central_wl))

            fit_wl = wave_array[line_idx - 8: line_idx + 8]
            fit_flux = spec[line_idx - 8: line_idx + 8]

            best_values = self.vor_gauss_fit(fit_wl,
                                             fit_flux,
                                             input_params)

            # assuming that the redshift measured in qfits is the
            # correct one - subtract the fitted centre and convert
            # to kms-1

            c = 2.99792458E5

            vel = c * ((best_values['center']
                        - central_wl) / central_wl)

            sig = c * ((best_values['sigma']) / central_wl)

            vel_dict[entry] = vel
            sig_dict[entry] = sig
            flux_dict[entry] = best_values['amplitude']

        plt.close('all')

        # have the velocity values for each bin now. Need to return and
        # assign to every pixel the correct velocity value and sigma value

        vel_list = []
        sig_list = []
        flux_list = []

        for entry in bin_arr:

            vel_list.append(vel_dict[entry])
            sig_list.append(sig_dict[entry])
            flux_list.append(flux_dict[entry])

        vel_list = np.array(vel_list)
        sig_list = np.array(sig_list)
        flux_list = np.array(flux_list)

        # and reshape to the 2D format

        cube_data_x = cubeOps(incube).data.shape[1]
        cube_data_y = cubeOps(incube).data.shape[2]

        vel_2d = vel_list.reshape((cube_data_x, cube_data_y))
        sig_2d = sig_list.reshape((cube_data_x, cube_data_y))
        flux_2d = flux_list.reshape((cube_data_x, cube_data_y))

        # plot the results

        vel_fig, vel_ax = plt.subplots(figsize=(18, 6), nrows=1, ncols=3)

        vel_ax[0].minorticks_on()
        vel_ax[1].minorticks_on()
        vel_ax[2].minorticks_on()

        # sometimes this throws a TypeError if hardly any data points
        try:

            vel_min, vel_max = np.nanpercentile(vel_2d, [10.0, 90.0])
            sig_min, sig_max = np.nanpercentile(sig_2d, [10.0, 90.0])
            flux_min, flux_max = np.nanpercentile(flux_2d, [10.0, 90.0])

        except TypeError:

            # origin of the error is lack of good S/N data
            # can set the max and min at whatever
            vel_min, vel_max = [-100, 100]
            sig_min, sig_max = [0, 100]
            flux_min, flux_max = [0, 5E-19]

        im_vel = vel_ax[0].imshow(vel_2d, aspect='auto',
                                  vmin=vel_min,
                                  vmax=vel_max,
                                  interpolation='nearest',
                                  cmap=plt.get_cmap('jet'))

        vel_ax[0].set_title(line + ' velocity')

        # add colourbar to each plot
        divider_vel = make_axes_locatable(vel_ax[0])
        cax_vel = divider_vel.append_axes('right', size='10%', pad=0.05)
        plt.colorbar(im_vel, cax=cax_vel)

        im_sig = vel_ax[1].imshow(sig_2d, aspect='auto',
                                  vmin=sig_min,
                                  vmax=sig_max,
                                  interpolation='nearest',
                                  cmap=plt.get_cmap('jet'))

        vel_ax[1].set_title(line + ' dispersion')

        # add colourbar to each plot
        divider_sig = make_axes_locatable(vel_ax[1])
        cax_sig = divider_sig.append_axes('right', size='10%', pad=0.05)
        plt.colorbar(im_sig, cax=cax_sig)

        im_vel = vel_ax[2].imshow(flux_2d, aspect='auto',
                                  vmin=flux_min,
                                  vmax=flux_max,
                                  interpolation='nearest',
                                  cmap=plt.get_cmap('jet'))

        vel_ax[2].set_title(line + ' flux')

        # add colourbar to each plot
        divider_vel = make_axes_locatable(vel_ax[2])
        cax_vel = divider_vel.append_axes('right', size='10%', pad=0.05)
        plt.colorbar(im_vel, cax=cax_vel)
        # vel_fig.colorbar(im)

        # plt.tight_layout()
        plt.show()
        return flux_2d

    def vor_pixel_stack(self, incube, bin_dict, stack):

        """
        Def: Take the bin_dict - look at each individual key, extract the
        spectra from the given set of pixel coordinates and combine these
        with the given stacking method. End product is a dictionary with the
        same unique keys as before but with one 2048 long spectrum as the
        value. This is a helper function for vor_output_fitting.

        Input:
                incube - the cube to which the pixel coordinates apply
                bin_dict - output from vor_output_fitting
                stack - stacking method set in vor_output_fitting

        Output:
                stack_dict - the resultant dictionary with bin numbers as
                            the keys and stacked spectra as values
                wave_array - the wavelength array corresponding to the flux
        """
        # open the incube

        cube_data = fits.open(incube)[1].data

        # grab the wavelength array from the incube

        wave_array = cubeOps(incube).wave_array

        # intialise new stacking dictionary

        stack_dict = dict()

        # loop around the bin_dict entries

        for entry in bin_dict:

            # initialise temporary stacking list

            temp_list = []

            # loop around the sets of coordinates for each entry

            for coords in bin_dict[entry]:

                temp_list.append(cube_data[:, coords[0], coords[1]])

            # stack the spectra with the chosen stacking method

            if stack == 'sum':

                stacked_flux = np.nansum(temp_list, axis=0)

            elif stack == 'average':

                stacked_flux = np.nanmean(temp_list, axis=0)

            elif stack == 'median':

                stacked_flux = np.nanmedian(temp_list, axis=0)

            # add the stacked flux to the stack_dict

            stack_dict[entry] = stacked_flux

        # return the wave_array and the stack_dict

        return wave_array, stack_dict

    def vor_gauss_fit(self, fit_wl, fit_flux, params):

        """
        Def:
        Performs simple gaussian fit, given initial dictionary of parameters
        and the input wavelength and input flux values.

        Input:
                fit_wl - wavelength of spectrum to fit
                fit_flux - flux of spectrum to fitsWavelength
                params - dictionary containing the keys
                        centre, sigma and amplitude which are initial guesses
                        at the gaussian parameters determined from an
                        integrated fit to the spectrum

        Output:
                fit_params - dictionary containing the best fit parameters
                            for each of the spectra
        """

        # construct gaussian model using lmfit

        gmod = GaussianModel()

        # set the initial parameter values

        pars = gmod.make_params()

        pars['center'].set(value=params['centre'],
                           min=params['centre'] - 0.0015,
                           max=params['centre'] + 0.0015)

        pars['sigma'].set(value=params['sigma'],
                          min=params['sigma'] - (0.5 * params['sigma']),
                          max=params['sigma'] + (0.5 * params['sigma']))

        pars['amplitude'].set(value=params['amplitude'])

        # perform the fit
        out = gmod.fit(fit_flux, pars, x=fit_wl)

        # print the fit report
        # print out.fit_report()

        # plot to make sure things are working
        fig, ax = plt.subplots(figsize=(14, 6))
        ax.plot(fit_wl, fit_flux, color='blue')
        ax.plot(fit_wl, out.best_fit, color='red')
        # plt.show()
        plt.close('all')

        return out.best_values

    def hb_metallicity(self, oiii_flux, hb_flux):

        """
        Def:
        take two 2d flux arrays, one of which is oiii and one of which is
        hb and then compute the 2d metallicity distribution from that. Uses
        the maiolino 2008 calibration.

        Input:
                oiii_flux - 2d oiii flux distribution
                hb_flux - 2d hb flux distribution

        Output: 
                met_2d - 2d metallicity plot
        """
        # take the ratio of the two 2d arrays
        overall_met = oiii_flux / hb_flux

        x_shape = overall_met.shape[0]
        y_shape = overall_met.shape[1]

        hb_met_array = np.empty(shape=(x_shape, y_shape))

        # initialise the coefficients, given in Maiolino 2008
        c_0_hb = 0.1549
        c_1_hb = -1.5031
        c_2_hb = -0.9790
        c_3_hb = -0.0297

        for i, xpix in enumerate(np.arange(0, x_shape, 1)):

            for j, ypix in enumerate(np.arange(0, y_shape, 1)):
                # print 'This is the number: %s' % overall_met[i, j]

                # if the number is nan, leave it as nan

                if np.isnan(overall_met[i, j]) \
                   or np.isinf(overall_met[i, j]) \
                   or (overall_met[i, j]) < 0:

                    hb_met_array[i, j] = np.nan

                # else subtract the log10(number) from
                # c_0_Hb and set up the polynomial from poly1D

                else:

                    c_0_hb_new = c_0_hb - np.log10(overall_met[i, j])

                    p = poly1d([c_3_hb, c_2_hb, c_1_hb, c_0_hb_new])
                    # print p.r
                    # the roots of the polynomial are given in units
                    # of metallicity relative to solar. add 8.69
                    # met_value = p.r[0] + 8.69
                    # if the root has an imaginary component, just take
                    # the real part
                    hb_met_array[i, j] = p.r[2].real + 8.69

        # plot the results

        fig, ax = plt.subplots(1, figsize=(14, 14))

        im = ax.imshow(hb_met_array,
                       aspect='auto',
                       vmin=7.5,
                       vmax=9.0,
                       interpolation='bicubic',
                       cmap=plt.get_cmap('jet'))

        ax.set_title('log([OIII] / Hb)')

        fig.colorbar(im)
        plt.show()
        plt.close('all')

    def voronoi_binning_by_line(self,
                                line,
                                incube,
                                redshift,
                                target_sn,
                                out_dir):

        """
        Def:
        Find the signal and noise columns for each pixel, as well as the pixel
        values themselves and apply the voronoi binning method.
        Write the output from the voronoi binning to an output file.
        This can then be fed directly into the above voronoi fitting method
        """
        # create cube object

        cube = cubeOps(incube)

        # use the cubeclass method to make the signal and noise maps

        signal_2d, noise_2d = cube.make_sn_map(line, redshift)

        # ravel to make into 1d vectors

        signal_1d = np.ravel(signal_2d)
        noise_1d = np.ravel(noise_2d)

        # make the coordinate arrays
        xbin_shape = signal_2d.shape[0]
        ybin_shape = signal_2d.shape[1]

        xbin = np.arange(0, xbin_shape, 1)
        ybin = np.arange(0, ybin_shape, 1)

        ybin, xbin = np.meshgrid(ybin, xbin)

        xbin = np.ravel(xbin)
        ybin = np.ravel(ybin)

        # now have everything required to run the voronoi_binning method

        # Perform the actual computation. The vectors
        # (binNum, xNode, yNode, xBar, yBar, sn, nPixels, scale)
        # are all generated in *output*
        #
        binNum, xNode, yNode, \
            xBar, yBar, sn, nPixels, \
            scale = self.voronoi_2d_binning(xbin,
                                            ybin,
                                            signal_1d,
                                            noise_1d,
                                            target_sn,
                                            plot=1,
                                            quiet=0)

        # Save to a text file the initial coordinates of each pixel together
        # with the corresponding bin number computed by this procedure.
        # binNum uniquely specifies the bins and for this reason it is the only
        # number required for any subsequent calculation on the bins.
        #
        np.savetxt(out_dir + '/voronoi_2d_binning_output.txt',
                   np.column_stack([xbin, ybin, binNum]),
                   fmt=b'%10.6f %10.6f %8i')
