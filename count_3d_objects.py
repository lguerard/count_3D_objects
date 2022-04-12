#@ File(label="Select the directory with your cropped images", style="directory") src_dir
#@ String(label="Extension for the images to look for", value="tif") filename_filter
#@ String(label="Channel of interest (if multiple, separate them by commas", value=3) ch_number
#@ Double(label="Intensity threshold for the spots", value=150) intensity_thresh

# ─── IMPORTS ────────────────────────────────────────────────────────────────────


from fiji.plugin.trackmate.detection import LogDetector, LogDetectorFactory
from fiji.plugin.trackmate import Model, Settings, TrackMate, SelectionModel, Logger
from fiji.plugin.trackmate.features import FeatureFilter as FeatureFilter
from fiji.plugin.trackmate.tracking.sparselap import SparseLAPTrackerFactory
from fiji.plugin.trackmate.tracking import LAPUtils
from fiji.plugin.trackmate.features.track import TrackDurationAnalyzer as TrackDurationAnalyzer
from fiji.plugin.trackmate.providers import SpotAnalyzerProvider
from fiji.plugin.trackmate.providers import EdgeAnalyzerProvider
from fiji.plugin.trackmate.providers import TrackAnalyzerProvider


from net.imglib2.img.display.imagej import ImageJFunctions

import os
import csv
import glob
import sys
from itertools import izip

from ij import IJ, ImagePlus, ImageStack, WindowManager as wm
from ij.plugin.frame import RoiManager
from ij.gui import PointRoi, WaitForUserDialog
from ij.measure import ResultsTable
from ij.process import ImageConverter
from ij.plugin import Duplicator

# Bioformats imports
from loci.plugins import BF
from loci.plugins.in import ImporterOptions

from collections import OrderedDict

# ─── VARIABLES ──────────────────────────────────────────────────────────────────

# ############################# #
# TRACKMATE DETECTION VARIABLES #
# ############################# #

# Radius for the spots /!\ NOT DIAMETER
radius     = 0.35
# Threshold
threshold  = 70.0
# Do Subpixel detection
doSubpixel = True
# Apply median before detection
doMedian   = False


# ─── FUNCTIONS ──────────────────────────────────────────────────────────────────


def checkForFiles(filepath):
    """Check if files are there no matter the extension

    Arguments:
        filepath {string} -- Path and name to check if exists

    Returns:
        bool -- Returns true if exists otherwise returns false
    """
    for filepath_object in glob.glob(filepath):
        if os.path.isfile(filepath_object):
            return True

    return False


def getFileList(directory, filteringString):
    """
    Returns a list containing the file paths in the specified directory
    path. The list is recursive (includes subdirectories) and will only
    include files whose filename contains the specified string.
    """
    files = []
    for (dirpath, dirnames, filenames) in os.walk(directory):
        # if out_dir in dirnames: # Ignore destination directory
            # dirnames.remove(OUT_SUBDIR)
        for f in filenames:
            if filteringString in f:
                files.append(os.path.join(dirpath, f))
    return files

def BFImport(indivFile):
    """
    Import the files using BioFormats.
    """
    options = ImporterOptions()
    options.setId(str(indivFile))
    options.setColorMode(ImporterOptions.COLOR_MODE_COMPOSITE)
    imps = BF.openImagePlus(options)
    return imps


def runTM(imp, rad, thresh, subpix, med, mean_int):

    #----------------------------
    # Create the model object now
    #----------------------------

    # Some of the parameters we configure below need to have
    # a reference to the model at creation. So we create an
    # empty model now.

    model = Model()

    # Send all messages to ImageJ log window.
    model.setLogger(Logger.IJ_LOGGER)


    #------------------------
    # Prepare settings object
    #------------------------

    settings = Settings(imp)

    spotAnalyzerProvider = SpotAnalyzerProvider(1)
    for key in spotAnalyzerProvider.getKeys():
        # print( key )
        settings.addSpotAnalyzerFactory( spotAnalyzerProvider.getFactory( key ) )

    # edgeAnalyzerProvider = EdgeAnalyzerProvider()
    # for  key in edgeAnalyzerProvider.getKeys():
    #     # print( key )
    #     settings.addEdgeAnalyzer( edgeAnalyzerProvider.getFactory( key ) )

    # trackAnalyzerProvider = TrackAnalyzerProvider()
    # for key in trackAnalyzerProvider.getKeys():
    #     # print( key )
    #     settings.addTrackAnalyzer( trackAnalyzerProvider.getFactory( key ) )

    # Configure detector - We use the Strings for the keys
    settings.detectorFactory = LogDetectorFactory()
    settings.detectorSettings = {
        'DO_SUBPIXEL_LOCALIZATION' : subpix,
        'RADIUS' : rad,
        'TARGET_CHANNEL' : 1,
        'THRESHOLD' : thresh,
        'DO_MEDIAN_FILTERING' : med,
    }

    # Add the filter on mean intensity
    # Here 'false' takes everything BELOW the mean_int value
    if (mean_int != 0):
        filter1 = FeatureFilter('MEAN_INTENSITY', mean_int, False)
        settings.addSpotFilter(filter1)


    # Configure tracker
    settings.trackerFactory = SparseLAPTrackerFactory()
    settings.trackerSettings = LAPUtils.getDefaultLAPSettingsMap()
    settings.trackerSettings['LINKING_MAX_DISTANCE'] = 10.0
    settings.trackerSettings['GAP_CLOSING_MAX_DISTANCE']=10.0
    settings.trackerSettings['MAX_FRAME_GAP']= 3


    settings.initialSpotFilterValue = 1

    # print(str(settings))


    #-------------------
    # Instantiate plugin
    #-------------------

    trackmate = TrackMate(model, settings)

    #--------
    # Process
    #--------

    ok = trackmate.checkInput()
    if not ok:
        sys.exit(str(trackmate.getErrorMessage()))

    ok = trackmate.process()
    if not ok:
        if "[SparseLAPTracker] The spot collection is empty." in str(trackmate.getErrorMessage()):
            return 0
        else:
            sys.exit(str(trackmate.getErrorMessage()))

    IJ.log("Found " + str(model.getSpots().getNSpots(True)) + " spots")

    return model.getSpots().getNSpots(True)



# ─── MAIN CODE ──────────────────────────────────────────────────────────────────

# Retrieve list of files
src_dir = str(src_dir)
files = getFileList(src_dir, filename_filter)

ch_number = [int(x) for x in ch_number.split(",")]
channels_array_text = ['C' + str(s) for s in ch_number]
channels_text = "_".join(channels_array_text)

results_dict = {}

# If the list of files is not empty
if files:
    # Lists for the result file
    results_dict['1. Filename'] = []
    results_dict['2. ROI Area'] = []

    outCSV = os.path.join(src_dir, "Results_" + channels_text +".csv")

    # For each file finishing with the filtered string
    for index, file in enumerate(files):
        # Get info for the files
        folder   = os.path.dirname(file)
        basename = os.path.basename(file)
        basename = os.path.splitext(basename)[0]

        # Import the file with BioFormats
        IJ.log("Currently opening " + basename + "...")
        # Actually not using BF since it can't read the ROIs somehow
        # TODO investigate
        # imps = BFImport(str(file))
        imp = IJ.openImage(str(file))
        # imp.show()
        # IJ.run("To ROI Manager", "");

        rm = RoiManager(False)

        roi = imp.getRoi()
        rm.addRoi(roi)

        if rm.getCount() == 0:
                IJ.log("Couldn't load the ROIs for this image. Check what happened.")
                continue

        roi_area = roi.getStatistics().area * imp.getCalibration().pixelWidth * imp.getCalibration().pixelHeight

        results_dict['1. Filename'].append(basename)
        results_dict['2. ROI Area'].append(roi_area)

        # Calculate the number of cells for each channel
        for current_channel in ch_number:

            if index == 0:
                results_dict['3. C' + str(current_channel) + ' cell count'] = []
                results_dict['4. C' + str(current_channel) + ' normalized count'] = []

            IJ.log("Currently working on channel " + str(current_channel))

            IJ.log("Now working on " + str(current_channel) +"...")
            imp_for_tm = Duplicator().run(imp, current_channel,
                                        current_channel, 1, imp.getNSlices(), 1, 1)
            rm.select(imp_for_tm, 0)
            IJ.run(imp_for_tm, "Clear Outside", "stack")
            IJ.run(imp_for_tm, "Select None", "")

            # Get the marker image with the peaks of cells using TrackMate
            cell_count_ch = runTM(
                imp_for_tm, radius, threshold, doSubpixel, doMedian, 0)

            IJ.log("Here : " + str(cell_count_ch))

            imp_for_tm.close()

            results_dict['3. C' + str(current_channel) + ' cell count'].append(cell_count_ch)
            results_dict['4. C' + str(current_channel) + ' normalized count'].append(cell_count_ch/roi_area)

        imp.close()

    results_dict = OrderedDict(sorted(results_dict.items()))

    with open(outCSV, 'wb') as f:
        writer = csv.writer(f)
        writer.writerow(
            results_dict.keys())
        writer.writerows(
            izip(*list(results_dict.values())))

        f.close()



IJ.log('###########################')
IJ.log('Script done, results saved in ' + outCSV)
IJ.log('###########################')
