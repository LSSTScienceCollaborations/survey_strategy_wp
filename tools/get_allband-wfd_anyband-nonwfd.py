########################################################################################################################
# The goal here is to read in the saved coadd data for all six bands for WFD and nonWFD/nonDD.
# For WFD, this script saves the pixels that have coverage in all six bands, while for nonWFD/nonDD, pixels with
# coverage in any band are saved.
#
# Humna Awan: humna.awan@rutgers.edu
#
########################################################################################################################
import os
import numpy as np
import pandas as pd
import time
import lsst.sims.maf.metricBundles as metricBundles
########################################################################################################################
time0 = time.time()
bands = ['u', 'g', 'r', 'i', 'z', 'y']
dbname = 'baseline2018a'
nside = 256
outdir = '/global/cscratch1/sd/awan/lsst_output/moar_area_output/'
# ------------------------------------------------------------------------------------------------
# path + tags to read in coadd data for the dbs; saved data has masks for each band.
file_yearTag = 'fullSurveyPeriod'
coadd_data_dir = '/global/cscratch1/sd/awan/lsst_output/coadd_output_noDith/'

# Read in the coadd data for the db
data_bundle = {}
for band in bands:
    # set up the folders
    folder = 'coaddM5Analysis_nside%s_withDustExtinction_'%(nside)
    folder += '0pixelRadiusForMasking_%sBand_%s_%s_directory/'%(band, dbname, file_yearTag)
    path = '%s/%s/unmaskedCoaddData/'%(coadd_data_dir, folder)
    # get the filenames
    filenames = [f for f in os.listdir(path) if f.endswith('.npz')]
    print('Reading %s from\n%s/unmaskedCoaddData.\n'%(filenames, folder))

    if len(filenames)>1:
        err = 'Have more than one npz file for %s band'%band
        err += ' for %s data: %s'%(yr_cut, filenames)
        raise ValueError(err)
    else:
        dither = filenames[0].split('%s_'%band)[-1]
        dither = dither.split('.npz')[0]
        mB = metricBundles.createEmptyMetricBundle()
        mB.read('%s/%s'%(path, filenames[0]))
        data_bundle['%s'%(band)]= mB

# find all-band footprint
all_band_pixels = None
for band in data_bundle:
    index = np.where((data_bundle[band].metricValues.mask == False) & \
                     (data_bundle[band].metricValues.data > 0))[0]
    # save the indices
    if all_band_pixels is None:
        # initate the list
        all_band_pixels = index
    else:
        # keep only the overlapping pixels
        all_band_pixels = list(set(all_band_pixels).intersection(index))

fname = '%s_all-band_wfd_pixels_nside%s.csv'%(dbname, nside)
pd.DataFrame({'pixels': all_band_pixels}).to_csv('%s/%s'%(outdir, fname), index=False)
print('Saved %s'%fname)

# ------------------------------------------------------------------------------------------------
# Read in the coadd data for the db
data_bundle = {}
for band in bands:
    # set up the folders
    folder = 'coaddm5depth_wdust_%s_nonwfd_nondd_%sband_nodither_nside%s'%(dbname, band, nside)
    path = '%s/%s/'%(outdir, folder)
    # get the filenames
    filenames = [f for f in os.listdir(path) if f.endswith('.npz')]
    print('Reading %s from\n%s.\n'%(filenames, folder))

    if len(filenames)>1:
        err = 'Have more than one npz file for %s band'%band
        err += ' for %s data: %s'%(yr_cut, filenames)
        raise ValueError(err)
    else:
        dither = filenames[0].split('%s_'%band)[-1]
        dither = dither.split('.npz')[0]
        mB = metricBundles.createEmptyMetricBundle()
        mB.read('%s/%s'%(path, filenames[0]))
        data_bundle['%s'%(band)]= mB

# find any-band footprint
any_band_pixels = None
for band in data_bundle:
    index = np.where((data_bundle[band].metricValues.mask == False) & \
                     (data_bundle[band].metricValues.data > 0))[0]
    # save the indices
    if any_band_pixels is None:
        # initate the list
        any_band_pixels = index
    else:
        # keep only the overlapping pixels
        any_band_pixels = np.hstack([any_band_pixels, index])

fname = '%s_any-band_nonwfd_nondd_pixels_nside%s.csv'%(dbname, nside)
pd.DataFrame({'pixels': np.unique(any_band_pixels)}).to_csv('%s/%s'%(outdir, fname), index=False)
print('Saved %s'%fname)

print('## Total time taken: %.2f min'%((time.time()-time0)/60.))