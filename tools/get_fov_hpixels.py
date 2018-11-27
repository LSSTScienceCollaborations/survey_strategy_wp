########################################################################################################################
# The goal here is to find and save the mapping between fieldIds and HEALpix pixels for WFD. We use pontus_2002 since
# covers the wider WFD footprint.
#
# Humna Awan: humna.awan@rutgers.edu
#
########################################################################################################################
import os
import pickle
import numpy as np
import healpy as hp
import time
import lsst.sims.maf
import lsst.sims.maf.db as db
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.stackers as stackers

########################################################################################################################
from optparse import OptionParser
parser = OptionParser()
parser.add_option('--nside', dest='nside', type='int',
			help='HEALPix resolution parameter.', default=256)
parser.add_option('--dbfile', dest='dbfile',
                  help='Path to the db file.',
                  default='/global/cscratch1/sd/awan/dbs_wp_unzipped/pontus_2002.db')
parser.add_option('--outdir', dest='outdir',
                  help='Path to the folder where all the output should be stored. The folder must already exist.',
                  default='/global/cscratch1/sd/awan/lsst_output/moar_area_output/')
parser.add_option('--nodith',
                  action='store_true', dest='nodith', default=False,
                  help= 'Use the tag to not include translational dithers. \
                          Otherwise random, per night dithers will be implemented for WFD.')
(options, args) = parser.parse_args()
nside = options.nside
dbfile = options.dbfile
outdir = options.outdir
nodith = options.nodith
########################################################################################################################
time0 = time.time()
# connect to the database
opsdb = db.OpsimDatabase(dbfile)
# get some of the metric data
simdata = opsdb.fetchMetricData(colnames=['fieldRA', 'fieldDec', 'night', 'fieldId', 'proposalId'],
                                sqlconstraint=None)
# ------------------------------------------------------------------------------------------------------
# figure out the pointing columns
if nodith:
    pointing_ra_col, pointing_dec_col = 'fieldRA', 'fieldDec'
else:
    # ------------------------------------------------------------------------------------------------------
    # add translational dithers: the stacker adds random dither per night column to simdatadata
    s = stackers.RandomDitherPerNightStacker(degrees=opsdb.raDecInDeg, randomSeed=42)
    simdata = s.run(simdata)

    # revert the non-WFD pointing RA, Dec to undithered ones
    prop_ids, prop_tags = opsdb.fetchPropInfo()
    # find the indices for the non-WFD visits
    ind = np.where(simdata['proposalId']!=prop_tags['WFD'])[0]
    # change the pointings for non-WFD to be undithered
    simdata['randomDitherPerNightRa'][ind] = simdata['fieldRA'][ind]
    simdata['randomDitherPerNightDec'][ind] = simdata['fieldDec'][ind]

    # assign the pointing columns
    pointing_ra_col, pointing_dec_col = 'randomDitherPerNightRa', 'randomDitherPerNightDec'
# ------------------------------------------------------------------------------------------------------
# set up the HEALPix slicer
slicer = slicers.HealpixSlicer(nside=nside, lonCol=pointing_ra_col, latCol=pointing_dec_col,
                               useCache=False, latLonDeg=opsdb.raDecInDeg)
# slice the data
slicer.setupSlicer(simdata)
# ------------------------------------------------------------------------------------------------------
# now loop over the HEALpix pixels and find all the correspondence between fieldIds and the pixels
pixels_in_fov = {}
for pixel in range(hp.nside2npix(nside=nside)):
    # find all the observations in this pixel
    ind_obs_in_pixel = slicer._sliceSimData(pixel)
    # find the fids associated with the all the visits to this pixels
    ids = simdata[ind_obs_in_pixel['idxs']]['fieldId']   # fieldIDs corresponding to pixel
    # now loop over all the unique ids and store the correspondence
    for uniq_id in np.unique(ids):
        if uniq_id not in pixels_in_fov.keys():
            pixels_in_fov[uniq_id] = []
        pixels_in_fov[uniq_id].append(pixel)
# ------------------------------------------------------------------------------------------------------
# add meta data
pixels_in_fov['meta'] = 'sims_maf verison: %s'%(lsst.sims.maf.__version__,)
# ------------------------------------------------------------------------------------------------------
# setup to save the filename
if nodith:
    dith_tag = 'undithered'
else:
    dith_tag = 'randomDitherPerNight'
# assemble the filename
filename = 'pixels_in_fov_%s_nside%s.pickle'%(dith_tag, nside)
# save the data
pickle.dump(pixels_in_fov, open( '%s/%s'%(outdir, filename) , "wb" ) )
print('## Saved data as %s\n'%filename)
print('## Total time taken: %.2f min'%((time.time()-time0)/60.))