########################################################################################################################
# The goal here is calcualte coaddm5 for nonWFD/nonDD region in the specified db for the specified band.
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
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.metricBundles as metricBundles
import lsst.sims.maf.maps as maps
########################################################################################################################
from optparse import OptionParser
parser = OptionParser()
parser.add_option('--nside', dest='nside', type='int',
			help='HEALPix resolution parameter.', default=256)
parser.add_option('--dbfile', dest='dbfile',
                  help='Path to the db file.',
                  default='/global/cscratch1/sd/awan/dbs_wp_unzipped/baseline2018a.db')
parser.add_option('--outdir', dest='main_outdir',
                  help='Path to the folder where all the output should be stored. The folder must already exist.',
                  default='/global/cscratch1/sd/awan/lsst_output/moar_area_output/')
parser.add_option('--band', dest='band',
                  help='Band to consider.',
                  default='i')
# read in the options.
(options, args) = parser.parse_args()
nside = options.nside
dbfile = options.dbfile
main_outdir = options.main_outdir
band = options.band
########################################################################################################################
time0 = time.time()
# connect to the database
opsdb = db.OpsimDatabase(dbfile)
dbname = dbfile.split('/')[-1].split('.db')[0]
# set up the outdir and resultsdb object
dither = 'nodither'
outdir = '%s/coaddm5depth_wdust_%s_nonwfd_nondd_%sband_%s_nside%s'%(main_outdir, dbname, band, dither, nside)
resultsDb = db.ResultsDb(outDir=outdir)
# set up the sql constraint
prop_ids, prop_tags = opsdb.fetchPropInfo()
sqlconstraint = 'proposalId != %s and proposalId != %s and filter=="%s"'%(prop_tags['WFD'][0],
                                                                          prop_tags['DD'][0],
                                                                          band)
# set up the slicer
slicer = slicers.HealpixSlicer(lonCol='fieldRA', latCol='fieldDec',
                               latLonDeg=opsdb.raDecInDeg, nside=nside, useCache=False)
# set up the metric
metric = metrics.ExgalM5(m5Col='fiveSigmaDepth', lsstFilter=band)
# set up the dustmap
dust_map = maps.DustMap(interp=False, nside=nside)
# set up the bundle
bundle = metricBundles.MetricBundle(metric=metric, slicer=slicer,
                                    sqlconstraint=sqlconstraint, mapsList=[dust_map])
# set up the bundlegroup
grp = metricBundles.MetricBundleGroup(bundleDict={0: bundle}, dbObj=opsdb, outDir=outdir,
                                      resultsDb=resultsDb, saveEarly=False)
# run the analysis
grp.runAll()

# save the data
fname = 'coaddm5_wdust_%sband_%s.npz'%(band, dither)
bundle.slicer.writeData('%s/%s'%(outdir, fname),
                        bundle.metricValues,
                        metricName = bundle.metric.name,
                        simDataName = bundle.runName,
                        constraint = bundle.constraint,
                        metadata = bundle.metadata,
                        displayDict = bundle.displayDict,
                        plotDict = bundle.plotDict)
print('Saved %s'%fname)
# remove the results db object
os.remove('%s/resultsDb_sqlite.db'%outdir)
print('## Total time taken: %.2f min'%((time.time()-time0)/60.))