#!/bin/bash
########################################################################################################################
# The goal here is calcualte coaddm5 for nonWFD/nonDD region in baseline2018a for all six bands.
#
# Humna Awan: humna.awan@rutgers.edu
#
########################################################################################################################
# source sims_maf
source /global/common/software/lsst/cori-haswell-gcc/stack/setup_current_sims.sh
setup lsst_sims

dbfile='/global/cscratch1/sd/awan/dbs_wp_unzipped/baseline2018a.db'
outdir='/global/cscratch1/sd/awan/lsst_output/moar_area_output/'

for band in u g r i z y
do
    python /global/homes/a/awan/LSST/lsstRepos/survey_strategy_wp/tools/get_coaddm5_nonwfd_nondd.py \
                    --dbfile=${dbfile} --outdir=${outdir} --band=${band} --nside=256 &
done