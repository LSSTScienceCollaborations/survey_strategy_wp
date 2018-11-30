########################################################################################################################
# The goal here is to find and save the mapping between fieldIds and HEALpix pixels for WFD. We use healpy query_disc
# routine to find all the pixels associated with each fieldId.
#
# This code saves two dictionaries: one with keys=fieldIds, values=corresponding pixels; and another with key=pixels,
# and values=correspodning fieldIds.
#
# Humna Awan: humna.awan@rutgers.edu
#
########################################################################################################################
import os
import pickle
import numpy as np
import healpy as hp
import pandas as pd
import time
from astropy import units as u
from astropy.coordinates import SkyCoord

########################################################################################################################
from optparse import OptionParser
parser = OptionParser()
parser.add_option('--nside', dest='nside', type='int',
			help='HEALPix resolution parameter.', default=256)
parser.add_option('--outdir', dest='outdir',
                  help='Path to the folder where all the output should be stored. The folder must already exist.',
                  default='/global/cscratch1/sd/awan/lsst_output/moar_area_output/')

(options, args) = parser.parse_args()
nside = options.nside
outdir = options.outdir

########################################################################################################################
time0 = time.time()

fids_file = '../field_list.csv'   # in this repo
all_fields = pd.read_csv(fids_file)
fov_radius = np.radians(1.75)

# now loop over the HEALpix pixels and find all the correspondence between fieldIds and the pixels
pixels_in_fid, fid_in_pixels = {}, {}
for fid, ra, dec in zip(all_fields.fieldId, all_fields.ra, all_fields.dec):
    c = SkyCoord(ra=ra*u.degree, dec= dec*u.degree)
    pixels = hp.query_disc(nside=nside, vec=c.cartesian.xyz,
                           radius=fov_radius, inclusive=False)
    if fid not in pixels_in_fid:
        pixels_in_fid[fid] = []
    # add this pixel to the fid-list
    pixels_in_fid[fid] += list(pixels)

    for pixel in pixels:
        # check if this pixel is in fid_in_pixels. if not, initiate the list of ids
        if pixel not in fid_in_pixels:
            fid_in_pixels[pixel] = []
        # add this fid to the pixel-list
        fid_in_pixels[pixel].append(fid)

for pixel in fid_in_pixels:
    fid_in_pixels[pixel] = np.unique(fid_in_pixels[pixel])
# ------------------------------------------------------------------------------------------------------
# setup to save the filename
# assemble the filename
filename = 'pixels_in_fieldid_nside%s.pickle'%(nside)
# save the data
pickle.dump(pixels_in_fid, open( '%s/%s'%(outdir, filename) , "wb" ) )
print('## Saved %s\n'%filename)

# assemble the filename
filename = 'fieldid_in_pixels_nside%s.pickle'%(nside)
# save the data
pickle.dump(fid_in_pixels, open( '%s/%s'%(outdir, filename) , "wb" ) )
print('## Saved %s\n'%filename)

print('## Total time taken: %.2f min'%((time.time()-time0)/60.))