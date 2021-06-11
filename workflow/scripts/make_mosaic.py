import numpy as np
import matplotlib.pyplot as plt
from reproject.mosaicking import find_optimal_celestial_wcs
from astropy import units as u
from reproject.mosaicking import reproject_and_coadd
from reproject import reproject_interp
import glob
from astropy.io import fits
from pathlib import Path
from tqdm import tqdm

import yaml
# import argparse

# parser = argparse.ArgumentParser()
# parser.add_argument("Config_filename")
# args = parser.parse_args()

# config_filename = Path(args.Config_filename)
config_filename = snakemake.input[0]

with open(config_filename, 'r') as stream:
    try:
        config = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)


all_files = glob.glob(config['individual_image_glob'])

all_hdus = []
for h in tqdm(all_files):
    all_hdus.append(fits.open(h)[1])

wcs_out, shape_out = find_optimal_celestial_wcs(all_hdus, resolution=config['image_resolution_arcsec'] * u.arcsec)

array, footprint = reproject_and_coadd(all_hdus, wcs_out, shape_out=shape_out, reproject_function=reproject_interp)

# Save a fits file
# Primary extention is the image, the other is the footprint of each individual image
new_hdu = fits.PrimaryHDU(data=array, header=wcs_out.to_header())
new_hdu_footprint = fits.ImageHDU(data=footprint, header=wcs_out.to_header())
hdulist = fits.HDUList([new_hdu, new_hdu_footprint])
#hdulist.writeto(config['final_image_filename'], overwrite=True)
hdulist.writeto(snakemake.output[0], overwrite=True)