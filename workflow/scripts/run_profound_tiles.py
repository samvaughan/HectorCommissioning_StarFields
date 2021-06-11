import numpy as np
import matplotlib.pyplot as plt
import subprocess
from pathlib import Path
import yaml
import argparse
import shlex
import glob
from astropy.io import fits
from astropy.wcs import WCS
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument("Config_filename")
args = parser.parse_args()

config_filename = Path(args.Config_filename)

with open(config_filename, 'r') as stream:
    try:
        config = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

all_images = glob.glob(config['individual_image_glob'])
for image in tqdm(all_images):

    hdu = fits.open(image)
    wcs = WCS(hdu[1].header)
    image_centre = wcs.pixel_to_world(*np.array(wcs.array_shape)/2.0)
    ra, dec = image_centre.ra.value, image_centre.dec.value

    input_image = image
    output_mask = Path(f"/Volumes/OtherFiles/Science/Panstarrs_Mosaic/SkyMasks/Small/segmap_StarField_{ra:.2f}_{dec:.2f}.fits")
    output_catalogue = Path("ProfoundCatalogues/profound_catalogue_{ra:.2f}_{dec:.2f}.csv")
    zeromag = config['profound_zero_mag']

    if output_mask.exists():
        continue
    else:
        bash_command = f'./profound_r_script.R {input_image} {output_mask} {output_catalogue} --mag_zero {zeromag}'

    print("Running Profound...")
    subprocess.check_call(shlex.split(bash_command))
    print("Done!")