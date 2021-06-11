import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.io import fits
import pandas as pd
from pathlib import Path
import yaml
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("Config_filename")
args = parser.parse_args()

config_filename = Path(args.Config_filename)

with open(config_filename, 'r') as stream:
    try:
        config = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)


catalogue_A = pd.read_csv(f"StarCatalogues/{config['final_star_catalogue_name_stem']}_targets.csv")
catalogue_B = pd.read_csv(f"StarCatalogues/{config['final_star_catalogue_name_stem']}_standards.csv")
catalogue_C = pd.read_csv(f"StarCatalogues/{config['final_star_catalogue_name_stem']}_guides.csv")

image_hdu = fits.open(config['final_image_filename'])
wcs = WCS(image_hdu[0].header)
image = image_hdu[0].data

fig, ax = plt.subplots(subplot_kw=dict(projection=wcs))
ax.imshow(image)
ax.scatter(catalogue_A.RA, catalogue_A.DEC, rasterized=True, transform=ax.get_transform("world"), marker='.', label='Target Stars')
ax.scatter(catalogue_B.RA, catalogue_B.DEC, rasterized=True, transform=ax.get_transform("world"), marker='.', label='Standard stars')
ax.scatter(catalogue_C.RA, catalogue_C.DEC, rasterized=True, transform=ax.get_transform("world"), marker='.', label='Guide stars')
