import pandas as pd
import numpy as np

"""
Make changes to the master star catalogue which we'll use as targets. Sometimes we want to select stars which are bunched up in the centre, sometimes equally spaced etc. This is based on the 'starfield_type' parameter in the config file
"""

config = snakemake.config

star_catalogue = snakemake.input.targets_file
df = pd.read_csv(star_catalogue)

if config['starfield_type'] == 'within_half_radius':

    mask = (df.RA - config['centre_ra'])**2 + (df.DEC - config['centre_dec'])**2 < 0.75**2
    df = df.loc[mask]

elif config['starfield_type'] == 'equally_spaced':
    selected_IDs = np.genfromtxt(snakemake.input.IDs_for_masking, dtype=int)
    mask = df.ID.isin(selected_IDs)
    assert mask.sum() == 19, "We must have 19 stars in our 'equally spaced' tile"

    df = df.loc[mask]

elif config['starfield_type'] == 'entire_catalogue':
    pass


print(f"Our final catalogue has {len(df)} targets")
df.to_csv(snakemake.output.final_catalogue_name, index=False)