import numpy as np
from astroquery.gaia import Gaia
from astroquery.utils.tap.core import TapPlus
import pandas as pd
import pandas_tools as P
import argparse
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u
from astropy.table import Table
import time
from pathlib import Path
import yaml
from astroquery.vizier import Vizier
from urllib.error import HTTPError
import scipy.spatial

config = snakemake.config

brightest_magnitude = 10
faintest_magnitude = 16

print(f"Selecting Stars for the {config['FIELD']} field")
print("\t* Running the Vizier query")

field_centre_coordinates = SkyCoord(ra=config['centre_ra'] * u.deg, dec=config['centre_dec'] * u.deg)
vquery = Vizier(columns=['objID', 'RAJ2000', 'DEJ2000',
                             'e_RAJ2000', 'e_DEJ2000',
                             'gmag', 'e_gmag',
                             'rmag', 'e_rmag',
                             'imag', 'e_imag',
                             'zmag', 'e_zmag',
                             'ymag', 'e_ymag'],
                    column_filters={"rmag":
                                    ("<%f" % 18)},
                                    row_limit=1000000)
r = vquery.query_region(field_centre_coordinates, radius=1.5 * u.deg, catalog="II/349/ps1")[0]
df_all_stars = r.to_pandas() 

print(f"\t* Selected {len(df_all_stars)} stars from Panstarrs")
# Rename some columns and add some new ones
df_all_stars.rename(columns=dict(RAJ2000='RA', DEJ2000='DEC', objID='ID', gmag='g_mean_psf_mag', rmag='r_mean_psf_mag', imag='i_mean_psf_mag', zmag='z_mean_psf_mag', ymag='y_mean_psf_mag'), inplace=True)

# remove stars which have close companions in the catalogue or chance alignments
# Do this using a KDTree and "query_ball_tree" https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.query_ball_tree.html#scipy.spatial.KDTree.query_ball_tree

max_sep = 30/60/60 #arcseconds
print(f'\t* Removing stars which have companions within {max_sep * 60 *60} arcseconds...')
X = np.column_stack((df_all_stars.RA.values, df_all_stars.DEC.values))
tree = scipy.spatial.KDTree(X)

# Query_ball_tree finds all points which are within r of another point
# It returns a list of lists- for each point, you get the point itself and any neighbours within r
# So we can then loop through this list and remove stars with more than 1 element in its list

neighbours = tree.query_ball_tree(tree, r=max_sep)

stars_to_keep = []
for neighbour in neighbours:
    if len(neighbour) == 1:
        stars_to_keep.append(neighbour[0])


panstarrs_df = df_all_stars.loc[stars_to_keep]
N_removed_stars = len(df_all_stars) - len(panstarrs_df)
print(f"\t\t...Done! Removed {N_removed_stars} stars")

# Now only keep stars within our bright/faint magnitude limits
panstarrs_df.dropna(subset=['g_mean_psf_mag', 'r_mean_psf_mag', 'i_mean_psf_mag', 'z_mean_psf_mag', 'y_mean_psf_mag'], inplace=True)
magnitude_mask = (panstarrs_df.r_mean_psf_mag > brightest_magnitude) & (panstarrs_df.r_mean_psf_mag < faintest_magnitude)
panstarrs_df = panstarrs_df.loc[magnitude_mask]

print(f"\t* Applying our magnitude limits: we now have {len(panstarrs_df)} stars")

print("\t* Cross matching with Gaia")
Gaia.login(credentials_file=snakemake.input.Gaia_credientials)
user_id = 'svaugh01'
gaia_table_name = f"field_{config['FIELD']}_stars"
full_qualified_table_name = f"user_{user_id}.{gaia_table_name}"

print("\t\tUploading the table to cross match")
Gaia.upload_table(upload_resource=Table.from_pandas(panstarrs_df), table_name=gaia_table_name)

print("\t\tSetting the RA and DEC flags")
Gaia.update_user_table(table_name=full_qualified_table_name,
                       list_of_changes=[["ra","flags","Ra"],
                                        ["dec","flags","Dec"]])

print("\t\tPerforming the cross-match")
xmatch_table_name = f'xmatch_{gaia_table_name}'

Gaia.cross_match(full_qualified_table_name_b=full_qualified_table_name,
                 full_qualified_table_name_a='gaiadr2.gaia_source',
                 results_table_name=xmatch_table_name, radius=1.0)
print("\t\tGetting the results")
crossmatch_query = f"""
SELECT c."separation", a."source_id", a."ref_epoch", a."pm", a."pmra", a."pmra_error", a."pmdec", a."pmdec_error", a."phot_g_mean_mag", a."phot_bp_mean_mag",  a."phot_rp_mean_mag", b."id", b."ra", b."dec", b."e_raj2000", b."e_dej2000", b."g_mean_psf_mag", b."e_gmag", b."r_mean_psf_mag", b."e_rmag", b."i_mean_psf_mag", b."e_imag", b."z_mean_psf_mag", b."e_zmag", b."y_mean_psf_mag", b."e_ymag" FROM gaiadr3.gaia_source AS a, {full_qualified_table_name} AS b, user_{user_id}.{xmatch_table_name} AS c WHERE (c.gaia_source_source_id = a.source_id AND c.{gaia_table_name}_{gaia_table_name}_oid = b.{gaia_table_name}_oid)
"""
# Put in a sleep here as otherwise we tend to get 500 server errors...
time.sleep(5)
job = Gaia.launch_job_async(query=crossmatch_query)
results = job.get_results()

df = results.to_pandas()

df.rename(columns=dict(ra='RA', dec='DEC', id='ID'), inplace=True)

print("\t\tRemoving cross-match duplicates and things without a proper motion")
duplicate_mask = df.duplicated(subset='ID', keep=False)
proper_motion_nan_mask = ~np.isfinite(df.pm)
df = df.loc[~(duplicate_mask | proper_motion_nan_mask)]
print(f"\n\t\tDone! We have {len(df)} stars")

df['priority'] = 1.0
# Add in columns for Re, z, Mstar and surface brightness
df['z'] = 0.0
df['Re'] = 1.0
df['Mstar'] = -99
df['GAL_MU_E_R'] = 19.0

#Account for Proper Motion
# We need a large value of distance in there to avoid an astropy bug- see https://github.com/astropy/astropy/issues/11747
skycoords = SkyCoord(ra=df.RA.values * u.deg, dec=df.DEC.values * u.deg, pm_ra_cosdec=df.pmra.values * u.mas/u.yr, pm_dec=df.pmdec.values * u.mas/u.yr, obstime=Time(df.ref_epoch.values, format='jyear'), distance=20*u.kpc)

print("\t* Applying the proper motion correction")
updated_skycoords = skycoords.apply_space_motion(Time(config['date_for_observations']))

df['RA_pmcorr'] = updated_skycoords.ra.value
df['DEC_pmcorr'] = updated_skycoords.dec.value

print("\t*Selecting Guide Stars")
# Select guide stars
guide_star_mask = (df.g_mean_psf_mag > 13.5) & (df.g_mean_psf_mag < 14.5)

print(f"\t\t... There are {guide_star_mask.sum()} guide stars")

if guide_star_mask.sum() > 1000:
    guide_stars = df.copy().loc[guide_star_mask].sample(1000, random_state=1)
else:
    guide_stars = df.copy().loc[guide_star_mask]
hexabundle_stars = df.copy().loc[(~guide_star_mask) & (df.r_mean_psf_mag > 14)]


# Select standard Stars
def standard_star_priority(df):

    # These values come from the "Stellar_Colour_Cuts/measure_F_star_colours.py" file
    X = np.sqrt(
    ((df.g_mean_psf_mag - df.r_mean_psf_mag) - 0.23)**2 +
    ((df.r_mean_psf_mag - df.i_mean_psf_mag) - 0.05)**2 + 
    ((df.i_mean_psf_mag - df.z_mean_psf_mag) - -0.01)**2 +
    ((df.i_mean_psf_mag - df.z_mean_psf_mag) - -0.01)**2 + 
    ((df.z_mean_psf_mag - df.y_mean_psf_mag) - -0.01)**2 +
    ((df.phot_bp_mean_mag - df.phot_rp_mean_mag) - 0.31)**2)

    return X

print("\tSelecting Standard Stars")
hexabundle_stars['StandardStar_X_Value'] = standard_star_priority(hexabundle_stars)

standard_star_mask = hexabundle_stars['StandardStar_X_Value'] < 2



if np.sum(standard_star_mask) > 300:
    standard_stars = hexabundle_stars.copy().loc[standard_star_mask]
    sorted_standards = standard_stars.sort_values('StandardStar_X_Value', ascending=True)
    standard_stars = sorted_standards.copy().iloc[:300]
else:
    standard_stars = hexabundle_stars.copy().loc[standard_star_mask]

print(f"\t\t... There are {len(standard_stars)} standards")
if (~standard_star_mask).sum() > 5000:
    target_stars = hexabundle_stars.copy().loc[~(hexabundle_stars.ID.isin(standard_stars.ID))].sample(1000, random_state=1)
else:
    target_stars = hexabundle_stars.copy().loc[~(hexabundle_stars.ID.isin(standard_stars.ID))]


print(f"\t\t... There are {len(target_stars)} target stars")
column_renamer = dict(g_mean_psf_mag='g_mag', r_mean_psf_mag='r_mag', i_mean_psf_mag='i_mag', z_mean_psf_mag='z_mag', y_mean_psf_mag='y_mag', phot_g_mean_mag='GAIA_g_mag', phot_bp_mean_mag='GAIA_bp_mag', phot_rp_mean_mag='GAIA_rp_mag', pmra='pmRA', pmdec='pmDEC')

target_stars.rename(columns=column_renamer, inplace=True)
standard_stars.rename(columns=column_renamer, inplace=True)
guide_stars.rename(columns=column_renamer, inplace=True)

target_stars.to_csv(snakemake.output.targets_file, index=False)
standard_stars.to_csv(snakemake.output.standards_file, index=False)
guide_stars.to_csv(snakemake.output.guides_file, index=False)


print("Done!")
