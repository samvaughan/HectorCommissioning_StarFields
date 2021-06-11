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

min_RA = config['centre_ra'] - config['star_catalogue_width_ra']
max_RA = config['centre_ra'] + config['star_catalogue_width_ra']

min_DEC = config['centre_dec'] - config['star_catalogue_width_dec']
max_DEC = config['centre_dec'] + config['star_catalogue_width_dec']


brightest_magnitude = 9
faintest_magnitude = 18

# This is the SQL query we use to select stars from APASS and their match to GAMA. 
panstarrs_sql_query = f"""SELECT TOP 100000
 
matchy.source_id,  
panstarrs.*,
main_gaia.ref_epoch,
main_gaia.ra,  
main_gaia.dec, 
main_gaia.pmra, 
main_gaia.pmdec, 
main_gaia.phot_g_mean_mag, 
main_gaia.phot_bp_mean_mag, 
main_gaia.phot_rp_mean_mag 
 
 
FROM gaiadr2.panstarrs1_best_neighbour AS matchy 
 
JOIN gaiadr2.gaia_source AS main_gaia ON main_gaia.source_id = matchy.source_id 
JOIN gaiadr2.panstarrs1_original_valid AS panstarrs ON panstarrs.obj_id = CAST(matchy.original_ext_source_id AS BIGINT)
 
WHERE (panstarrs.ra > {min_RA}) AND (panstarrs.ra < {max_RA}) AND (panstarrs.dec < {max_DEC}) AND (panstarrs.dec > {min_DEC}) AND (main_gaia.pmra  > -100) AND (main_gaia.pmra < 100) AND (main_gaia.pmdec > -100) AND (main_gaia.pmdec < 100) AND (panstarrs.r_mean_psf_mag > {brightest_magnitude}) AND (panstarrs.r_mean_psf_mag < {faintest_magnitude})
"""

print("Running the GAIA/PANSTARRS query...")
job = Gaia.launch_job_async(panstarrs_sql_query)
r = job.get_results()
print("...Done!")
print(f"Query selected {len(r)} stars")

df = r.to_pandas() 

# Rename some columns and add some new ones
df.rename(columns=dict(ra='RA', dec='DEC', source_id='ID'), inplace=True)
df.dropna(inplace=True)
df['priority'] = 1.0

#Account for Proper Motion
# We need a large value of distance in there to avoid an astropy bug- see https://github.com/astropy/astropy/issues/11747
skycoords = SkyCoord(ra=df.RA.values * u.deg, dec=df.DEC.values * u.deg, pm_ra_cosdec=df.pmra.values * u.mas/u.yr, pm_dec=df.pmdec.values * u.mas/u.yr, obstime=Time(df.ref_epoch.values, format='jyear'), distance=20*u.kpc)
updated_skycoords = skycoords.apply_space_motion(Time(config['date_for_observations']))

df['RA_pmcorr'] = updated_skycoords.ra.value
df['DEC_pmcorr'] = updated_skycoords.dec.value


# Select guide stars
guide_star_mask = df.g_mean_psf_mag < 14.5

guide_stars = df.loc[guide_star_mask]
hexabundle_stars = df.loc[(~guide_star_mask) & (df.r_mean_psf_mag > 14)]


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

hexabundle_stars['StandardStar_X_Value'] = standard_star_priority(hexabundle_stars)

standard_star_mask = hexabundle_stars['StandardStar_X_Value'] < 0.5
standard_stars = hexabundle_stars.loc[standard_star_mask]
target_stars = hexabundle_stars.loc[~standard_star_mask].sample(1000)

column_renamer = dict(g_mean_psf_mag='g_mag', r_mean_psf_mag='r_mag', i_mean_psf_mag='i_mag', z_mean_psf_mag='z_mag', y_mean_psf_mag='y_mag', phot_g_mean_mag='GAIA_g_mag', phot_bp_mean_mag='GAIA_bp_mag', phot_rp_mean_mag='GAIA_rp_mag', pmra='pmRA', pmdec='pmDEC')

target_stars.rename(columns=column_renamer, inplace=True)
standard_stars.rename(columns=column_renamer, inplace=True)
guide_stars.rename(columns=column_renamer, inplace=True)

target_stars.to_csv(f"StarCatalogues/{config['final_star_catalogue_name_stem']}_targets.csv", index=False)
standard_stars.to_csv(f"StarCatalogues/{config['final_star_catalogue_name_stem']}_standards.csv", index=False)
guide_stars.to_csv(f"StarCatalogues/{config['final_star_catalogue_name_stem']}_guides.csv", index=False)





# apass_sql_query = f"""SELECT TOP 100000
 
# matchy.source_id,  
# apass.*,  
# main_gaia.ra,  
# main_gaia.dec, 
# main_gaia.pmra, 
# main_gaia.pmdec, 
# main_gaia.phot_g_mean_mag, 
# main_gaia.phot_bp_mean_mag, 
# main_gaia.phot_rp_mean_mag 
 
 
# FROM gaiadr2.apassdr9_best_neighbour AS matchy 
 
# JOIN gaiadr2.gaia_source AS main_gaia ON main_gaia.source_id = matchy.source_id 
# JOIN external.apassdr9 AS apass ON apass.recno = CAST(matchy.original_ext_source_id AS BIGINT) 
 
# WHERE (apass.raj2000 > 212.5) AND (apass.raj2000 < 222.5) AND (apass.dej2000 < 2) AND (apass.dej2000 > -1) AND (main_gaia.pmra  > -15) AND (main_gaia.pmra < 15) AND (main_gaia.pmdec > -15) AND (main_gaia.pmdec < 15) AND (apass.r_mag > 12) AND (apass.r_mag < 15.5)
# """

# print("Running the GAIA/APASS query...")
# job = Gaia.launch_job_async(apass_sql_query)
# r = job.get_results()
# print("...Done!")
# print(f"Query selected {len(r)} stars")

# df = r.to_pandas() 
# df.rename(columns=dict(ra='RA', dec='DEC', source_id='ID'), inplace=True)
# df.dropna(inplace=True)
# df['priority'] = 1.0
# #df[['Mstar', 'z', 'GAL_MU_0_U', 'GAL_MU_R_at_3Re', 'GAL_MU_E_I', 'Re', 'g_m_i', 'GAL_MU_E_Z', 'GAL_MU_E_G', 'GAL_MU_R_at_2Re', 'GAL_MU_E_R', 'SersicIndex_r', 'Dingoflag', 'GAL_MU_0_Z', 'MassHIpred', 'GAL_MAG_I', 'GAL_MU_0_G', 'Ellipticity_r', 'IFU_diam_2Re', 'GAL_MAG_G', 'WALLABYflag', 'GAL_MU_0_I', 'GAL_MU_0_R', 'GAL_MU_E_U']] = 0.0
# # Split the dataframe into three unique parts. One of these will be 'targets', the others will be 'guide' stars, the third will be 'standard' stars
# A, B, C = np.split(df.sample(n=len(df) - len(df) % 3, random_state=42), 3)

# n_samples = len(A)
# if len(A) > 2000:
#     n_samples = 2000

# A.sample(n_samples).to_csv("catalogues/RA_225_stars_A.csv")
# B.to_csv("catalogues/RA_225_stars_B.csv")
# C.to_csv("catalogues/RA_225_stars_C.csv")






# apass_sql_query = f"""SELECT TOP 100000
 
# matchy.source_id,  
# apass.*,  
# main_gaia.ra,  
# main_gaia.dec, 
# main_gaia.pmra, 
# main_gaia.pmdec, 
# main_gaia.phot_g_mean_mag, 
# main_gaia.phot_bp_mean_mag, 
# main_gaia.phot_rp_mean_mag
 
 
# FROM gaiadr2.apassdr9_best_neighbour AS matchy 
 
# JOIN gaiadr2.gaia_source AS main_gaia ON main_gaia.source_id = matchy.source_id 
# JOIN external.apassdr9 AS apass ON apass.recno = CAST(matchy.original_ext_source_id AS BIGINT) 
 
# WHERE (apass.raj2000 > 339.1) AND (apass.raj2000 < 350.9) AND (apass.dej2000 < -32) AND (apass.dej2000 > -34) AND (main_gaia.pmra  > -15) AND (main_gaia.pmra < 15) AND (main_gaia.pmdec > -15) AND (main_gaia.pmdec < 15) AND (apass.r_mag > 12) AND (apass.r_mag < 15.5)
# """

# print("Running the GAIA/APASS query...")
# job = Gaia.launch_job_async(apass_sql_query)
# r = job.get_results()
# print("...Done!")
# print(f"Query selected {len(r)} stars")

# df = r.to_pandas()
# df.rename(columns=dict(ra='RA', dec='DEC', source_id='ID'), inplace=True)
# df.dropna(inplace=True)
# df['priority'] = 1.0
# #df[['Mstar', 'z', 'GAL_MU_0_U', 'GAL_MU_R_at_3Re', 'GAL_MU_E_I', 'Re', 'g_m_i', 'GAL_MU_E_Z', 'GAL_MU_E_G', 'GAL_MU_R_at_2Re', 'GAL_MU_E_R', 'SersicIndex_r', 'Dingoflag', 'GAL_MU_0_Z', 'MassHIpred', 'GAL_MAG_I', 'GAL_MU_0_G', 'Ellipticity_r', 'IFU_diam_2Re', 'GAL_MAG_G', 'WALLABYflag', 'GAL_MU_0_I', 'GAL_MU_0_R', 'GAL_MU_E_U']] = 0.0
# # Split the dataframe into three unique parts. One of these will be 'targets', the others will be 'guide' stars, the third will be 'standard' stars
# A, B, C = np.split(df.sample(n=len(df) - len(df) % 3, random_state=42), 3)

# n_samples = len(A)
# if len(A) > 2000:
#     n_samples = 2000

# A.sample(n_samples).to_csv("catalogues/RA_345_stars_A.csv")
# B.to_csv("catalogues/RA_345_stars_B.csv")
# C.to_csv("catalogues/RA_345_stars_C.csv")
