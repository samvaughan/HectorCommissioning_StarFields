from hop import pipeline
from pathlib import Path
import yaml
import argparse
import shlex

parser = argparse.ArgumentParser()
parser.add_argument("Config_filename")
args = parser.parse_args()

config_filename = Path(args.Config_filename)

with open(config_filename, 'r') as stream:
    try:
        config = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)


# HP = pipeline.HectorPipe(config_filename=config_filename, Profit_files_location=Path("/Volumes/OtherFiles/Science/Panstarrs_Mosaic/SkyMasks/Mosaic/").expanduser())
HP = pipeline.HectorPipe(config_filename=config_filename, Profit_files_location=Path("SkyMasks/Mosaic/").expanduser())

HP.load_input_catalogue()

HP.df_targets['COMPLETED'] = False
HP.df_targets['Tile_number'] = -999
HP.df_guide_stars['type'] = 2
HP.df_standard_stars['type'] = 1

# Add the columns about how many observations we need to complete for each galaxy
# The remaining obsverations column will count down to 0
HP.df_targets['N_observations_to_complete'] = 1
HP.df_targets['remaining_observations'] = HP.df_targets['N_observations_to_complete'].copy()

HP.df_targets["Re"] = 1.0
HP.df_targets['GAL_MU_E_R'] = 19
HP.df_targets['Mstar'] = -99

HP.tile_field(configure_tiles=True, apply_distortion_correction=True, check_sky_fibres=True, date="2021 07 10 14:00") # Time in UTC
HP.allocate_hexabundles_for_single_tile(0) 
