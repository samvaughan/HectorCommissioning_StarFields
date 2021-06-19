from hop import pipeline
from pathlib import Path
import yaml
import argparse
import shlex

parser = argparse.ArgumentParser()
parser.add_argument("Config_filename")
parser.add_argument("tile_number")
args = parser.parse_args()

config_filename = Path(args.Config_filename)
tile_number = int(args.tile_number)

with open(config_filename, 'r') as stream:
    try:
        config = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)


print(f"Applying the configuration code to Tile {tile_number}")
# HP = pipeline.HectorPipe(config_filename=config_filename, Profit_files_location=Path("/Volumes/OtherFiles/Science/Panstarrs_Mosaic/SkyMasks/Mosaic/").expanduser())
HP = pipeline.HectorPipe(config_filename=config_filename, Profit_files_location=Path("results/SkyMasks/Mosaic/").expanduser())

HP.load_input_catalogue()

# HP.tile_field(configure_tiles=False, apply_distortion_correction=True, check_sky_fibres=True, date=config['date_for_observations']) # Time in UTC
HP.apply_config_code_to_tile(tile_number)
HP.allocate_hexabundles_for_single_tile(tile_number)
