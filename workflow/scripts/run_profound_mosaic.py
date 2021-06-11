import numpy as np
import matplotlib.pyplot as plt
import subprocess
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


input_image = config['final_image_filename']
output_mask = config['profound_mask_name']
output_catalogue = config['profound_catalogue_name']
zeromag = config['profound_zero_mag']

bash_command = f'./profound_r_script.R {input_image} {output_mask} {output_catalogue} --mag_zero {zeromag}'

print("Running Profound...")
subprocess.check_call(shlex.split(bash_command))
print("Done!")