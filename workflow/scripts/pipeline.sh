set -e 

config_file=Config_Files/Field_300_m20.yaml

python3 make_mosaic.py $config_file
python3 select_star_fields.py $config_file
python3 run_profound.py $config_file