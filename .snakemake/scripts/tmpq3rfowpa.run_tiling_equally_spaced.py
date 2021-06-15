
######## snakemake preamble start (automatically inserted, do not edit) ########
import sys; sys.path.extend(['/opt/anaconda3/envs/snakemake/lib/python3.9/site-packages', '/Users/samvaughan/Science/Hector/Targets/HectorInputCatalogues/Commissioning/StarFields/workflow/scripts']); import pickle; snakemake = pickle.loads(b'\x80\x04\x954\r\x00\x00\x00\x00\x00\x00\x8c\x10snakemake.script\x94\x8c\tSnakemake\x94\x93\x94)\x81\x94}\x94(\x8c\x05input\x94\x8c\x0csnakemake.io\x94\x8c\nInputFiles\x94\x93\x94)\x81\x94(\x8c(config/Field_300_m20_equally_spaced.yaml\x94\x8c1results/StarCatalogues/equally_spaced_300_m22.csv\x94\x8cAresults/StarCatalogues/Panstarrs_match_Gaia_300_m20_standards.csv\x94\x8c>results/StarCatalogues/Panstarrs_match_Gaia_300_m20_guides.csv\x94\x8c5results/SkyMasks/Mosaic/segmap_StarField_300_-22.fits\x94e}\x94(\x8c\x06_names\x94}\x94(\x8c\x0fconfig_filename\x94K\x00N\x86\x94\x8c\x0ctargets_file\x94K\x01N\x86\x94\x8c\x0estandards_file\x94K\x02N\x86\x94\x8c\x0bguides_file\x94K\x03N\x86\x94\x8c\x10profound_skymask\x94K\x04N\x86\x94u\x8c\x12_allowed_overrides\x94]\x94(\x8c\x05index\x94\x8c\x04sort\x94eh\x1e\x8c\tfunctools\x94\x8c\x07partial\x94\x93\x94h\x06\x8c\x19Namedlist._used_attribute\x94\x93\x94\x85\x94R\x94(h$)}\x94\x8c\x05_name\x94h\x1esNt\x94bh\x1fh"h$\x85\x94R\x94(h$)}\x94h(h\x1fsNt\x94bh\x12h\nh\x14h\x0bh\x16h\x0ch\x18h\rh\x1ah\x0eub\x8c\x06output\x94h\x06\x8c\x0bOutputFiles\x94\x93\x94)\x81\x94\x8cBresults/TilingOutputs/300_m20_equally_spaced/completed_tiling.done\x94a}\x94(h\x10}\x94\x8c\x13completed_flag_file\x94K\x00N\x86\x94sh\x1c]\x94(h\x1eh\x1feh\x1eh"h$\x85\x94R\x94(h$)}\x94h(h\x1esNt\x94bh\x1fh"h$\x85\x94R\x94(h$)}\x94h(h\x1fsNt\x94bh5h2ub\x8c\x06params\x94h\x06\x8c\x06Params\x94\x93\x94)\x81\x94}\x94(h\x10}\x94h\x1c]\x94(h\x1eh\x1feh\x1eh"h$\x85\x94R\x94(h$)}\x94h(h\x1esNt\x94bh\x1fh"h$\x85\x94R\x94(h$)}\x94h(h\x1fsNt\x94bub\x8c\twildcards\x94h\x06\x8c\tWildcards\x94\x93\x94)\x81\x94}\x94(h\x10}\x94h\x1c]\x94(h\x1eh\x1feh\x1eh"h$\x85\x94R\x94(h$)}\x94h(h\x1esNt\x94bh\x1fh"h$\x85\x94R\x94(h$)}\x94h(h\x1fsNt\x94bub\x8c\x07threads\x94K\x01\x8c\tresources\x94h\x06\x8c\tResources\x94\x93\x94)\x81\x94(K\x01K\x01e}\x94(h\x10}\x94(\x8c\x06_cores\x94K\x00N\x86\x94\x8c\x06_nodes\x94K\x01N\x86\x94uh\x1c]\x94(h\x1eh\x1feh\x1eh"h$\x85\x94R\x94(h$)}\x94h(h\x1esNt\x94bh\x1fh"h$\x85\x94R\x94(h$)}\x94h(h\x1fsNt\x94bheK\x01hgK\x01ub\x8c\x03log\x94h\x06\x8c\x03Log\x94\x93\x94)\x81\x94}\x94(h\x10}\x94h\x1c]\x94(h\x1eh\x1feh\x1eh"h$\x85\x94R\x94(h$)}\x94h(h\x1esNt\x94bh\x1fh"h$\x85\x94R\x94(h$)}\x94h(h\x1fsNt\x94bub\x8c\x06config\x94}\x94(\x8c\x08filename\x94\x8c(config/Field_300_m20_equally_spaced.yaml\x94\x8c\x05FIELD\x94\x8c\x03300\x94\x8c\tcentre_ra\x94M,\x01\x8c\ncentre_dec\x94J\xea\xff\xff\xff\x8c\x17star_catalogue_width_ra\x94K\x01\x8c\x18star_catalogue_width_dec\x94K\x01\x8c\x1efinal_star_catalogue_name_stem\x94\x8c\x1cPanstarrs_match_Gaia_300_m20\x94\x8c\x15date_for_observations\x94\x8c\x102021-07-10 14:00\x94\x8c\x12target_output_name\x94\x8c?results/StarCatalogues/Panstarrs_match_Gaia_300_m20_targets.csv\x94\x8c\x13guide_star_filename\x94\x8c>results/StarCatalogues/Panstarrs_match_Gaia_300_m20_guides.csv\x94\x8c\x16standard_star_filename\x94\x8cAresults/StarCatalogues/Panstarrs_match_Gaia_300_m20_standards.csv\x94\x8c\x15individual_image_glob\x94\x8c;/Volumes/OtherFiles/Science/Panstarrs_Mosaic/300_m20/*.fits\x94\x8c\x14final_image_filename\x94\x8c-results/Mosaics/Panstarrs_Mosiac_300_m20.fits\x94\x8c\x17image_resolution_arcsec\x94G?\xf8\x00\x00\x00\x00\x00\x00\x8c\x12profound_mask_name\x94\x8c5results/SkyMasks/Mosaic/segmap_StarField_300_-22.fits\x94\x8c\x17profound_catalogue_name\x94\x8c9results/ProfoundCatalogues/profound_catalogue_300_-22.csv\x94\x8c\x11profound_zero_mag\x94G@8\xc2\x8f\\(\xf5\xc3\x8c\nchosen_IDs\x94\x8c?resources/IDs_for_tile_selection/300_m22_equally_spaced_IDs.txt\x94\x8c\x0estarfield_type\x94\x8c\x0eequally_spaced\x94\x8c\tSourceCat\x94\x8c\x14300_m20_comissioning\x94\x8c\x14final_catalogue_name\x94\x8c1results/StarCatalogues/equally_spaced_300_m22.csv\x94\x8c\x14output_filename_stem\x94\x8c\x14300_m20_comissioning\x94\x8c\routput_folder\x94\x8c,results/TilingOutputs/300_m20_equally_spaced\x94\x8c\x0bfresh_start\x94\x88\x8c\tMAX_TRIES\x94M\xe7\x03\x8c\x0fconfigure_field\x94\x88\x8c\x0btiling_type\x94\x8c\x06greedy\x94\x8c\x0fallocation_type\x94\x8c\x06random\x94\x8c\tproximity\x94K\xdc\x8c\x04Nsel\x94K\x13\x8c\x1aN_targets_per_Hector_field\x94K\x13\x8c\x0bNsel_guides\x94Kd\x8c\x0eNsel_standards\x94Kd\x8c\x10TwoDF_FOV_radius\x94G?\xf0\x00\x00\x00\x00\x00\x00\x8c\x0eexclusion_zone\x94K\xdc\x8c\x17Hector_FOV_outer_radius\x94G?\xed\xd1N;\xcd5\xa8\x8c\x17Hector_FOV_inner_radius\x94G\x00\x00\x00\x00\x00\x00\x00\x00\x8c\x15fraction_to_reobserve\x94G\x00\x00\x00\x00\x00\x00\x00\x00\x8c\x11Rescale_proximity\x94\x88\x8c\x18fill_spares_with_repeats\x94\x88\x8c\x1ecolumns_for_target_tile_saving\x94]\x94(\x8c\x02ID\x94\x8c\x02RA\x94\x8c\x03DEC\x94\x8c\x05g_mag\x94\x8c\x05r_mag\x94\x8c\x05i_mag\x94\x8c\x05z_mag\x94\x8c\x05y_mag\x94\x8c\nGAIA_g_mag\x94\x8c\x0bGAIA_bp_mag\x94\x8c\x0bGAIA_rp_mag\x94\x8c\x04pmRA\x94\x8c\x05pmDEC\x94\x8c\x02Re\x94\x8c\nGAL_MU_E_R\x94\x8c\x05Mstar\x94e\x8c\x1dcolumns_for_guide_tile_saving\x94]\x94(\x8c\x02ID\x94\x8c\x02RA\x94\x8c\x03DEC\x94\x8c\x05r_mag\x94\x8c\x04type\x94\x8c\x04pmRA\x94\x8c\x05pmDEC\x94eu\x8c\x04rule\x94\x8c\x19run_tiling_equally_spaced\x94\x8c\x0fbench_iteration\x94N\x8c\tscriptdir\x94\x8ch/Users/samvaughan/Science/Hector/Targets/HectorInputCatalogues/Commissioning/StarFields/workflow/scripts\x94ub.'); from snakemake.logging import logger; logger.printshellcmds = False; __real_file__ = __file__; __file__ = '/Users/samvaughan/Science/Hector/Targets/HectorInputCatalogues/Commissioning/StarFields/workflow/scripts/run_tiling_equally_spaced.py';
######## snakemake preamble end #########
from hop import pipeline
from pathlib import Path
import yaml
import argparse
import shlex

# parser = argparse.ArgumentParser()
# parser.add_argument("Config_filename")
# args = parser.parse_args()

# config_filename = Path(args.Config_filename)

# with open(config_filename, 'r') as stream:
#     try:
#         config = yaml.safe_load(stream)
#     except yaml.YAMLError as exc:
#         print(exc)
config_filename = snakemake.input.config_filename

# HP = pipeline.HectorPipe(config_filename=config_filename, Profit_files_location=Path("/Volumes/OtherFiles/Science/Panstarrs_Mosaic/SkyMasks/Mosaic/").expanduser())
skyfiles_location = Path(snakemake.input.profound_skymask).parent
HP = pipeline.HectorPipe(config_filename=config_filename, Profit_files_location=skyfiles_location)

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

import ipdb; ipdb.set_trace()

HP.tile_field(configure_tiles=True, apply_distortion_correction=True, check_sky_fibres=True, date="2021 07 10 14:00") # Time in UTC
HP.allocate_hexabundles_for_single_tile(0) 

Path(snakemake.output.completed_flag_file).touch()