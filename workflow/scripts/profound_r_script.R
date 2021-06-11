#!/usr/bin/env Rscript


library("Rfits")
library("ProFound")
library("argparse")

parser = ArgumentParser(description='Run Profound on an input image')
parser$add_argument('image_filename', type="character",
                    help='FITS image')
parser$add_argument('output_mask_name', type="character",
                    help='Output name for the sky mask')
parser$add_argument('output_catalogue_name', type="character",
                    help='Output name for the Profound catalogue')
parser$add_argument('--mag_zero', type="double",
                    help='Zero magnitude in the image', default=24.76)
args = parser$parse_args()

image_filename = args$image_filename
output_mask_name = args$output_mask_name
output_catalogue_name = args$output_catalogue_name
mag_zero = args$mag_zero

image = Rfits_read_image(image_filename, header=TRUE, ext=1)

profound = profoundProFound(image, magzero=mag_zero, plot=FALSE, verbose=TRUE)


#writeFITSim(profound$segim, file=output_mask_name, header=image$header)
Rfits_write_image(profound$segim, filename = output_mask_name)
Rfits_write_header(filename=output_mask_name, keyvalues=image$keyvalues, keycomments=image$keycomments)
write.csv(profound$segstats, file=output_catalogue_name)