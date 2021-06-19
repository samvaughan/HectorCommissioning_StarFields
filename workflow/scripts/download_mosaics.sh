#!/usr/bin/env bash

set -e

skycell=$1
directory_prefix=$2

for i in $(seq -f "%03g" 0 99)
do
	wget -nc --directory-prefix=$directory_prefix http://ps1images.stsci.edu/rings.v3.skycell/$skycell/$i/rings.v3.skycell.$skycell.$i.stk.r.unconv.fits
done