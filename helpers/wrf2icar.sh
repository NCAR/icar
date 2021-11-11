#!/bin/sh
# Script to remove bottom level in selected variables in met_em files
# also rotates wind field to be earth relative instead of grid relative as in met files
#
# Original from Marie Pontoppidan 29/11/2017
# Converted to shell script Ethan Gutmann 26 Jan. 2018

if [ "$#" -lt 2 ]; then
  echo "">&2
  echo "Usage: $0 input_directory output_directory [prefix] [output_file] [temp_prefix]">&2
  echo ""                                                                               >&2
  echo "        If output_file is specified, prefix must be too."                       >&2
  echo "        If temp_prefix is specified, output_file must be too."                  >&2
  echo "        Prefix defaults to 'wrfout'"                                            >&2
  echo "        temp_prefix defaults to 'temp_'"                                        >&2
  echo "        If output_file is specified, all outputs are ncrcated into it."         >&2
  exit 1
fi

if ! [ -d "$1" ]; then
  echo "$1 is not a directory" >&2
  exit 1
fi

input_dir=$1
output_dir=$2

if [ "$#" -ge 3 ]; then
    prefix=$3
else
    prefix='wrfout'
fi

if [ "$#" -eq 5 ]; then
    TEMP=$5
else
    TEMP='temp_'
fi

if [ ${prefix::6} == 'met_em' ]; then
    varlist_3d="PRES,TT,RH,GHT"
    varlist_2d="ALBEDO12M,XLAT_M,XLONG_M,GREENFRAC,HGT_M,LANDMASK,ST,SM,TAVGSFC,SKINTEMP,SOILTEMP"
    uvar="UU"
    vvar="VV"
    general_3d="TT"
    vertical_dim="num_metgrid_levels"
elif [ ${prefix::6} == 'wrfout' ]; then
    varlist_3d="P,PB,T,QVAPOR,PH,PHB"
    varlist_2d="XLAT,XLONG,HGT,LANDMASK,TSLB,SMOIS,TSK,SWDOWN,GLW,PREC_ACC_C"
    uvar="U"
    vvar="V"
    general_3d="T"
    vertical_dim="bottom_top"
fi

sinalpha=SINALPHA
cosalpha=COSALPHA

echo $varlist_3d
echo $varlist_2d
echo ${general_3d},${uvar},${vvar},${cosalpha},${sinalpha}

mkdir -p $output_dir

for f in ${input_dir}/${prefix}*; do
    filename=`basename $f`
    output_file=${output_dir}/${filename}
    echo $filename $output_file

    # One of the primary reasons for this script is to rotate the wind field into earth-relative coordinates
    # instead of having them as model grid relative.
    # First make a file that only has the necessary variables in it
    ncks -v ${general_3d},${uvar},${vvar},${cosalpha},${sinalpha} $f ${TEMP}general.nc

    # Create a U wind variable on the mass grid
    ncap2 -A -s "UU_M=${general_3d}" ${TEMP}general.nc
    # Compute U on the mass grid destaggered from the wind grid
    ncap2 -A -s "UU_M(:,:,:,:)=(${uvar}(:,:,:,0:-2)+${uvar}(:,:,:,1:))/2" ${TEMP}general.nc 2>/dev/null

    # Now repeat for V
    # Create a V wind variable on the mass grid
    ncap2 -A -s "VV_M=${general_3d}" ${TEMP}general.nc
    # Compute V on the mass grid destaggered from the wind grid
    ncap2 -A -s "VV_M(:,:,:,:)=(${vvar}(:,:,0:-2,:)+${vvar}(:,:,1:,:))/2" ${TEMP}general.nc 2>/dev/null

    # Now perform the rotation, recompute
    ncap2 -s "UU_M=UU_M*${cosalpha}-VV_M*${sinalpha}" ${TEMP}general.nc ${TEMP}uvfinal.nc
    ncap2 -A -s "VV_M=VV_M*${cosalpha}+UU_M*${sinalpha}" ${TEMP}general.nc ${TEMP}uvfinal.nc

    if [ ${prefix::6} == 'met_em' ]; then
        ncks -d $vertical_dim,1,-1 -v UU_M,VV_M ${TEMP}uvfinal.nc ${TEMP}new.nc
        ncks -d $vertical_dim,1,-1 -v $varlist_3d $f ${TEMP}lev.nc
    else
        ncks -v UU_M,VV_M ${TEMP}uvfinal.nc ${TEMP}new.nc
        ncks -v $varlist_3d $f ${TEMP}lev.nc
    fi
    ncks -v $varlist_2d $f $output_file
    ncks -A ${TEMP}new.nc $output_file
    ncks -A ${TEMP}lev.nc $output_file

    rm ${TEMP}*

    ncatted -a units,UU_M,m,c,"m s-1" $output_file
    ncatted -a units,VV_M,m,c,"m s-1" $output_file
    ncatted -a description,UU_M,m,c,"earth-rotated wind" $output_file
    ncatted -a description,VV_M,m,c,"earth-rotated wind" $output_file

    if [ ${prefix::6} == 'met_em' ]; then
        ncatted -a scale_factor,RH,c,f,0.01 $output_file
    fi

    # For WRF output files we handle a number of issues to make the files smaller and more intuitive
    # Note that ICAR can actually handle all of these internally, so this is not strictly necessary
    if [ ${prefix::6} == 'wrfout' ]; then

        # First, combine the P and PB (perturbation and base pressure) to save space
        ncap2 -A -s "P=P+PB" $output_file
        ncatted -a description,P,m,c,"pressure" $output_file

        # Second, combine the PH and PHB (perturbation and base geopotential) to save space and convert to height
        ncap2 -A -s "PH=(PH+PHB)/9.81" $output_file
        ncatted -a description,PH,m,c,"height" $output_file
        ncatted -a units,PH,m,c,"m" $output_file

        # Third, add 300 to T (because WRF output files have 300 subtracted)
        ncap2 -A -s "T=T+300" $output_file

        # Fourth, remove the now unnecessary PB and PHB variables
        ncks -x -v PB,PHB $output_file ${output_file}.temp
        mv ${output_file}.temp $output_file

        # Finally, we get rid of the time dimension to XLAT and XLONG to save more space
        # Create a new file without lat and long (which otherwise have a time dimension)
        ncks -C -x -v XLAT,XLONG $output_file ${output_file}.temp
        # Put lat and long WITHOUT a time dimension into a temp file
        ncwa -d Time,0 -a Time -v XLAT,XLONG $output_file ${output_file}.temp2
        # then recombine lat/lon back into the main file so there is no time dimension to them
        mv ${output_file}.temp $output_file
        ncks -A ${output_file}.temp2 $output_file
        rm ${output_file}.temp2

    fi
done


if [ "$#" -ge 4 ]; then
    echo "Concatenating all files together."
    ncrcat ${output_dir}/${prefix}* $4
fi
