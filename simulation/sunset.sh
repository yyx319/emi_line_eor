#!/usr/bin/env bash
module load intel-fortran/16.0.3
module load intel/2019u5
module load openmpi/intel/11.1/1.4.2

sunset_dir=/home/yy503/Desktop/sunset
sim_dir=/data/ERCblackholes2/smartin/Dwarf1 
pp_dir=/data/ERCblackholes3/yuxuan/dwarf_data/post_processing 
# $pp_dir/RTCRiMHD+SfFb/output_00031/

$sunset_dir/sunset.out -inp $sim_dir/RTCRiMHD+SfFb/output_00050 -out ./sunset_op_u -xmi 0.51556216 -xma 0.51800819 -ymi 0.47619844 -yma 0.47864447 -zmi 0.48841226 -zma 0.49085829 -amx 0.0 -amy 0.0 -amz 1.0 -nvx 0.0 -nvy 1.0 -nvz 0.0 -bnd u_prime -dst $sunset_dir/DataFiles/kext_albedo_WD_MW_3.1_60.txt -nzf 0 -met 1 -m2d 0 









