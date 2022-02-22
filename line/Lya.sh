#!/usr/bin/env bash
module load intel-fortran/16.0.3
module load intel/2019u5
module load openmpi/intel/11.1/1.4.2

#1 sim_name: RTCRiMHD+SfFb
#2 machine name: calx089, 151, 156

run=parallel
for snapnum in {41..46..1}
  do
    python gen_conf_file.py $1 $snapnum

    echo $1
    echo $snapnum
    echo $2

    if [[ $2 == 'calx089' ]]; then
        rascas_dir=/data/ERCblackholes3/yuxuan/rascas-sdpo_calx089
        nmpi=30
    elif [[ $2 == 'calx151' ]]; then
        rascas_dir=/data/ERCblackholes3/yuxuan/rascas-sdpo_calx151
        nmpi=20
    elif [[ $2 == 'calx156' ]]; then
        rascas_dir=/home/yy503/Desktop/rascas-sdpo
        nmpi=70
    fi

    $rascas_dir/f90/CreateDomDump Lya/CDD.conf

    $rascas_dir/f90/LyaPhotonsFromGas Lya/ppic_col.conf
    $rascas_dir/f90/LyaPhotonsFromGas Lya/ppic_rec.conf
    echo finish ppic

    if [[ $run == 'parallel' ]]; then
        mpiexec -n $nmpi $rascas_dir/f90/rascas Lya/RaS_col.conf
        mpiexec -n $nmpi $rascas_dir/f90/rascas Lya/RaS_rec.conf
    else
        $rascas_dir/f90/rascas-serial Lya/RaS_col.conf
        $rascas_dir/f90/rascas-serial Lya/RaS_rec.conf
    fi;
  done
