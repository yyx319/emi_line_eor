############################
#
#
# python gen_analysis_script.py $sim_name
############################

from logging.handlers import NTEventLogHandler
import os
import sys
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import rcParams
import matplotlib.cm as cm
import matplotlib.colors as colors


import yt
import importlib #importlib.reload(module)

sys.path.append('/data/ERCblackholes3/yuxuan')

from scipy import stats
sys.path.append('/home/yy503/Desktop/simulation_code/ramses_tools/HaloMakerRoutines/')
import NewHaloMakerRoutines as hmr
sys.path.append('/home/yy503/Desktop/simulation_code/ramses_tools/MakeMovies')
import movie_utils
from support_functions import *
#import basic_analysis
import plot_util as pu
import astropy.units as u; import astropy.constants as c
m_p = c.m_p.cgs.value
cm2pc = (u.cm).to(u.pc)
cm2kpc = (u.cm).to(u.kpc)
s2yr = (u.s).to(u.yr)
s2Gyr = (u.s).to(u.Gyr)
g2Msun = (u.g).to(u.M_sun)


mode=['cube'] #cube map sfr
cube_region='zoom'

sim_name=sys.argv[1]


snapnum_a, mainhalo_pos_a, mainhalo_vel_a, gashalo_pos_a, gashalo_vel_a = pu.get_halo_info(sim_name)

nvar = pu.variable_info(sim_name)

typ_array = [1, nvar-2] #np.arange(1, nvar+1)


sim_dir = '/data/ERCblackholes2/smartin/Dwarf1'
pp_dir = '/data/ERCblackholes3/yuxuan/dwarf_data/post_processing'
ramses_dir = '/home/yy503/Desktop/simulation_code/ramses_rtcrmhd'

op_a2 = np.arange(51, 61, 1) # snapshots for which we perform analysis 

# write ramses_first_look.sh
with open('ramses_first_look.sh', 'w') as f:
    f.write('sim_dir=%s \n'%sim_dir )
    f.write('ramses_dir=%s \n'%ramses_dir )
    f.write('pp_dir=%s \n'%pp_dir)
    for snapnum in op_a2:
        print(mode)
        f.write('# %d \n'%snapnum )
        snapnum = int(snapnum)
        dat_dir = '%s/output_000%d'%(sim_name, snapnum )
        lfac, dfac, tfac, redshift, redshiftnum, main_halo_pos, main_halo_vel, track_pos_halo, gas_halo_pos, gas_halo_vel = pu.get_sim_info(sim_name, snapnum )
        print(snapnum)
        if cube_region == 'zoom':
            size = 6/lfac/cm2kpc # 6 kpc
            lma = 18
        elif cube_region == 'full':
            size = 100/lfac/cm2kpc # 100 kpc box
            lma = 14

        xmin, ymin, zmin = track_pos_halo - size/2
        xmax, ymax, zmax = track_pos_halo + size/2
        if os.path.isdir('%s/%s'%(pp_dir, sim_name) )==0:
            f.write('mkdir %s/%s \n'%(pp_dir, sim_name) )
        if os.path.isdir('%s/%s/output_000%d'%(pp_dir, sim_name, snapnum) )==0:
            f.write('mkdir %s/%s/output_000%d \n'%(pp_dir, sim_name, snapnum) )
        
        # for each snapshot generate map or cube for 
        for typ in typ_array:
            typ = int(typ)        
            if 'cube' in mode:
                f.write('$ramses_dir/utils/f90/amr2cube -inp $sim_dir/%s -out $pp_dir/%s/output_000%d/%s_cube_%d.dat -xmi %.8f -xma %.8f \
-ymi %.8f -yma %.8f -zmi %.8f -zma %.8f -lma %d -typ %d \n'
                        %(dat_dir, sim_name, snapnum, cube_region, typ, xmin, xmax, ymin, ymax, zmin, zmax, lma, typ) )
            if 'map' in mode:
                for incl_ang in ['z']:
                    f.write('$ramses_dir/utils/f90/amr2map -inp $sim_dir/%s -out $pp_dir/%s/output_000%d/%s_map_%d_%s.dat -dir %s -xmi %.8f -xma %.8f \
-ymi %.8f -yma %.8f -zmi %.8f -zma %.8f -lma %d -typ %d \n'
                        %(dat_dir, sim_name, snapnum, cube_region, typ, incl_ang, incl_ang, xmin, xmax, ymin, ymax, zmin, zmax, lma, typ) )
        if 'sfr' in mode:
            f.write('$ramses_dir/utils/f90/part2sfr -inp $sim_dir/%s -out $pp_dir/%s/output_000%d/%s_sfr.dat -xmi %.8f -xma %.8f -ymi %.8f -yma %.8f -zmi %.8f -zma %.8f \n'
                %(dat_dir, sim_name, snapnum, cube_region, xmin, xmax, ymin, ymax, zmin, zmax) )
            
    f.write('ls $pp_dir/%s/output_000%d/'%(sim_name, snapnum) )
    f.close()