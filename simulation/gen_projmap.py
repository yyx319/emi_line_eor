# generate projection map for multiple snapshots and LOS

from operator import index
import sys
sys.path.append('/home/yy503/Desktop/rascas/py')
import jphot as jp
import lya_utils as lya
import numpy as np
from matplotlib import pyplot as plt
import os
from matplotlib.colors import LogNorm
from matplotlib import rcParams
import matplotlib.cm as cm
import matplotlib.colors as colors
from scipy import stats
import importlib #importlib.reload(module)

sys.path.append('/home/yy503/Desktop/simulation_code/ramses_tools/HaloMakerRoutines/')
import NewHaloMakerRoutines as hmr
sys.path.append('/home/yy503/Desktop/simulation_code/ramses_tools/MakeMovies')
import movie_utils

import astropy.units as u; import astropy.constants as c
m_p = c.m_p.cgs.value
cm2pc = (u.cm).to(u.pc)
cm2kpc = (u.cm).to(u.kpc)
s2yr = (u.s).to(u.yr)
s2Gyr = (u.s).to(u.Gyr)
g2Msun = (u.g).to(u.M_sun)

sys.path.append('/home/yy503/Desktop/emi_line_eor/simulation')
import plot_util as pu


cube_region='zoom'
size=6
f_H = 0.76
f_He = 0.24

sim_name = sys.argv[1]
op_idx_a = [51, 52, 53]

nvar = pu.variable_info(sim_name)
w = np.linspace(-6,6, 100)
anal_flag=['HI']

for op_idx in op_idx_a:
    emi_dir = '/data/ERCblackholes3/yuxuan/emi_line/%s/output000%s'%(sim_name, op_idx)
    dat_dir = '/data/ERCblackholes3/yuxuan/dwarf_data/post_processing/%s/output_000%s'%(sim_name, op_idx)
    # read cube: axis0 z axis1 y axis0 x
    lfac, dfac, tfac, redshift, redshiftnum, main_halo_pos, main_halo_vel = pu.get_sim_info(sim_name, op_idx )
    for i in [1,nvar-2]:
        filename='/data/ERCblackholes3/yuxuan/dwarf_data/post_processing/%s/\
output_000%s/%s_cube_%d.dat'%(sim_name, op_idx, cube_region, i)
        locals()['cube%d'%i], info = movie_utils.readcube_from_file( filename )
        #locals()['cube%d'%i][ locals()['cube%d'%i]<0 ] = 0
    cube_fHII = eval( 'cube%d'%(nvar-2) )
    nx = info['nx']
    dx = size/nx/cm2kpc # in cm
    cube1 = cube1*dfac
    cube_mass = cube1*dx**3*g2Msun # in M_sun

    if 'RT' in sim_name:
        cube_nHI = cube1/m_p*f_H*(1-cube_fHII)
        cube_HI_mass = cube_nHI*m_p*dx**3*g2Msun


    nx, ny, nz=np.shape(cube_nHI)
    x_3d,y_3d,z_3d = np.meshgrid(np.linspace(-6,6,nx),np.linspace(-6,6,ny),np.linspace(-6,6,nz), indexing='ij')
    
    LOS_a = np.loadtxt('%s/LOS.txt'%emi_dir)
    for i, LOS in enumerate(LOS_a):
        print('generating projection map for LOS %d'%i)
        wx, wy = pu.cal_projection_pos(x_3d, y_3d, z_3d, LOS )
        if 'HI' in anal_flag:
            nHI_p = stats.binned_statistic_2d(wx.flatten(), wy.flatten(), cube_nHI.flatten(), bins=[w,w])
            np.savetxt('%s/%s_HI_dist_LOS_%d.txt'%(dat_dir, op_idx, i), nHI_p.statistic)
        if 'Ndens' in anal_flag:
            nHI_p = stats.binned_statistic_2d(wx.flatten(), wy.flatten(), cube1.flatten(), bins=[w,w])
            np.savetxt('%s/%s_Ndens_dist_LOS_%d.txt'%(dat_dir, op_idx, i), nHI_p.statistic)