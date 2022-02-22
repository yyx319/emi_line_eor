# general import
import os
import sys
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import rcParams
import matplotlib.cm as cm
import matplotlib.colors as colors

sys.path.append('/data/ERCblackholes3/yuxuan')
from scipy import stats

sys.path.append('/home/yy503/Desktop/simulation_code/ramses_tools/HaloMakerRoutines/')
import NewHaloMakerRoutines as hmr
sys.path.append('/home/yy503/Desktop/simulation_code/ramses_tools/MakeMovies')
import movie_utils
from support_functions import *
#import basic_analysis

import astropy.units as u; import astropy.constants as c
m_p = c.m_p.cgs.value
cm2pc = (u.cm).to(u.pc)
cm2kpc = (u.cm).to(u.kpc)
s2yr = (u.s).to(u.yr)
s2Gyr = (u.s).to(u.Gyr)
g2Msun = (u.g).to(u.M_sun)


dat_dir = '/data/ERCblackholes3/yuxuan/dwarf_data'
#dat_dir = '/data/ERCblackholes2/smartin/Dwarf1'
prj_dir = '/home/yy503/Desktop/emi_line_eor'
from matplotlib.colors import LogNorm
import pynbody
import osyris

import healpy
import plot_util as pu
import track_halo 
import imp



sim_name='RT+SfFb'
op_idx = 42
op_idx_a, main_halo_pos_a, main_halo_vel_a = pu.get_halo_info(sim_name)
lfac, dfac, tfac, redshift, redshiftnum, main_halo_pos, main_halo_vel= pu.get_sim_info(sim_name, op_idx)

plt.subplot(1,2,1)
plt.plot(op_idx_a, main_halo_pos_a)
plt.subplot(1,2,2)
plt.plot(op_idx_a, main_halo_vel_a*lfac/tfac/1e5) #km/s






# read cube: axis0 z axis1 y axis0 x
cube_region='zoom'
for i in [1,8]:
    filename='/data/ERCblackholes3/yuxuan/dwarf_data/post_processing/%s/\
output_000%s/%s_cube_%d.dat'%(sim_name, op_idx, cube_region, i)
    locals()['cube%d'%i], info = movie_utils.readcube_from_file( filename )
    #locals()['cube%d'%i][ locals()['cube%d'%i]<0 ] = 0

size=6
nx = info['nx']
dx = size/nx/cm2kpc # in cm

f_H = 0.76
f_He = 0.24
cube1 = cube1*dfac
cube_mass = cube1*dx**3*g2Msun # in M_sun
cube2 = cube2*lfac/tfac/1e5 # in km/s
cube3 = cube3*lfac/tfac/1e5
cube4 = cube4*lfac/tfac/1e5


if 'RT' in sim_name:
    #cube12 = cube12*dfac*lfac**2/tfac**2
    cube_nHI = cube1/m_p*f_H*(1-cube8)
    cube_HI_mass = cube_nHI*m_p*dx**3*g2Msun
    cube_ndens = cube1/m_p * ( f_H*(1+cube15) + f_He/4*(1+cube16+2*cube17) )
    cube_T = cube12 / cube_ndens /c.k_B.cgs.value


############
# map of column density , temperature, 
############
plt.figure(figsize=(22,12))
plt.subplot(2,3,1)
NHI = np.sum(cube_nHI*dx, axis=0 )
plt.imshow( NHI, norm=LogNorm(), origin='lower',extent=[-3,3, -3,3])
cbar = plt.colorbar()
cbar.set_label(r'$N_{\rm HI}$ [cm$^{-2}$]')
plt.contour(NHI, [1e17, 2e20], colors =['k','red'], origin='lower',extent=[-3,3, -3,3] )


plt.subplot(2,3,2)
xHI_m = np.sum( (1-cube8)*cube_nHI, axis=0 )/np.sum( cube_nHI, axis=0 )
plt.imshow( xHI_m, norm=LogNorm(),origin='lower', extent=[-3,3, -3,3])
cbar = plt.colorbar()
cbar.set_label(r'$x_{\rm HI, m}$')
plt.xlabel(r'z [kpc]')
plt.ylabel(r'y [kpc]')
plt.contour(NHI, [1e17, 2e20], colors =['k','red'], origin='lower',extent=[-3,3, -3,3] )

plt.subplot(2,3,3)
xHI_v = np.mean( 1-cube8, axis=0 )
plt.imshow( xHI_v, norm=LogNorm(), origin='lower',extent=[-3,3, -3,3])
cbar = plt.colorbar()
cbar.set_label(r'$x_{\rm HI, V}$')
plt.contour(NHI, [1e17, 2e20], colors =['k','red'], origin='lower',extent=[-3,3, -3,3] )

plt.subplot(2,3,4)
plt.imshow( xHI_m/xHI_v, norm=LogNorm(vmin=1, vmax=1e2), origin='lower',extent=[-3,3, -3,3])
cbar = plt.colorbar()
cbar.set_label(r'$x_{\rm HI, m}/x_{\rm HI, V}$')

plt.contour(NHI, [1e17, 2e20], colors =['k','red'], origin='lower',extent=[-3,3, -3,3] )

plt.subplot(2,3,5)
#plt.imshow( np.sum(cube_T*cube_ndens, axis=0 )/np.sum(cube_ndens, axis=0) , norm=LogNorm(), 
#           extent=[-size/2, size/2, -size/2, size/2], origin='lower')
#cbar = plt.colorbar()
#cbar.set_label(r'$T$ [K]')
#plt.contour(NHI, [1e17], colors =['k'], extent=[-size/2, size/2, -size/2, size/2] )

plt.subplot(2,3,6)
#N = np.sum(cube_ndens*dx, axis=0 )
#plt.imshow( N, norm=LogNorm(), extent=[-size/2, size/2, -size/2, size/2], origin='lower')
#cbar = plt.colorbar()
#cbar.set_label(r'$N$ [cm$^{-2}$]')

plt.savefig('../figure/%s/zoom_HI_dist_%s.pdf'%(sim_name, op_idx) )








# velocity map and PV diagram
## slice
inter=30
x = np.linspace(-3,3, nx)
cube_x, cube_y, cube_z = np.meshgrid( x,x,x, indexing='ij')
plt.rcParams['font.size']='20'

y_bin = np.linspace(-3,3,60)
v_bin = np.linspace(-150, 150, 100)

plt.figure(figsize=(28,18))

axis=0
width = 1.0
loc_a = [-1.5, 0, 1.5]

plt.subplot(2,3,1)
idx = np.where( (x>loc_a[0]-1/2*width ) & (x<loc_a[0]+1/2*width ) )[0]
plt.imshow( np.mean(cube_ndens[idx,:,:], axis=axis) , extent=[-3,3,-3,3], norm=LogNorm(), origin='lower' )
cbar = plt.colorbar()
cbar.set_label(r'n [cm$^{-3}$]')
plt.quiver( x[::inter],x[::inter], np.mean(cube4[idx, ::inter, ::inter], axis=axis), \
           np.mean( cube3[idx, ::inter, ::inter], axis=axis) )
plt.xlabel(r'z [kpc]')
plt.ylabel(r'y [kpc]')
plt.title(r'%.1f kpc'%loc_a[0] )


plt.subplot(2,3,4)
m_pv = stats.binned_statistic_2d(cube_y[idx,:,:].flatten(), cube4[idx,:,:].flatten(),\
                              cube_mass[idx,:,:].flatten(), 'sum', [y_bin, v_bin])
                                                                       
plt.imshow(m_pv.statistic, norm=LogNorm(), origin='lower', extent=[-150, 150, -3, 3], aspect='auto')
cbar = plt.colorbar()
cbar.set_label(r'mass [g]')
plt.xlabel(r'$v_z$ [km/s]')
plt.ylabel(r'y [kpc]')

plt.subplot(2,3,2)
idx = np.where( (x>loc_a[1]-1/2*width ) & (x<loc_a[1]+1/2*width ) )[0]
plt.imshow( np.mean( cube_ndens[idx,:,:], axis=axis), extent=[-3,3,-3,3], norm=LogNorm(), origin='lower' )
plt.colorbar()
plt.quiver( x[::inter],x[::inter], np.mean(cube4[idx, ::inter, ::inter], axis=axis ), \
           np.mean(cube3[idx, ::inter, ::inter], axis=axis ) )
plt.title(r'%.1f kpc'%loc_a[1] )

plt.subplot(2,3,5)
# slice with fixed x
m_pv = stats.binned_statistic_2d(cube_y[idx,:,:].flatten(), cube4[idx,:,:].flatten(),\
                              cube_mass[idx,:,:].flatten(), 'sum', [y_bin, v_bin])
                                                                       
plt.imshow(m_pv.statistic, norm=LogNorm(), origin='lower', extent=[-150, 150, -3, 3], aspect='auto')
plt.colorbar()

plt.subplot(2,3,3)
idx = np.where( (x>loc_a[2]-1/2*width ) & (x<loc_a[2]+1/2*width ) )[0]
plt.imshow( np.mean(cube_ndens[idx,:,:], axis=axis), extent=[-3,3,-3,3], norm=LogNorm(), origin='lower' )
plt.colorbar()
plt.quiver( x[::inter],x[::inter], np.mean(cube4[idx, ::inter, ::inter], axis=axis),\
            np.mean( cube3[idx, ::inter, ::inter], axis=axis ) )
plt.title(r'%.1f kpc'%loc_a[2] )


plt.subplot(2,3,6)
m_pv = stats.binned_statistic_2d(cube_y[idx,:,:].flatten(), cube4[idx,:,:].flatten(),\
                              cube_mass[idx,:,:].flatten(), 'sum', [y_bin, v_bin])
                                                                       
plt.imshow(m_pv.statistic, norm=LogNorm(), origin='lower', extent=[-150, 150, -3, 3], aspect='auto')
plt.colorbar()
# PV plot
plt.savefig('figure/%s/pp_pv_plot_%s.pdf'%(sim_name, op_idx) )
 


  
#
# radial profile 
#
vxm, _, _ = stats.binned_statistic(r, vx*mass, 'sum', r_bin)
vxsqm, _, _ = stats.binned_statistic(r, vx**2*mass, 'sum', r_bin)
E_vx = vxm/m
E_vx2 = vxsqm/m
sigma_vx = np.sqrt(E_vx2-E_vx**2)

vym, _, _ = stats.binned_statistic(r, vy*mass, 'sum', r_bin)
vysqm, _, _ = stats.binned_statistic(r, vy**2*mass, 'sum', r_bin)
E_vy = vym/m
E_vy2 = vysqm/m
sigma_vy = np.sqrt(E_vy2-E_vy**2)

vzm, _, _ = stats.binned_statistic(r, vz*mass, 'sum', r_bin)
vzsqm, _, _ = stats.binned_statistic(r, vz**2*mass, 'sum', r_bin)
E_vz = vzm/m
E_vz2 = vzsqm/m
sigma_vz = np.sqrt(E_vz2-E_vz**2)

plt.plot(r_bin[1:], sigma_vx/1e5, label = 'x')
plt.plot(r_bin[1:], sigma_vy/1e5, label = 'y')
plt.plot(r_bin[1:], sigma_vz/1e5, label = 'z')
plt.plot(r_bin[1:], np.sqrt( sigma_vx**2 + sigma_vy**2 + sigma_vz**2 )/1e5, label = sim_name)
plt.xlabel(r'r [kpc]')
plt.ylabel(r'$\sigma_{\rm g}$ [km/s]')
plt.legend()
plt.savefig('%s/figure/%s/%s_radial_sigma_gas.png'%(prj_dir, sim_name, sim_name)  )
plt.clf()
