from lib2to3.pygram import python_grammar_no_print_statement
import os
import sys
import numpy as np 
import matplotlib.pyplot as plt 
import yt
sys.path.append('/data/ERCblackholes3/yuxuan')
import astropy.units as u; import astropy.constants as c
from scipy import stats
sys.path.append('/home/yy503/Desktop/simulation_code/ramses_tools/HaloMakerRoutines/')
import NewHaloMakerRoutines as hmr
sys.path.append('/home/yy503/Desktop/simulation_code/ramses_tools/MakeMovies')
import movie_utils
from support_functions import *

cm2pc = (u.cm).to(u.pc)
cm2kpc = (u.cm).to(u.kpc)
s2yr = (u.s).to(u.yr)
s2Gyr = (u.s).to(u.Gyr)
g2Msun = (u.g).to(u.M_sun)



dat_dir='/data/ERCblackholes3/yuxuan/dwarf_data'
sim_dir = '/data/ERCblackholes2/smartin/Dwarf1'

def variable_info(sim_name):
    nvar = 7 # number of basic variable: density, vx, vy, vz, P_th, pass_sca1, pass_sca2
    if 'MHD' in sim_name:
        nvar += 6 # B_left (xyz), B_right (xyz)
    if 'RT' in sim_name:
        nvar += 3 # pass_sca3,4,5
    if 'CR' in sim_name:
        nvar += 1 # P_CR
    return nvar

def LOS2angle(LOS_a):
    # test: [1,0,0] -> [pi/2, 0], [0,1,0] -> [pi/2, pi/2]
    theta_a = []
    phi_a = []
    for LOS in LOS_a:
        sintheta = np.sqrt(LOS[0]**2 + LOS[1]**2 )
        theta = np.arctan2( sintheta, LOS[2] )
        phi = np.arctan2( LOS[1], LOS[0] )
        theta_a.append(theta)
        phi_a.append(phi)
    return theta_a, phi_a

def cal_projection_pos(x_mesh, y_mesh, z_mesh, LOS):
    # projection map along arbitrary LOS

    #unit vector for the axis, defined in the same way as kobs_perp_1,2 in module_mock.f90, 
    # k_perp1 = k_obs cross xhat then normalize to unit; k_perp2 = kobs cross k_perp1
    if np.array_equiv(LOS, [1,0,0]): 
        k_perp1=[1,0,0]
        k_perp2=[0,1,0]
    else:
        k_perp1 = np.cross(LOS, [1,0,0] )
        k_perp1 = k_perp1/np.sqrt( k_perp1.dot(k_perp1) )
        k_perp2 = np.cross(LOS, k_perp1)
        # assign projected coordinate sx sy for each cell
    print(k_perp1, k_perp2)
    # projected coordinate
    wx = x_mesh*k_perp1[0]+ y_mesh*k_perp1[1]+ z_mesh*k_perp1[2]
    wy = x_mesh*k_perp2[0]+ y_mesh*k_perp2[1]+ z_mesh*k_perp2[2]
    return wx, wy


def cal_sph_angles(x_mesh, y_mesh, z_mesh, center):
    # given the center, calculate two angles in spherical coordinate of cells. 
    x_mesh -= center[0]
    y_mesh -= center[1]
    z_mesh -= center[2]
    r_mesh = np.sqrt( x_mesh**2 + y_mesh**2 + z_mesh**2  )
    theta_mesh = np.arctan( np.sqrt( x_mesh**2+y_mesh**2) / z_mesh )
    phi_mesh = np.arctan2( y_mesh, x_mesh )
    return r_mesh, theta_mesh, phi_mesh

##################################
# get basic info from simulation #
##################################

# time evolution of snapshots in a simulation 
def get_halo_info(sim_name):
    # read 
    op_idx_a = np.loadtxt('%s/post_processing/%s/op_index.txt'%(dat_dir, sim_name) )
    main_halo_pos_a = np.loadtxt('%s/post_processing/%s/main_halo_pos.txt'%(dat_dir, sim_name) )
    main_halo_vel_a = np.loadtxt('%s/post_processing/%s/main_halo_vel.txt'%(dat_dir, sim_name) )
    gas_halo_pos_a = np.loadtxt('%s/post_processing/%s/gas_halo_pos.txt'%(dat_dir, sim_name) )
    gas_halo_vel_a = np.loadtxt('%s/post_processing/%s/gas_halo_vel.txt'%(dat_dir, sim_name) )
    return op_idx_a, main_halo_pos_a, main_halo_vel_a, gas_halo_pos_a, gas_halo_vel_a

def get_sim_info(sim_name, op_idx):
    op_idx_a, main_halo_pos_a, main_halo_vel_a, gas_halo_pos_a, gas_halo_vel_a = get_halo_info(sim_name)
    idx = np.where(op_idx_a==int(op_idx) )[0][0]

    main_halo_pos = main_halo_pos_a[idx]
    main_halo_vel = main_halo_vel_a[idx]
    gas_halo_pos = gas_halo_pos_a[idx]
    gas_halo_vel = gas_halo_vel_a[idx]

    hydro_file_descriptor = '%s/%s/output_000%s/hydro_file_descriptor.txt'%(sim_dir, sim_name, op_idx)
    with open(hydro_file_descriptor, 'r') as f:
        print(f.read())
        
    info_file = '%s/%s/output_000%s/info_000%s.txt'%(sim_dir, sim_name, op_idx, op_idx)
    lfac, dfac, tfac, redshift, redshiftnum = movie_utils.read_infofile(info_file)

    return lfac, dfac, tfac, redshift, redshiftnum, main_halo_pos, main_halo_vel, gas_halo_pos, gas_halo_vel



def ID_dmhalo(part_ID, part_Z):
    dm_ID = part_ID[part_Z==0]
    n, bins, patches = plt.hist(dm_ID, bins=100)
    
    return n, bins, patches









prj_dir = '/home/yy503/Desktop/emi_line_eor'



def read_data_yt(sim_name, op_idx, center, radius, mode='complex'):
    load_fields = dlff.define_loading_fields_funct(sim_name)
    print(load_fields)
    filename = '%s/%s/output_%s/info_%s.txt'%(dat_dir, sim_name, op_idx, op_idx)
    extra_fields = [('tform', 'd'), ('metal', 'd'), ('imass', 'd') ]
    ds = yt.load(filename, fields = load_fields, extra_particle_fields = extra_fields)
    redshift = ds.current_redshift
    time = (ds.current_time*ds.time_unit*s2Gyr).value

    #tform = vir_sph[('io', 'tform')]
    part_x = vir_sph[('all', 'particle_position_x')].value * (13.43e3/(1+ds.current_redshift) ) - center[0]
    part_y = vir_sph[('all', 'particle_position_y')].value * (13.43e3/(1+ds.current_redshift) ) - center[1]
    part_z = vir_sph[('all', 'particle_position_z')].value * (13.43e3/(1+ds.current_redshift) ) - center[2]
    part_r = np.sqrt( part_x**2 + part_y**2 + part_z**2 ) # in kpc
    part_mass = vir_sph[('all', 'particle_mass')]
    part_age = vir_sph[('all', 'particle_age')]
    #part_r = vir_sph[('all', 'particle_radius')]
    part_cpi = vir_sph[('all', 'particle_position_cylindrical_radius')]
    part_Z = vir_sph[('all', 'particle_metallicity')]
    part_ID = vir_sph[('all', 'particle_index') ]

    # star particle 
    idx_star = np.where(part_Z!=0)[0]
    cpi_star = part_cpi[idx_star].value*cm2kpc/(1+redshift) #in kpc
    age_star = part_age[idx_star].value*s2Gyr # in Gyr
    mass_star = part_mass[idx_star].value*g2Msun # in Msun
 

    if 'RT' in sim_name:
        return ds, redshift, time, dens, temp, mass, vol, vx, vy, vz, vr, \
            x, y, z, r, Z, xHII, part_ID, part_x, part_y, part_z, \
            part_r, part_mass, part_age, part_Z, cpi_star, age_star, mass_star
    else:
        return ds, redshift, time, dens, temp, mass, vol, vx, vy, vz, vr, \
            x, y, z, r, Z, part_ID, part_x, part_y, part_z, \
            part_r, part_mass, part_age, part_Z, cpi_star, age_star, mass_star



# radial profile of number of cells
# dens prof (volume weighted)


def mdot(vr, mass, r, r_bin, v_halo, geometry='sphere', in_out='outflow', v_cut=0):
    # mass outflow rate 
    v_cut*=1e5
    vr = vr - v_halo
    dL = r_bin[1:]-r_bin[:-1]
    r_a = 1/2*(r_bin[1:]+r_bin[:-1])
    mdot_a = np.zeros( len(dL) )
    
    print(v_cut, vr, dL, r_a)
    for i in range( len(dL) ):
        if in_out == 'outflow':
            idx = np.where( (r>r_bin[i]) & (r<r_bin[i+1]) & (vr>v_cut) )[0]
        if in_out == 'inflow':
            idx = np.where( (r>r_bin[i]) & (r<r_bin[i+1]) & (vr<v_cut) )[0]
        mdot_a[i] = np.sum( vr[idx]*mass[idx] ) /dL[i] # in gcm/s/kpc
        
    mdot_a = (mdot_a*u.g*u.cm/u.s/u.kpc).to(u.M_sun/u.yr) # in M_sun/yr
    return r_a, mdot_a
  



def KS_relation(cpi_star, age_star, mass_star, t_ave, cpi, mass, r_bin):
    young_idx = np.where( age_star<t_ave/1000 )[0]
    mtot,_,_ = stats.binned_statistic( cpi_star[young_idx], mass_star[young_idx], 'sum', r_bin )
    S_a = np.pi* (r_bin[1:]**2 - r_bin[:-1])
    sfr_2d = mtot/t_ave/1e6/S_a # in M_sun/yr/kpc^2

    mass = mass*g2Msun
        
    mtot_gas,_,_ = stats.binned_statistic( cpi, mass, 'sum', r_bin )
    surf_dens = mtot_gas/S_a/1e6 # M_sun/pc2
    return sfr_2d, surf_dens

def rotation_curve(r, mass, part_r, part_mass, r_bin):
    m, _, _ = stats.binned_statistic(r, mass, 'sum', r_bin)
    m_part, _, _ = stats.binned_statistic(part_r, part_mass, 'sum', r_bin)
    
    M_r = np.zeros( len(r_bin) ) 
    M_r_g = np.zeros( len(r_bin) )
    M_r_p = np.zeros( len(r_bin) )
    for i in range( len(r_bin) ): 
        M_r[i] = np.sum( m[:i] + m_part[:i] )
        M_r_g[i] = np.sum( m[:i] )
        M_r_p[i] = np.sum( m_part[:i] )

    v_c = np.sqrt(c.G*M_r[1:]*u.g/r_bin[1:]/u.kpc).to(u.km/u.s)
    v_c_g = np.sqrt(c.G*M_r_g[1:]*u.g/r_bin[1:]/u.kpc).to(u.km/u.s)
    v_c_p = np.sqrt(c.G*M_r_p[1:]*u.g/r_bin[1:]/u.kpc).to(u.km/u.s)
    
    r_a = 1/2*(r_bin[1:]+r_bin[:-1])
    return r_a, v_c_g, v_c_p, v_c







    