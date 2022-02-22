from operator import ge
import sys
sys.path.append('/home/yy503/Desktop/rascas/py')
import os
import numpy as np
from random import seed
from random import random
import astropy.units as u; import astropy.constants as c
sys.path.append('/home/yy503/Desktop/emi_line_eor/simulation')
import plot_util as pu

m_p = c.m_p.cgs.value
cm2pc = (u.cm).to(u.pc)
cm2kpc = (u.cm).to(u.kpc)
s2yr = (u.s).to(u.yr)
s2Gyr = (u.s).to(u.Gyr)
g2Msun = (u.g).to(u.M_sun)

def write_params_CDD(file, dat_dir, sim_dir, snapnum, main_halo_pos, rsp, overwrite='F'):
    f = open(file,'w')
    f.write('[CreateDomDump] \n')
    f.write('DomDumpDir = %s/ \n'%dat_dir)
    f.write('repository = %s/  \n'%sim_dir)
    f.write('snapnum = %s \n'%snapnum)
    f.write('comput_dom_type = sphere \n')
    f.write('comput_dom_pos = %f %f %f \n'%(main_halo_pos[0], main_halo_pos[1], main_halo_pos[2]) )
    f.write('comput_dom_rsp = %e \n'%rsp )
    f.write('decomp_dom_type = sphere \n')
    f.write('decomp_dom_ndomain = 1 \n')
    f.write('decomp_dom_xc = %f \n'%main_halo_pos[0] ) 
    f.write('decomp_dom_yc = %f \n'%main_halo_pos[1] )
    f.write('decomp_dom_zc = %f \n'%main_halo_pos[2] )
    f.write('decomp_dom_rsp = %e \n'%rsp )
    f.write('verbose = T \n')
    f.write('reading_method = hilbert \n')      
    f.write('[gas_composition] \n') 
    f.write('f_ion = 1.0000000000000000e-02 \n')       
    f.write('Zref = 5.0000000000000001e-03 \n')
    f.write('gas_overwrite = %s \n'%overwrite) 
    f.write('verbose = T \n')

    f.write('[ramses] \n')                               
    f.write('self_shielding = F \n')
    f.write('ramses_rt = T \n') 
    f.write('verbose = T \n') 
    f.write('use_initial_mass = T \n')  
    f.write('cosmo = T \n')
    f.write('use_proper_time = T \n') 
    f.write('read_rt_variables = T \n')
    if 'RTCRiMHD+SfFb' in sim_dir:
        f.write('  itemp  = 12 \n')
        f.write('  imetal = 14 \n')
        f.write('  ihii   = 15 \n')
        f.write('  iheii  = 16 \n')
        f.write('  iheiii = 17 \n')

    elif 'RT+SfFb' in sim_dir:
        f.write('  itemp  =  5\n')
        f.write('  imetal =  7\n')
        f.write('  ihii   =  8\n')
        f.write('  iheii  =  9\n')
        f.write('  iheiii =  10\n')

    elif 'RTiMHD+SfFb' in sim_dir:
        f.write('  itemp  = 11 \n')
        f.write('  imetal = 13 \n')
        f.write('  ihii   = 14 \n')
        f.write('  iheii  = 15 \n')
        f.write('  iheiii = 16 \n')
    f.close()


## LyaPhotonFromGas
def write_params_LyaPhotonFromGas(file, dat_dir, sim_dir, snapnum, main_halo_pos, rsp, nphotons, rec='False', col='True'):
    f = open(file,'w')
    f.write('[LyaPhotonsFromGas] \n')
    if rec=='True':
        f.write('  outputfileRec = %s/rec.IC \n'%dat_dir ) 
    elif col=='True':
        f.write('  outputfileCol = %s/col.IC \n'%dat_dir ) 
    f.write('  repository = %s/ \n'%sim_dir )
    f.write('  snapnum = %s \n'%snapnum )
    f.write('  emission_dom_type = sphere \n')
    f.write('  emission_dom_pos = %f %f %f \n'%(main_halo_pos[0], main_halo_pos[1], main_halo_pos[2])  )
    f.write('  emission_dom_rsp = %e \n'%rsp )
    if rec=='True':
        f.write('  doRecombs = True \n')
        f.write('  doColls = False \n')
    elif col=='True':
        f.write('  doRecombs = False \n')
        f.write('  doColls = True \n')
    f.write('  tcool_resolution = 5.0000000000000000e+00 \n') 
    f.write('  nPhotonPackets = %d \n'%nphotons )                    
    f.write('  ranseed = -100 \n')
    f.write('  verbose = T \n')

    f.write('[ramses] \n')      
    if 'RT' in sim_dir:
        f.write('  self_shielding = F \n')
        f.write('  ramses_rt = T \n') 
        f.write('  verbose = T \n') 
        f.write('  use_initial_mass = T \n')  
        f.write('  cosmo = T \n')
        f.write('  use_proper_time = T \n') 
        f.write('  read_rt_variables = T \n')
    else:
        f.write('  self_shielding = F \n')
        f.write('  ramses_rt = F \n') 
        f.write('  verbose = T \n') 
        f.write('  use_initial_mass = T \n')  
        f.write('  cosmo = T \n')
        f.write('  use_proper_time = T \n') 
        f.write('  read_rt_variables = F \n')

    if 'RTCRiMHD+SfFb' in sim_dir:
        f.write('  itemp  = 12 \n')
        f.write('  imetal = 14 \n')
        f.write('  ihii   = 15 \n')
        f.write('  iheii  = 16 \n')
        f.write('  iheiii = 17 \n')

    elif 'RT+SfFb' in sim_dir:
        f.write('  itemp  =  5\n')
        f.write('  imetal =  7\n')
        f.write('  ihii   =  8\n')
        f.write('  iheii  =  9\n')
        f.write('  iheiii =  10\n')

    elif 'RTiMHD+SfFb' in sim_dir:
        f.write('  itemp  = 11 \n')
        f.write('  imetal = 13 \n')
        f.write('  ihii   = 14 \n')
        f.write('  iheii  = 15 \n')
        f.write('  iheiii = 16 \n')
    f.close()




def write_param_PFSM(file,lmin,lmax,nphot,outputFile,l0=1190.,beta=0.):
    f = open(file,'w')
    f.write('[PhotonsFromSourceModel]\n')
    f.write('# input / output parameters \n')
    f.write('outputfile = %s \n'%(outputFile))
    f.write('# computational domain parameters \n')
    f.write('source_type = pointlike \n')
    f.write('source_pos  = 0.5 0.5 0.5 \n')
    f.write('spec_type = PowLaw \n')
    f.write('spec_powlaw_lmin_Ang = %e \n'%(lmin))
    f.write('spec_powlaw_lmax_Ang = %e \n'%(lmax))
    f.write('spec_powlaw_l0_Ang = %e \n'%(l0))
    f.write('spec_powlaw_beta = %e \n'%(beta))
    f.write('# miscelaneous parameters \n')
    f.write('nphot           = %i \n'%(nphot))
    f.write('ranseed         = -100 \n')
    f.write('verbose         = T \n')
    f.close()
    
def write_params_Ra(file, dat_dir, PhotICFile,outputFile, ndir, mock_name, overwrite='F', serial=False):
    f = open(file,'w')
    if serial == True:
        f.write('[RASCAS-serial]\n')
    if serial == False:
        f.write('[RASCAS]\n')
    f.write('DomDumpDir = %s \n'%dat_dir)
    f.write('DomDumpFile  = %s/domain_decomposition_params.dat\n'%dat_dir )
    f.write('PhotonICFile = %s/%s\n'%(dat_dir, PhotICFile))
    f.write('fileout      = %s/%s\n'%(dat_dir, outputFile))
    f.write('nbundle = 10 \n')
    f.write('verbose      = T\n')

    f.write('[worker] \n') 
    f.write('verbose = T \n') 
    f.write('[master]\n') 
    f.write('verbose = T\n')
    f.write('restart = F\n')
    f.write('[gas_composition] \n')
    f.write('# overwrite parameters\n')
    f.write('gas_overwrite       = %s\n'%overwrite)
    f.write('fix_nhi   = 0.0d0\n')
    f.write('fix_vth   = 1.0d4\n')
    f.write('fix_vel   = 0.0d0\n')
    f.write('fix_box_size_cm = 3e21 \n' )
    f.write('f_ion = 1.0000000000000000e-02 \n') 
    f.write('Zref = 5.0000000000000001e-03 \n')
    f.write('verbose = T \n') 
    f.write('[HI] \n')
    f.write('isotropic = F \n') 
    f.write('recoil = T \n')
    f.write('HI_core_skip = T \n')                                                       
    f.write('xcritmax = 10.0000000000000000e+00 \n')                                      
    f.write('[dust] \n')
    f.write('albedo = 3.2000000000000001e-01 \n')                                        
    f.write('g_dust = 7.2999999999999998e-01 \n')  
    f.write('dust_model = SMC \n') 
    f.write('# miscelaneous parameters\n')
    f.write('verbose             = T\n')                                                                                    
    f.write('[mock] \n')
    f.write('nDirections = %d \n'%ndir)                                                                    
    f.write('mock_parameter_file = Lya/mockparams.conf \n')   
    f.write('mock_outputfilename = %s/%s \n'%(dat_dir, mock_name) )
    f.close()

def generate_LOS_a(ndir=100):
    #seed(1)
    LOS_a = [ [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]] # 3 axis
    for i in range(ndir-3):
        theta = random()*np.pi
        phi = random()*2.*np.pi
        LOS_a.append([ np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta) ] )
    LOS_a = np.array(LOS_a)
    return LOS_a

def write_mockparams(file, LOS_a, main_halo_pos, flux_aper = 1.5228222e-02, 
                    spec_npix=1000, spec_aperture = 1.5228222e-02, spec_lmin=1210., spec_lmax=1222.,
                    image_npix=1000, image_side=1.5228222e-02,
                    cube_lbda_npix=0, cube_image_npix=0, cube_lmin=0, cube_lmax=0, cube_side=0 ):
    f = open(file, 'w')
    for LOS in LOS_a:
        f.write('%f %f %f \n'%(LOS[0], LOS[1], LOS[2]) ) 
        f.write(' %f %f %f \n'%(main_halo_pos[0], main_halo_pos[1], main_halo_pos[2] ) ) 
        f.write(' %e \n'%flux_aper)  
        f.write(' %d %f %f %f\n'%(spec_npix, spec_aperture, spec_lmin, spec_lmax) ) 
        f.write(' %d %f \n'%(image_npix, image_side) )    
        f.write(' %d %d %f %f %f \n'%(cube_lbda_npix, cube_image_npix, cube_lmin, cube_lmax, cube_side) )
    f.close()



sim_name=sys.argv[1]
snapnum=sys.argv[2]
print(sim_name, snapnum)
op_idx_a, main_halo_pos_a, main_halo_vel_a, gas_halo_pos_a, gas_halo_vel_a = pu.get_halo_info(sim_name)
lfac, dfac, tfac, redshift, redshiftnum, main_halo_pos, main_halo_vel, gas_halo_pos, gas_halo_vel = pu.get_sim_info(sim_name, snapnum)
sys.path.append('/home/yy503/Desktop/emi_line_eor/simulation')


rsp =  3./lfac/cm2kpc # 6 pkpc region 
nphotons = 10000

sim_dir = '/data/ERCblackholes2/smartin/Dwarf1/%s'%sim_name
dat_dir = '/data/ERCblackholes3/yuxuan/emi_line/%s/output000%s'%(sim_name, snapnum)


os.makedirs('/data/ERCblackholes3/yuxuan/emi_line/%s'%sim_name, exist_ok = True)
os.makedirs(dat_dir, exist_ok = True)

ndir=9
LOS_a = generate_LOS_a(ndir=ndir)
np.savetxt('%s/LOS.txt'%dat_dir, LOS_a) 



write_params_CDD('Lya/CDD.conf', dat_dir, sim_dir, snapnum, gas_halo_pos, rsp, overwrite='F')
write_params_LyaPhotonFromGas('Lya/ppic_col.conf', dat_dir, sim_dir, snapnum, gas_halo_pos, rsp, nphotons, rec='False', col='True')

write_mockparams('Lya/mockparams.conf' , LOS_a, gas_halo_pos, flux_aper = 2*rsp, 
                    spec_npix=1000, spec_aperture = 2*rsp, spec_lmin=1210., spec_lmax=1222.,
                    image_npix=1000, image_side= 2*rsp,
                    cube_lbda_npix=0, cube_image_npix=0, cube_lmin=0, cube_lmax=0, cube_side=0 )
write_params_Ra('Lya/RaS_col.conf', dat_dir, 'col.IC', 'col.res', ndir=ndir, mock_name='col', overwrite='F')

write_params_LyaPhotonFromGas('Lya/ppic_rec.conf', dat_dir, sim_dir, snapnum, gas_halo_pos, rsp, nphotons, rec='True', col='False')
write_params_Ra('Lya/RaS_rec.conf', dat_dir, 'rec.IC', 'rec.res', ndir=ndir, mock_name='rec', overwrite='F')

print('conf file generated')

''' 
# generate PhotICs around 1190.A and run RASCAS
l90 = 1190.42
dl = 0.01
write_param_PFSM('params_PFSM_1190.dat',l90-dl,l90+dl,500000,'PhotICs_1190.dat')
write_params_RaS('params_RaS_1190.dat','PhotICs_1190.dat','result_1190_3e13.dat',boxsize_cm=3e13)

# generate PhotICs around 1193 A and run RASCAS
l93 = 1193.28
dl = 0.01
write_param_PFSM('params_PFSM_1193.dat',l93-dl,l93+dl,500000,'PhotICs_1193.dat')
write_params_RaS('params_RaS_1193.dat','PhotICs_1193.dat','result_1193_3e13.dat',boxsize_cm=3e13)

os.system('./gen_file.sh')
'''


