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
import plot_util as pu

cm2pc = (u.cm).to(u.pc)
cm2kpc = (u.cm).to(u.kpc)
s2yr = (u.s).to(u.yr)
s2Gyr = (u.s).to(u.Gyr)
g2Msun = (u.g).to(u.M_sun)


dat_dir = '/data/ERCblackholes2/smartin/Dwarf1'

pp_dir = '/data/ERCblackholes3/yuxuan/dwarf_data/post_processing'
prj_dir = '/home/yy503/Desktop/emi_line_eor'


redshift_a = []
time_a = []

# starting point
sim_name_a1 = ['HD+SfFb', 'RT+SfFb', 'RTiMHD+SfFb', 'RTCRiMHD+SfFb']
op_idx_final_a1 = ['00056','00058','00053','00066'] 
center_final_a1  = np.array([ [0.5218194249249059, 0.4734236787303123, 0.48469643197929374], #HD+SfFb_00056
  [0.5218773023297338, 0.4734778679807518, 0.4844564553481938],   #RT+SfFb_00058
  [0.5219084718820158, 0.47347109473437143, 0.48448011400515445], #RTiMHD+SfFb_00053   
  [0.5216487828425217, 0.4735734139372305, 0.4846249604429559] ]) #RTCRiMHD+SfFb_00066



def extract_main_halo_pos(sim_name, parttype):
	# parttype: darkmatter, gasystars
	redshift_a2 = []
	op_idx_a2 = []
	main_halo_pos_a2 = []
	main_halo_vel_a2 = []

	#isim = sim_name_a1.index(sim_name)
	#op_idx_final = int( op_idx_final_a1[isim] )
	op_idx_final=54
	op_idx_start=30
    # save op_idx_a and main halo position
	op_idx_a = np.loadtxt('%s/%s/op_index.txt'%(pp_dir, sim_name) )
	main_halo_pos_a = np.loadtxt('%s/%s/main_halo_pos.txt'%(pp_dir, sim_name) )

	for i, op_idx in enumerate( np.arange(op_idx_final, op_idx_start, -1) ):
		op_idx = int(op_idx); print(op_idx)
		op_idx_a2.append(op_idx)

		# setting halo center and halo radius. tracking DM halo backward in time; tracking gas halo with calculated DM halo.  
		if parttype=='darkmatter':
			if i==0:
				xzoom=[0.5, 0.5, 0.5]; rzoom=0.1
			elif i>0:
				xzoom=main_halo_pos_a2[i-1]; rzoom=1e-2
		elif parttype=='gasystars':
			idx = np.where(op_idx_a==int(op_idx) )[0][0]
			main_halo_pos = main_halo_pos_a[idx]
			xzoom = main_halo_pos
			rzoom=1e-2

		

		all_halos = hmr.read_halos(ID=op_idx,read_route=dat_dir+'/'+sim_name,
									parttype=parttype, verbose=1, xzoom=xzoom, rzoom=rzoom )

		halo_mass = []
		for i in range( len(all_halos) ):
			halo_mass.append(all_halos[i].halo_mass )
		main_halo_pos_a2.append( all_halos[ np.argmax( halo_mass) ].halo_pos )
		main_halo_vel_a2.append( all_halos[ np.argmax( halo_mass) ].halo_vel )


	#zipped = zip(op_idx_a2, main_halo_pos_a2, main_halo_vel_a2 )
	# store the position for most massive halo
	#np.savetxt('%s/%s/main_halo.csv'%(pp_dir, sim_name), zipped)
	#np.savetxt(main_halo_vel_a2, '%s/%s/main_halo_vel.txt'% )
	np.savetxt('%s/%s/op_index.txt'%(pp_dir, sim_name), op_idx_a2, )

	if parttype=='darkmatter':
		np.savetxt('%s/%s/main_halo_pos.txt'%(pp_dir, sim_name), main_halo_pos_a2,  )
		np.savetxt('%s/%s/main_halo_vel.txt'%(pp_dir, sim_name), main_halo_vel_a2,  )

	if parttype=='gasystars':
		np.savetxt('%s/%s/gas_halo_pos.txt'%(pp_dir, sim_name), main_halo_pos_a2,  )
		np.savetxt('%s/%s/gas_halo_vel.txt'%(pp_dir, sim_name), main_halo_vel_a2,  )
	return redshift_a2, op_idx_a2, main_halo_pos_a2, main_halo_vel_a2

sim_name = sys.argv[1]

extract_main_halo_pos(sim_name, parttype='gasystars')




