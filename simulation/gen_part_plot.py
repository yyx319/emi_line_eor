# sunset stellar image 
 

# rotation curve
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

locals()['v_c_%d'%i ] = v_c


plt.plot(r_bin[1:], v_c, label = sim_name)
plt.plot(r_bin[1:], v_c_g, '--', label = 'gas')
plt.plot(r_bin[1:], v_c_p, '.', label = 'dm+star' )
plt.xlabel(r'r [kpc]')
plt.ylabel(r'V_c [km/s]')
plt.legend()
plt.savefig('%s/figure/%s/%s_rotation_curve.png'%(prj_dir, sim_name, sim_name)  )
plt.clf()
