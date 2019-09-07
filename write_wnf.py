import numpy as np
import analysis_wamit

# Reading Output parameters

[params,axis,vol,cb,cg,rest_coef,nome_out] = analysis_wamit.output_params()

[g,ulen,rho,water_depth,water_depth_aux] = params
[xbody,ybody,zbody,phibody] = axis
[xvol,yvol,zvol] = vol
[xb,yb,zb] = cb
[xg,yg,zg] = cg
[c33,c34,c35,c44,c45,c46,c55,c56] = rest_coef

# Reading Force file 
[mass,damp,rest_coef_ext] = analysis_wamit.read_frc()

nomesDof = ['SURGE','SWAY','HEAVE','ROLL','PITCH','YAW']
    
# Reading RAO
[rao, rao_phase, per, inc, dof, arq4d, rao_c] = analysis_wamit.raos()

# Reading Wave Forces
[wforce, wforce_phase, arq2d] = analysis_wamit.wave_forces()
    
# Reading Drift forces (Momentum)
[wdforce, wdforce_phase, arq8d] = analysis_wamit.drift_forces_momentum()

# Reading Added Mass and Potential Damping
[added_mass, pot_damp, dof1, arq1d, added_mass_matrix, pot_damp_matrix] = analysis_wamit.added_mass_pot_damping()

# Write parameters in ship.wnf
wnf = open('ship.wnf','w')
wnf.write('%OUT_FILE\n')
wnf.write('"'+nome_out+'"\n\n')
wnf.write('%NUMBER_OF_BODIES\n')
wnf.write('1\n\n')
wnf.write('%BODY\n')
wnf.write('1\n\n')
wnf.write('%SCALE\n')
wnf.write(str(ulen)+'\n\n')
wnf.write('%WATER_DENSITY\n')
wnf.write(str(rho)+'\n\n')
wnf.write('%GRAVITY\n')
wnf.write(str(g)+'\n\n')
wnf.write('%DISPLACEMENT_VOLUME\n')
#wnf.write('{:.1f}'.format(np.mean([xvol,yvol,zvol]))+'\n\n') #mean volume
wnf.write('{:.1f}'.format(xvol)+'\n\n')
wnf.write('%DISPLACEMENT_WEIGHT\n')
wnf.write('{:.1f}'.format(mass[0][0])+'\n\n')
wnf.write('%FLOAT_CENTER\n')
wnf.write('{:.6f}'.format(xb)+' '+'{:.6f}'.format(yb)+' '+'{:.6f}'.format(zb)+'\n\n')
wnf.write('%WAMIT_AXIS\n')
wnf.write('{:.6f}'.format(xbody)+' '+'{:.6f}'.format(ybody)+' '+'{:.6f}'.format(zbody)+' '+'{:.6f}'.format(phibody)+'\n\n')
wnf.write('%CENTER_OF_GRAVITY\n')
wnf.write('{:.6f}'.format(xg)+' '+'{:.6f}'.format(yg)+' '+'{:.6f}'.format(zg)+'\n\n')
wnf.write('%DEPTH\n')

if water_depth_aux==1:
    wnf.write(water_depth+'\n\n')
else:
    wnf.write('{:f}'.format(water_depth)+'\n\n')
    
wnf.write('%GLOBAL_MASS\n')

for i in mass:
    wnf.write('    '.join('{:.5e}'.format(x) for x in i)+'\n')
    
wnf.write('\n')
wnf.write('%EXTERNAL_DAMPING\n')

for i in damp:
    wnf.write('    '.join('{:.5e}'.format(x) for x in i)+'\n')
    
wnf.write('\n')
wnf.write('%NUMBER_FREQUENCIES\n')
wnf.write(str(per.size)+'\n\n')
wnf.write('%NUMBER_HEADINGS\n')
wnf.write(str(inc.size)+'\n\n')
tx_inc = 12*' '

for x in inc:
    aux = '{:.0f}'.format(x)
    n_carac = len(aux)
    tx_inc = tx_inc + (16-n_carac)*' ' + aux
    
#RAOs
for j in range(dof.size):
    wnf.write('%RAO_' + nomesDof[j] + '_AMPLITUDE\n')
    wnf.write(tx_inc + '\n')
    for i in range(per.size):
        tx_rao = '{:.6e}'.format(2*np.pi/per[i])
        for x in rao[j][i]:
            aux = '{:.6e}'.format(x)
            n_carac = len(aux) 
            tx_rao = tx_rao + (16-n_carac)*' ' + aux 
        wnf.write(tx_rao + '\n')
    wnf.write('\n')
    wnf.write('%RAO_' + nomesDof[j] + '_PHASE\n')
    wnf.write(tx_inc + '\n')
    
    for i in range(per.size):
        tx_rao = '{:.6e}'.format(2*np.pi/per[i])
        for x in rao_phase[j][i]:
            aux = '{:.6e}'.format(x)
            n_carac = len(aux) 
            tx_rao = tx_rao + (16-n_carac)*' ' + aux 
        wnf.write(tx_rao + '\n')
    wnf.write('\n')

#EXCITING_WAVE_FORCE
for j in range(dof.size):
    wnf.write('%EXCITING_WAVE_FORCE_' + nomesDof[j] + '_AMPLITUDE\n')
    wnf.write(tx_inc + '\n')
    
    for i in range(per.size):
        tx_rao = '{:.6e}'.format(2*np.pi/per[i])
        for x in wforce[j][i]:
            aux = '{:.6e}'.format(x)
            n_carac = len(aux) 
            tx_rao = tx_rao + (16-n_carac)*' ' + aux 
        wnf.write(tx_rao + '\n')
    wnf.write('\n')
    wnf.write('%EXCITING_WAVE_FORCE_' + nomesDof[j] + '_PHASE\n')
    wnf.write(tx_inc + '\n')
    
    for i in range(per.size):
        tx_rao = '{:.6e}'.format(2*np.pi/per[i])
        for x in wforce_phase[j][i]:
            aux = '{:.6e}'.format(x)
            n_carac = len(aux) 
            tx_rao = tx_rao + (16-n_carac)*' ' + aux 
        wnf.write(tx_rao + '\n')
    wnf.write('\n')

#SLOW_DRIFT_FORCE
for j in range(dof.size):
    wnf.write('%SLOW_DRIFT_FORCE_' + nomesDof[j] + '_AMPLITUDE\n')
    wnf.write(tx_inc + '\n')
    
    for i in range(per.size):
        tx_rao = '{:.6e}'.format(2*np.pi/per[i])
        for x in wdforce[j][i]:
            aux = '{:.6e}'.format(x)
            n_carac = len(aux) 
            tx_rao = tx_rao + (16-n_carac)*' ' + aux 
        wnf.write(tx_rao + '\n')
        
    wnf.write('\n')
    wnf.write('%SLOW_DRIFT_FORCE_' + nomesDof[j] + '_PHASE\n')
    wnf.write(tx_inc + '\n')
    
    for i in range(per.size):
        tx_rao = '{:.6e}'.format(2*np.pi/per[i])
        for x in wdforce_phase[j][i]:
            aux = '{:.6e}'.format(x)
            n_carac = len(aux) 
            tx_rao = tx_rao + (16-n_carac)*' ' + aux 
        wnf.write(tx_rao + '\n')
    wnf.write('\n')
    
#ADDED_MASS
tx_ad = 12*' '

for x in dof1:
    tx_ad = tx_ad + 10*' ' + 'A(' + '{:d}'.format(int(x[0])) + ',' '{:d}'.format(int(x[1])) + ')'
wnf.write('\n')
wnf.write('%ADDED_MASS\n')
wnf.write(tx_ad + '\n')

for i in range(per.size):
    tx_rao = '{:.6e}'.format(2*np.pi/per[i])
    for x in added_mass[i]:
        aux = '{:.6e}'.format(x)
        n_carac = len(aux) 
        tx_rao = tx_rao + (16-n_carac)*' ' + aux 
    wnf.write(tx_rao + '\n')
    
#POTENTIAL_DAMPING
tx_pd = 12*' '

for x in dof1:
    tx_pd = tx_pd + 10*' ' + 'B(' + '{:d}'.format(int(x[0])) + ',' '{:d}'.format(int(x[1])) + ')'

wnf.write('\n')
wnf.write('%POTENTIAL_DAMPING_COEFFICIENTS\n')
wnf.write(tx_pd + '\n')

for i in range(per.size):
    tx_rao = '{:.6e}'.format(2*np.pi/per[i])
    for x in pot_damp[i]:
        aux = '{:.6e}'.format(x)
        n_carac = len(aux) 
        tx_rao = tx_rao + (16-n_carac)*' ' + aux 
    wnf.write(tx_rao + '\n')

#REST_Coefs
wnf.write('\n')
wnf.write('%REST_COEFS\n')

wnf.write('{:.6e}'.format(rest_coef[0]))

for x in rest_coef[1:]:
    aux = '{:.6e}'.format(x)
    n_carac = len(aux)
    wnf.write((16-n_carac)*' ' + aux)
    
wnf.write('\n')
wnf.write('\n')
wnf.write('%END\n')

wnf.close()