import numpy as np
import analysis_wamit
import sys

if len(sys.argv)>1:
    dt = sys.argv[1]
else:
    dt = ''

val=dt

while (val != 'm') & (val != 'p') & (val != 'c'):
    val = input('Drift Analysis: [''m'',''p'',''c'']: ')
    print(val != 'm')
    print(val != 'p')
    print(val != 'c')
    
dt = val


# Reading Output parameters

p_out = analysis_wamit.output_params()

[params,axis,vol,cb,cg,rest_coef,name_out] = p_out

[g,ulen, rho, water_depth, water_depth_aux, NBODY] = params

# Reading RAO
[rao, rao_phase, per, inc, dof, arq4d, rao_c] = analysis_wamit.raos(param_out=p_out)

# Reading Wave Forces
[wforce, wforce_phase, arq2d] = analysis_wamit.wave_forces(param_out=p_out)
    
# Reading Drift forces ('m' - Momentum, 'p' - Pressure, 'c' - control surface)




[wdforce, wdforce_phase, arq8d] = analysis_wamit.drift_forces(drift_analysis_type = dt, param_out=p_out)

# Reading Added Mass and Potential Damping
[added_mass, pot_damp, dof1, arq1d, added_mass_matrix, pot_damp_matrix]= analysis_wamit.added_mass_pot_damping(param_out=p_out)

dof=np.reshape(dof,(NBODY,6))


for ii in range(NBODY):
    [xbody,ybody,zbody,phibody] = axis[ii]
    [xvol,yvol,zvol] = vol[ii]
    [xb,yb,zb] = cb[ii]
    [xg,yg,zg] = cg[ii]
    [c33,c34,c35,c44,c45,c46,c55,c56] = rest_coef[ii]

    # Reading Force file 
    [mass,damp,rest_coef_ext] = analysis_wamit.read_frc()
    
    # Names of Degrees of Freedom
    naemsDof_aux = ['SURGE','SWAY','HEAVE','ROLL','PITCH','YAW']
    namesDof=[]
    for nn in range(NBODY):
        for nd in naemsDof_aux:
            namesDof.append(nd)
        
    # Write parameters in ship.wnf
    name_wnf = 'ship{:d}.wnf'.format(ii + 1)
    wnf = open(name_wnf,'w')
    wnf.write('%OUT_FILE\n')
    wnf.write('"'+name_out+'"\n\n')
    wnf.write('%NUMBER_OF_BODIES\n')
    wnf.write(str(NBODY) + '\n\n')
    wnf.write('%BODY\n')
    wnf.write(str(ii+1) + '\n\n')
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
    wnf.write('{:.1f}'.format(mass[ii][0][0])+'\n\n')
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
    
    for i in mass[ii]:
        wnf.write('    '.join('{:.5e}'.format(x) for x in i)+'\n')
        
    wnf.write('\n')
    wnf.write('%EXTERNAL_DAMPING\n')
    
    for i in damp[ii]:
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
    for j in range(dof[ii].size):
        wnf.write('%RAO_' + namesDof[j] + '_AMPLITUDE\n')
        wnf.write(tx_inc + '\n')
        for i in range(per.size):
            tx_rao = '{:.6e}'.format(2*np.pi/per[i])
            for x in rao[ int(dof[ii][j])-1 ][i]:
                aux = '{:.6e}'.format(x)
                n_carac = len(aux) 
                tx_rao = tx_rao + (16-n_carac)*' ' + aux 
            wnf.write(tx_rao + '\n')
        wnf.write('\n')
        wnf.write('%RAO_' + namesDof[j] + '_PHASE\n')
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
    for j in range(dof[ii].size):
        wnf.write('%EXCITING_WAVE_FORCE_' + namesDof[j] + '_AMPLITUDE\n')
        wnf.write(tx_inc + '\n')
        
        for i in range(per.size):
            tx_rao = '{:.6e}'.format(2*np.pi/per[i])
            for x in wforce[int(dof[ii][j])-1][i]:
                aux = '{:.6e}'.format(x)
                n_carac = len(aux) 
                tx_rao = tx_rao + (16-n_carac)*' ' + aux 
            wnf.write(tx_rao + '\n')
        wnf.write('\n')
        wnf.write('%EXCITING_WAVE_FORCE_' + namesDof[j] + '_PHASE\n')
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
    for j in range(dof[ii].size):
        wnf.write('%SLOW_DRIFT_FORCE_' + namesDof[j] + '_AMPLITUDE\n')
        wnf.write(tx_inc + '\n')
        
        for i in range(per.size):
            tx_rao = '{:.6e}'.format(2*np.pi/per[i])
            for x in wdforce[int(dof[ii][j])-1][i]:
                aux = '{:.6e}'.format(x)
                n_carac = len(aux) 
                tx_rao = tx_rao + (16-n_carac)*' ' + aux 
            wnf.write(tx_rao + '\n')
            
        wnf.write('\n')
        wnf.write('%SLOW_DRIFT_FORCE_' + namesDof[j] + '_PHASE\n')
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
    sz = int(dof1.shape[0]/NBODY)
    dof1_aux = dof1[ii*sz:ii*sz+sz]
    
    for x in dof1_aux:
        str_aux = 'A(' + '{:d}'.format(int(x[0])) + ',' '{:d}'.format(int(x[1])) + ')'
        blank_spaces = 16 - len(str_aux)
        tx_ad = tx_ad + blank_spaces*' ' + str_aux
    wnf.write('\n')
    wnf.write('%ADDED_MASS\n')
    wnf.write(tx_ad + '\n')
    
    for i in range(per.size):
        tx_rao = '{:.6e}'.format(2*np.pi/per[i])
        for x in added_mass[i][ii*sz:ii*sz+sz]:
            aux = '{:.6e}'.format(x)
            n_carac = len(aux) 
            tx_rao = tx_rao + (16-n_carac)*' ' + aux 
        wnf.write(tx_rao + '\n')
        
    #POTENTIAL_DAMPING
    tx_pd = 12*' '
    
    for x in dof1_aux:
        str_aux = 'B(' + '{:d}'.format(int(x[0])) + ',' '{:d}'.format(int(x[1])) + ')'
        blank_spaces = 16 - len(str_aux)
        tx_pd = tx_pd + blank_spaces*' ' + str_aux
    
    wnf.write('\n')
    wnf.write('%POTENTIAL_DAMPING_COEFFICIENTS\n')
    wnf.write(tx_pd + '\n')
    
    for i in range(per.size):
        tx_rao = '{:.6e}'.format(2*np.pi/per[i])
        for x in pot_damp[i][ii*sz:ii*sz+sz]:
            aux = '{:.6e}'.format(x)
            n_carac = len(aux) 
            tx_rao = tx_rao + (16-n_carac)*' ' + aux 
        wnf.write(tx_rao + '\n')
    
    #REST_Coefs
    wnf.write('\n')
    wnf.write('%REST_COEFS\n')
    wnf.write('{:.6e}'.format(rest_coef[ii][0]))
    
    for x in rest_coef[ii][1:]:
        aux = '{:.6e}'.format(x)
        n_carac = len(aux)
        wnf.write((16-n_carac)*' ' + aux)
        
    wnf.write('\n')
    wnf.write('\n')
    wnf.write('%END\n')    
    wnf.close()

print('')
print(' * Generated WNF file(s)' )
print('')
