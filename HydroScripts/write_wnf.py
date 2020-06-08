import numpy as np

def wnf_sim(p_out, a, w, r, d, dt):
    
    [params, axis, vol, cb, cg, rest_coef, nome_out, GMt, GMl, M, Bvisc, C, Cext] = p_out

    [g,ulen, rho, water_depth, water_depth_aux, NBODY] = params

    if NBODY > 1 and dt == 'm':
        raise NameError('For multiple bodies the drift analysis must be Pressure or Control Surface!')

    [added_mass, pot_damp, dof1, arq1d, added_mass_matrix, pot_damp_matrix] = a

    [rao, rao_phase, per, inc, dof, arq4d, rao_c] = r

    [wforce, wforce_phase, arq2d] = w

    [wdforce, wdforce_phase, arq8d] = d

    dof=np.reshape(dof, (int(np.max(dof)/6), 6))

    # [M,Bvisc,Cext] = analysis_wamit.read_mmx()
    
    for ii in range(len(dof)):
        [xbody,ybody,zbody,phibody] = axis[ii]
        [xvol,yvol,zvol] = vol[ii]
        [xb,yb,zb] = cb[ii]
        [xg,yg,zg] = cg[ii]
        [c33,c34,c35,c44,c45,c46,c55,c56] = rest_coef[ii]

        # Reading Force file 
        # [mass,damp,rest_coef_ext] = analysis_wamit.read_frc()
        pi = ii*6
        pf = ii*6+6
        mass = M[pi:pf, pi:pf]
        damp = Bvisc[pi:pf, pi:pf]
        
        # Names of Degrees of Freedom
        namesDof_aux = ['SURGE','SWAY','HEAVE','ROLL','PITCH','YAW']
        namesDof=[]
        for nn in range(NBODY):
            for nd in namesDof_aux:
                namesDof.append(nd)

        sufix = {'m':'_mom', 'p':'_press', 'c': '_csf'}

        # Write parameters in ship.wnf
        name_wnf = 'ship{:d}'.format(ii + 1)  + '_sim' + sufix[dt] + '.wnf'
        wnf = open(name_wnf,'w')
        wnf.write('%OUT_FILE\n')
        wnf.write('"' + nome_out + '"\n\n')
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
            n_dec = 0
            z = x % 1
            while z != 0:
                n_dec = n_dec + 1
                z = z*10 % 1
                # print(n_dec)
                # print(z)

            aux = '{:.' + '{:d}'.format(n_dec) +'f}'
            # print(aux)
            aux = aux.format(x)
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
                for x in rao_phase[ int(dof[ii][j])-1 ][i]:
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
                for x in wforce[ int(dof[ii][j])-1 ][i]:
                    aux = '{:.6e}'.format(x)
                    n_carac = len(aux) 
                    tx_rao = tx_rao + (16-n_carac)*' ' + aux 
                wnf.write(tx_rao + '\n')
            wnf.write('\n')
            wnf.write('%EXCITING_WAVE_FORCE_' + namesDof[j] + '_PHASE\n')
            wnf.write(tx_inc + '\n')
            
            for i in range(per.size):
                tx_rao = '{:.6e}'.format(2*np.pi/per[i])
                for x in wforce_phase[ int(dof[ii][j])-1 ][i]:
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
                for x in wdforce[ int(dof[ii][j])-1 ][i]:
                    aux = '{:.6e}'.format(x)
                    n_carac = len(aux) 
                    tx_rao = tx_rao + (16-n_carac)*' ' + aux 
                wnf.write(tx_rao + '\n')
                
            wnf.write('\n')
            wnf.write('%SLOW_DRIFT_FORCE_' + namesDof[j] + '_PHASE\n')
            wnf.write(tx_inc + '\n')
            
            for i in range(per.size):
                tx_rao = '{:.6e}'.format(2*np.pi/per[i])
                for x in wdforce_phase[ int(dof[ii][j])-1 ][i]:
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
        print(' * Generated ' + name_wnf + ' file(s)' )

def wnf_tpn(p_out, a, w, r, d, dt):
   
    [params, axis, vol, cb, cg, rest_coef, nome_out, GMt, GMl, M, Bvisc, C, Cext] = p_out

    [g,ulen, rho, water_depth, water_depth_aux, NBODY] = params

    if NBODY > 1 and dt == 'm':
        raise NameError('For multiple bodies the drift analysis must be Pressure or Control Surface!')

    [added_mass, pot_damp, dof1, arq1d, added_mass_matrix, pot_damp_matrix] = a

    [rao, rao_phase, per, inc, dof, arq4d, rao_c] = r

    [wforce, wforce_phase, arq2d] = w

    [wdforce, wdforce_phase, arq8d] = d

    dof=np.reshape(dof, (int(np.max(dof)/6), 6))

    for ii in range(len(dof)):
        [xbody,ybody,zbody,phibody] = axis[ii]
        [xvol,yvol,zvol] = vol[ii]
        [xb,yb,zb] = cb[ii]
        [xg,yg,zg] = cg[ii]
        [c33,c34,c35,c44,c45,c46,c55,c56] = rest_coef[ii]

        # Reading Force file 
        # [mass,damp,rest_coef_ext] = analysis_wamit.read_frc()
        pi = ii*6
        pf = ii*6+6
        mass = M[pi:pf, pi:pf]
        damp = Bvisc[pi:pf, pi:pf]
        
        # Names of Degrees of Freedom
        namesDof_aux = ['SURGE','SWAY','HEAVE','ROLL','PITCH','YAW']
        namesDof=[]
        for nn in range(NBODY):
            for nd in namesDof_aux:
                namesDof.append(nd)

        sufix = {'m':'_mom', 'p':'_press', 'c': '_csf'}

        # Write parameters in ship.wnf
        name_wnf = 'ship{:d}'.format(ii + 1)  + '_tpn' + sufix[dt] + '.wnf'
        wnf = open(name_wnf,'w')
        wnf.write('%VERSION\n')
        wnf.write('2.1.4\n\n')
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
        wnf.write('%WAMIT_AXIS\n')
        wnf.write('{:.4f}'.format(xbody)+' '+'{:.4f}'.format(ybody)+' '+'{:.4f}'.format(zbody)+' '+'{:.1f}'.format(phibody)+'\n\n')
        wnf.write('%PHIBODY\n')
        wnf.write('{:.1f}'.format(phibody)+'\n\n') 
        wnf.write('%SYMMETRY\n')
        wnf.write('0\n\n')
        wnf.write('%CENTER_OF_BUOYANCY\n')
        wnf.write('{:.6f}'.format(xb)+' '+'{:.6f}'.format(yb)+' '+'{:.6f}'.format(zb)+'\n\n')
        wnf.write('%CENTER_OF_GRAVITY\n')
        wnf.write('{:.6f}'.format(xg)+' '+'{:.6f}'.format(yg)+' '+'{:.6f}'.format(zg)+'\n\n')
        wnf.write('%DEPTH\n')
        # wnf.write('\n\n')
        if water_depth_aux==1:
            wnf.write(water_depth+'\n\n')
        else:
            wnf.write('{:f}'.format(water_depth)+'\n\n')

        wnf.write('%VOLUME\n')
        wnf.write('{:.0f}'.format(np.mean([xvol,yvol,zvol]))+'\n\n') #mean volume
        
        wnf.write('%VOLUMEX\n')
        wnf.write('{:.0f}'.format(xvol)+'\n\n')
        wnf.write('%VOLUMEY\n')
        wnf.write('{:.0f}'.format(yvol)+'\n\n')
        wnf.write('%VOLUMEZ\n')
        wnf.write('{:.0f}'.format(zvol)+'\n\n')

        # wnf.write('%DISPLACEMENT_WEIGHT\n')
        # wnf.write('{:.1f}'.format(mass[ii][0][0])+'\n\n')
        
        wnf.write('%MASS_AND_INERTIA\n')

        for i in mass:
            wnf.write('    '.join('{:.5e}'.format(x) for x in i)+'\n')
            
        wnf.write('\n')

        wnf.write('%EXTERNAL_GLOBAL_DAMPING\n')
        for i in damp:
            wnf.write('    '.join('{:.5e}'.format(x) for x in i)+'\n')
            
        wnf.write('\n')
        wnf.write('%NUMBER_FREQUENCIES\n')
        wnf.write(str(per.size)+'\n')
        wnf.write('%NUMBER_HEADINGS\n')
        wnf.write(str(inc.size)+'\n')
        tx_inc = 12*' '
        
        for x in inc:
            n_dec = 0
            z = x % 1
            while z != 0:
                n_dec = n_dec + 1
                z = z*10 % 1
                # print(n_dec)
                # print(z)

            aux = '{:.' + '{:d}'.format(n_dec) +'f}'
            # print(aux)
            aux = aux.format(x)
            n_carac = len(aux)
            tx_inc = tx_inc + (16-n_carac)*' ' + aux
            
        #RAOs
        for j in range(dof[ii].size):
            wnf.write('%RAO_' + namesDof[j] + '_MOD\n')
            wnf.write(tx_inc + '\n')
            for i in range(per.size):
                tx_rao = '{:.5e}'.format(2*np.pi/per[i])
                for x in rao[ int(dof[ii][j])-1 ][i]:
                    aux = '{:.5e}'.format(x)
                    n_carac = len(aux) 
                    tx_rao = tx_rao + (16-n_carac)*' ' + aux 
                wnf.write(tx_rao + '\n')
            wnf.write('\n')
            wnf.write('%RAO_' + namesDof[j] + '_PHA\n')
            wnf.write(tx_inc + '\n')
            
            for i in range(per.size):
                tx_rao = '{:.5e}'.format(2*np.pi/per[i])
                for x in rao_phase[ int(dof[ii][j])-1 ][i]:
                    aux = '{:.5e}'.format(x)
                    n_carac = len(aux) 
                    tx_rao = tx_rao + (16-n_carac)*' ' + aux 
                wnf.write(tx_rao + '\n')
            wnf.write('\n')
        
        #EXCITING_WAVE_FORCE
        for j in range(dof[ii].size):
            wnf.write('%HASKIND_EXCITING_FORCES_' + namesDof[j] + '_MOD\n')
            wnf.write(tx_inc + '\n')
            
            for i in range(per.size):
                tx_rao = '{:.5e}'.format(2*np.pi/per[i])
                for x in wforce[ int(dof[ii][j])-1 ][i]:
                    aux = '{:.5e}'.format(x)
                    n_carac = len(aux) 
                    tx_rao = tx_rao + (16-n_carac)*' ' + aux 
                wnf.write(tx_rao + '\n')
            wnf.write('\n')
            wnf.write('%HASKIND_EXCITING_FORCES_' + namesDof[j] + '_PHA\n')
            wnf.write(tx_inc + '\n')
            
            for i in range(per.size):
                tx_rao = '{:.5e}'.format(2*np.pi/per[i])
                for x in wforce_phase[ int(dof[ii][j])-1 ][i]:
                    aux = '{:.5e}'.format(x)
                    n_carac = len(aux) 
                    tx_rao = tx_rao + (16-n_carac)*' ' + aux 
                wnf.write(tx_rao + '\n')
            wnf.write('\n')
        
        #SLOW_DRIFT_FORCE
        for j in range(dof[ii].size):
            wnf.write('%DRIFT_FORCES_' + namesDof[j] + '_MOD\n')
            wnf.write(tx_inc + '\n')
            
            for i in range(per.size):
                tx_rao = '{:.5e}'.format(2*np.pi/per[i])
                for x in wdforce[ int(dof[ii][j])-1 ][i]:
                    aux = '{:.5e}'.format(x)
                    n_carac = len(aux) 
                    tx_rao = tx_rao + (16-n_carac)*' ' + aux 
                wnf.write(tx_rao + '\n')
                
            wnf.write('\n')
            wnf.write('%DRIFT_FORCES_' + namesDof[j] + '_PHA\n')
            wnf.write(tx_inc + '\n')
            
            for i in range(per.size):
                tx_rao = '{:.5e}'.format(2*np.pi/per[i])
                for x in wdforce_phase[ int(dof[ii][j])-1 ][i]:
                    aux = '{:.5e}'.format(x)
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
            tx_rao = '{:.5e}'.format(2*np.pi/per[i])
            for x in added_mass[i][ii*sz:ii*sz+sz]:
                aux = '{:.5e}'.format(x)
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
        wnf.write('%POTENTIAL_DAMPING\n')
        wnf.write(tx_pd + '\n')
        
        for i in range(per.size):
            tx_rao = '{:.5e}'.format(2*np.pi/per[i])
            for x in pot_damp[i][ii*sz:ii*sz+sz]:
                aux = '{:.5e}'.format(x)
                n_carac = len(aux) 
                tx_rao = tx_rao + (16-n_carac)*' ' + aux 
            wnf.write(tx_rao + '\n')
        
        #REST_Coefs
        wnf.write('\n')
        wnf.write('%HYDROSTATIC_RESTORING\n')
        wnf.write('{:.5e}'.format(rest_coef[ii][0]))
        
        for x in rest_coef[ii][1:]:
            aux = '{:.5e}'.format(x)
            n_carac = len(aux)
            wnf.write((16-n_carac)*' ' + aux)
            
        wnf.write('\n')
        wnf.write('\n')
        wnf.write('%END\n')    
        wnf.close()

        print('')
        print(' * Generated ' + name_wnf + ' file(s)' )
        