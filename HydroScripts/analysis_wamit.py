import numpy as np
import re
import math
from plot_curves import plot_curves
from tabulate import tabulate

def read_ULEN():
    nome_out = 'force.out'
    arq_out_aux = open(nome_out,'r')
    arq_out = arq_out_aux.readlines()
    arq_out_aux.close()
    padrao = "[+-]?\\d+(?:\\.\\d+)?(?:[eE][+-]?\\d+)?"
        
    for x in arq_out:   
        if ('Gravity:' in x)==True:
            [g,ulen] = [float(i) for i in re.findall(padrao,x)]
    return ulen

def verify_NBODY():
    nome_out = 'force.out'
    arq_out_aux = open(nome_out,'r')
    arq_out = arq_out_aux.readlines()
    arq_out_aux.close()
    padrao = "[+-]?\\d+(?:\\.\\d+)?(?:[eE][+-]?\\d+)?"
    
    NBODY_Test = []
    NBODY = 1
    
    for x in arq_out: 
        if ('Body number: N='  in x)==True:
            aux = re.findall(padrao,x)
            #print(int(aux[0]))
            NBODY_Test.append(int(aux[0]))
            #print('Achou N = {:d}'.format(int(aux[0])))
    
    if len(NBODY_Test) > 0:
        NBODY = max(NBODY_Test)
    print('')       
    print('NBODY={:d}'.format(NBODY))
    return NBODY
    
def output_params(): 
    # Reading running output parameters
    nome_out = 'force.out'
    arq_out_aux = open(nome_out,'r')
    arq_out = arq_out_aux.readlines()
    arq_out_aux.close()
    padrao = "[+-]?\\d+(?:\\.\\d+)?(?:[eE][+-]?\\d+)?"
    
    NBODY = verify_NBODY()
    
    params_aux = []
    axis = []
    vol = []
    cb = []
    cg = []
    rest_coef = []           
    for x in arq_out:   
        if ('Gravity:' in x)==True:
            [g,ulen] = [float(i) for i in re.findall(padrao,x)]
        if ('Water depth:' in x)==True:
            if len(re.findall(padrao,x))==1:
                [rho] = [float(i) for i in re.findall(padrao,x)]
                water_depth='infinite'
                water_depth_aux=1
            else:
                [water_depth,rho] = [float(i) for i in re.findall(padrao,x)]
                water_depth_aux=0
            params_aux.append([g, ulen, rho, water_depth, water_depth_aux, NBODY])
        if ('XBODY =' in x)==True:
            [xbody,ybody,zbody,phibody] = [float(i) for i in re.findall(padrao,x)]
            axis.append([xbody,ybody,zbody,phibody])
        if ('Volumes (VOLX,VOLY,VOLZ):' in x)==True:
            [xvol,yvol,zvol] = [float(i) for i in re.findall(padrao,x)]
            vol.append([xvol,yvol,zvol])
        if ('Center of Buoyancy (Xb,Yb,Zb):' in x)==True:
            [xb,yb,zb] = [float(i) for i in re.findall(padrao,x)]
            cb.append([xb,yb,zb])
        if ('C(3,3),C(3,4),C(3,5):' in x)==True:
            [_,_,_,_,_,_,c33,c34,c35] = [float(i) for i in re.findall(padrao,x)]
        if ('C(4,4),C(4,5),C(4,6):' in x)==True:
            [_,_,_,_,_,_,c44,c45,c46] = [float(i) for i in re.findall(padrao,x)]
        if ('C(5,5),C(5,6):' in x)==True:
            [_,_,_,_,c55,c56] = [float(i) for i in re.findall(padrao,x)]
            rest_coef.append([c33*g*rho*(ulen**2),c34*g*rho*(ulen**3),c35*g*rho*(ulen**3),
                              c44*g*rho*(ulen**4),c45*g*rho*(ulen**4),c46*g*rho*(ulen**4),
                              c55*g*rho*(ulen**4),c56*g*rho*(ulen**4)])
        if ('Center of Gravity  (Xg,Yg,Zg):' in x)==True:
            [xg,yg,zg] = [float(i) for i in re.findall(padrao,x)]
            cg.append([xg,yg,zg])
                     
    params = params_aux[0]
    axis = axis[0:NBODY]
    vol = vol[0:NBODY]
    cb = cb[0:NBODY]
    cg = cg[0:NBODY]
    rest_coef = rest_coef[0:NBODY]
    
    dof_rest_coef =[[3,3],[3,4],[3,5],[4,4],[4,5],[4,6],[5,5],[5,6]]

    # frc_out = read_frc() #substitute by read_mmx after
    [M, Bvisc, Cext] = read_mmx()
    C = []
    # Cext = []
    # M = []
    Mass = []
    # Bvisc = []
    GMt = []
    GMl = []
    for ii in range(NBODY):
        C.append(np.zeros((6,6)))
        cont=0
        for x in dof_rest_coef:
            C[ii][x[0]-1,x[1]-1] = rest_coef[ii][cont]
            cont+=1
        
        # M.append(frc_out[0][ii])
        # Bvisc.append(frc_out[1][ii])
        Mass.append(M[ (ii)*6 , (ii)*6])
        print("Massa = %f" % Mass[ii])
        # Cext.append(frc_out[2][ii])
    
        for i in [3,4]:
            if i == 3:
                # GMt.append((C[ii][i,i]+Cext[ii][i,i])/(Mass[ii]*g)) 
                GMt.append((C[ii][i,i])/(Mass[ii]*g)) 
            elif i == 4:
                # GMl.append((C[ii][i,i]+Cext[ii][i,i])/(Mass[ii]*g))
                GMl.append((C[ii][i,i])/(Mass[ii]*g))     
    print('')     
    print('** WAMIT OUTPUT PARAMETERS - HYDROSCRIPTS **')
    print('') 
    print('g = {:.2f} m/s^2'.format(g))
    print('ULEN = {:.2f}'.format(ulen))
    if water_depth_aux == 1:
        print('Water Depth = ' + water_depth)
    else:
        print('Water Depth = {:.2f} m'.format(water_depth))
    for ii in range(NBODY):
        print('Body N = {:d}'.format(ii+1))
        print('  Vols = ' + '[' + ', '.join(["{:.2f}".format(v) for v in vol[ii]]) + '] m^3')
        print('  Mass = {:.2f} t'.format(Mass[ii]))
        print('  CoB = ' + '[' + ', '.join(["{:.2f}".format(v) for v in cb[ii]]) + '] m')
        print('  CoG = ' + '[' + ', '.join(["{:.2f}".format(v) for v in cg[ii]]) + '] m')
        print('  Wamit Axis = ' + '[' + ', '.join(["{:.2f}".format(v) for v in axis[ii]]) + '] m')      
        print('  GMt = {:.2f} m'.format(GMt[ii]))
        print('  GMl = {:.2f} m'.format(GMl[ii]))
    
    return [params, axis, vol, cb, cg, rest_coef, nome_out, GMt, GMl, M, Bvisc, C, Cext]

def read_mmx():
   
    name_mmx = 'force.mmx'
    arq_mmx_aux = open(name_mmx, 'r')
    arq_mmx = arq_mmx_aux.readlines()
    arq_mmx_aux.close()

    N_rows = len(arq_mmx)
    
    cont = 1
    pos = []
    for x in arq_mmx:
        # print(cont)
        # print(x)
        if 'External force matrices:' in x:
            pos.append(cont)
        cont = cont+1
    # print(pos)

    cont = 1
    MMi = []
    BBi = []
    KKi = []
    tamanho_m = []
    for ii in range(len(pos)):
        arq_mmx_aux = open(name_mmx, 'r')
        # I, J, MM_aux, BB_aux, KK_aux = np.loadtxt(arq_mmx_aux, skiprows=pos+1, unpack=True, max)
        if len(pos) == cont:
            I, J, MM_aux, BB_aux, KK_aux = np.genfromtxt(arq_mmx_aux, skip_header=pos[ii]+1, unpack=True)
        else:
            I, J, MM_aux, BB_aux, KK_aux = np.genfromtxt(arq_mmx_aux, skip_header=pos[ii]+1, unpack=True, skip_footer=N_rows-(pos[ii+1]-6))

        arq_mmx_aux.close()

        N_dof = int(np.max(I))

        MMi.append(np.zeros((N_dof, N_dof)))
        BBi.append(np.zeros((N_dof, N_dof)))
        KKi.append(np.zeros((N_dof, N_dof)))

        for i, j, mm, bb, kk in zip(I, J, MM_aux, BB_aux, KK_aux):
            MMi[cont-1][int(i)-1, int(j)-1] = mm
            BBi[cont-1][int(i)-1, int(j)-1] = bb
            KKi[cont-1][int(i)-1, int(j)-1] = kk

        tamanho_m.append(len(KKi[cont-1]))
        # print("M - Body %d" % cont)
        # print(tabulate(MMi[cont-1], floatfmt=".2e", tablefmt="fancy_grid"))
        # print("B - Body %d" % cont)
        # print(tabulate(BBi[cont-1], floatfmt=".2e", tablefmt="fancy_grid"))
        # print("K - Body %d" % cont)
        # print(tabulate(KKi[cont-1], floatfmt=".2e", tablefmt="fancy_grid"))

        cont = cont+1

    N_total = np.sum(tamanho_m)
    MM = np.zeros((N_total, N_total))
    BB = np.zeros((N_total, N_total))
    KK = np.zeros((N_total, N_total))
    idx = 0
    for mm, bb, kk, tt in zip(MMi, BBi, KKi, tamanho_m):
        MM[idx:idx+tt, idx:idx+tt] = mm
        BB[idx:idx+tt, idx:idx+tt] = bb
        KK[idx:idx+tt, idx:idx+tt] = kk
        idx = idx + tt


    return [MM, BB, KK]


def read_frc():
    # Read main file FRC to search for FRC files names
    name_main_frc = 'force.frc'
    arq_frc_aux = open(name_main_frc,'r')
    arq_main_frc = arq_frc_aux.readlines()
    arq_frc_aux.close()
    
    nome_main_frc=[]
    for x in arq_main_frc[1:]:
        if '.frc' in x:
            pos=x.find('.frc')+4
            nome_main_frc.append(x[0:pos])
    
    mass=[]
    damp=[]
    rest_coef_ext=[]
            
    for ii in range(len(nome_main_frc)):
        padrao = "[+-]?\\d+(?:\\.\\d+)?(?:[eE][+-]?\\d+)?"
        # Reading Force file
        nome_frc = nome_main_frc[ii]
        arq_frc_aux = open(nome_frc,'r')
        arq_frc = arq_frc_aux.readlines()
        arq_frc_aux.close()
        
        mass.append([])
        for x in arq_frc[5:11]:
            mass[ii].append([float(i) for i in re.findall(padrao,x)])
            
        damp.append([])
        for x in arq_frc[12:12+6]:
            damp[ii].append([float(i) for i in re.findall(padrao,x)])
    
        rest_coef_ext.append([])
        for x in arq_frc[12+7:12+7+6]:
            rest_coef_ext[ii].append([float(i) for i in re.findall(padrao,x)])
            
        mass[ii] = np.array(mass[ii])
        damp[ii] = np.array(damp[ii])
        rest_coef_ext[ii] = np.array(rest_coef_ext[ii])
    
    return [mass,damp,rest_coef_ext]

def raos(plota=0, dof_plot=[1,2,3,4,5,6], inc_plot=[0,45,90,135,180], multi_fig=False, T_lim = [0, 25], param_out=[]):
    # from matplotlib.ticker import FormatStrFormatter
    # from scipy import interpolate
    if not param_out:
        param_out = output_params()

    # Inputs (after must be imported from a configuration file)
    arq4 = np.loadtxt('force.4')
    ULEN = param_out[0][1]
    NBODY = param_out[0][5]
    
    # Column:  0-Period, 1-Incidence angle, 2-DOF, 3-Amp, 4-Phase, 5-Real, 6-Imag.
    # OPTN.4:    PER    BETA    I    Mod(ξi)    Pha(ξi)    Re(ξi)    Im(ξi)

    # Unique with no sort
    _, idx = np.unique(arq4[:, 0], return_index=True)
    per = np.array([arq4[index, 0] for index in sorted(idx)])
    
#    # identify period to remove from analysis
    # pos_per=[]
    # if np.array(remove_per).size > 0:
    #     print('Removendo períodos:')
    #     print(remove_per)
    #     delta_p = .1
    #     for p in remove_per:
    #         pos_per.append(np.logical_not(np.logical_and(per > (p-delta_p), per < (p+delta_p))))
            
    #     pos_per = np.prod(np.array(pos_per),0)==1
    #     per = per[pos_per]

    inc = np.unique(arq4[:, 1])
    dof = np.unique(arq4[:, 2])

    dof_aux = np.arange(1, dof.max()+1)
    dof_aux = dof_aux.reshape((-1, 6))
    dof_aux = dof_aux[:, [3, 4, 5]]
    dof_aux = dof_aux.reshape((-1, 1))

    dim = np.ones(arq4.shape)

    # dimensionalization of arq4 file
    # arq4: non-dimensional file
    # arq4d: dimensional file
    pos = []

    # positions for dimensioning
    for ii in dof_aux:
        pos.append(arq4[:, 2] == ii)
    pos = np.sum(pos, axis=0) == 1
    dim[np.ix_(pos, [3, 5, 6])] = 1/(ULEN**1)
    arq4d = arq4*dim

    # RAO in complex form
    rao_c_aux = arq4d[:, 5] + arq4d[:, 6] * 1j

    # Function to interpolate the RAO in the given incidences
    
    rao = []
    rao_phase = []
    rao_c =[]
    for ii in dof:
        aux = []
        aux2 = []
        aux3 = []
        for jj in per:
            aux.append(arq4d[(arq4d[:, 2] == ii) & (arq4d[:, 0] == jj), 3])
            aux2.append(arq4d[(arq4d[:, 2] == ii) & (arq4d[:, 0] == jj), 4])
            aux3.append(rao_c_aux[(arq4d[:, 2] == ii) & (arq4d[:, 0] == jj)])
        rao.append(aux)
        rao_phase.append(aux2)
        rao_c.append(np.array(aux3))
        
    # plots
    if plota==1:      
        plot_curves('rao', arq_n_d=arq4d, per=per, dof_plot=dof_plot, inc_plot=inc_plot, NBODY=NBODY,multi_fig=multi_fig, T_lim=T_lim)

    print('')
    print(' * Response Amplitude Operator')

    return [rao, rao_phase, per, inc, dof, arq4d, rao_c]

def wave_forces(plota=0,dof_plot=[1,2,3,4,5,6],inc_plot=[0,45,90,135,180],multi_fig=False, T_lim = [0, 25], param_out=[]):

    if not param_out:
        param_out = output_params()

    arq2 = np.loadtxt('force.2')
    ULEN = param_out[0][1]
    NBODY = param_out[0][5]
    # Unique with no sort
    _, idx = np.unique(arq2[:, 0], return_index=True)
    per = np.array([arq2[index, 0] for index in sorted(idx)])

    inc = np.unique(arq2[:, 1])
    dof = np.unique(arq2[:, 2])

    dim = np.ones(arq2.shape)

    dof_aux = np.arange(1, dof.max()+1)
    dof_aux = dof_aux.reshape((-1, 6))
    dof_aux = dof_aux[:, [3, 4, 5]]
    dof_aux = dof_aux.reshape((-1, 1))
    
    pos = []

    # positions for dimensioning
    for ii in dof_aux:
        pos.append(arq2[:, 2] == ii)
    pos = np.sum(pos, axis=0) == 1
    dim[np.ix_(pos, [3, 5, 6])] = 1.025*9.80665*(ULEN**3)

    dof_aux = np.arange(1, dof.max()+1)
    dof_aux = dof_aux.reshape((-1, 6))
    dof_aux = dof_aux[:, [0, 1, 2]]
    dof_aux = dof_aux.reshape((-1, 1))
    
    pos = []

    # positions for dimensioning
    for ii in dof_aux:
        pos.append(arq2[:, 2] == ii)
    pos = np.sum(pos, axis=0) == 1
    dim[np.ix_(pos, [3, 5, 6])] = 1.025*9.80665*(ULEN**2)

    arq2d = arq2*dim
    
    wforce = []
    wforce_phase = []
    for ii in dof:
        aux = []
        aux2 = []
        for jj in per:
            aux.append(arq2d[(arq2d[:, 2] == ii) & (arq2d[:, 0] == jj), 3])
            aux2.append(arq2d[(arq2d[:, 2] == ii) & (arq2d[:, 0] == jj), 4])
        wforce.append(aux)
        wforce_phase.append(aux2)
        
    # plots
    if plota==1:      
        plot_curves(tipo='wf', arq_n_d=arq2d, per=per, dof_plot=dof_plot, inc_plot=inc_plot, NBODY=NBODY, multi_fig=multi_fig, T_lim=T_lim)

    print('')
    print(' * Wave forces')  
    return [wforce, wforce_phase, arq2d]

def drift_forces(plota=0, drift_analysis_type = 'm', dof_plot=[1,2,6], inc_plot=[0,45,90,135,180], multi_fig=False, T_lim = [0, 25], param_out=[]):
    
    dt_arq_name = {'m': 'force.8', 'p': 'force.9', 'c': 'force.7'}
    dt_message = {'m': 'Momentum', 'p': 'Pressure', 'c': 'Control Surface'}
    arq_name = dt_arq_name[drift_analysis_type]
    print('')  
    print(' * Drift Analysis: ' + dt_message[drift_analysis_type])
    
    if not param_out:
        param_out = output_params()

    arq8 = np.loadtxt(arq_name)
    ULEN = param_out[0][1]
    NBODY = param_out[0][5]
    # Unique with no sort
    _, idx = np.unique(arq8[:, 0], return_index=True)
    per = np.array([arq8[index, 0] for index in sorted(idx)])

    inc = np.unique(arq8[:, 1])
    dof = np.unique(arq8[:, 3])
    dof = dof[dof>0]

    dim = np.ones(arq8.shape)

    dof_aux = np.arange(1, dof.max()+1)
    dof_aux = dof_aux.reshape((-1, 6))
    dof_aux = dof_aux[:, [3, 4, 5]]
    dof_aux = dof_aux.reshape((-1, 1))
    
    pos = []

    # positions for dimensioning
    for ii in dof_aux:
        pos.append(arq8[:, 3] == ii)
    pos = np.sum(pos, axis=0) == 1
    
    dim[np.ix_(pos, [4, 6, 7])] = 1.025*9.80665*(ULEN**2)

    dof_aux = np.arange(1, dof.max()+1)
    dof_aux = dof_aux.reshape((-1, 6))
    dof_aux = dof_aux[:, [0, 1, 2]]
    dof_aux = dof_aux.reshape((-1, 1))
    
    pos = []

    # positions for dimensioning
    for ii in dof_aux:
        pos.append(arq8[:, 3] == ii)
    pos = np.sum(pos, axis=0) == 1
    dim[np.ix_(pos, [4, 6, 7])] = 1.025*9.80665*(ULEN**1)

    arq8d = arq8*dim
    
    wdforce = []
    wdforce_phase = []
    new_dof = np.arange(1, dof.max()+1)
    
    for ii in new_dof:
        aux = []
        aux2 = []
        for jj in per:
            if np.isin(ii,dof):
                aux.append(arq8d[(arq8d[:, 3] == ii) & (arq8d[:, 0] == jj), 4])
                aux2.append(arq8d[(arq8d[:, 3] == ii) & (arq8d[:, 0] == jj), 5])
            else:
                aux.append(np.zeros(len(inc)))
                aux2.append(np.zeros(len(inc)))
                
        wdforce.append(aux)
        wdforce_phase.append(aux2)
        
    # plots
    if plota==1:      
        plot_curves(tipo='mdf', arq_n_d=arq8d, per=per, dof_plot=dof_plot, inc_plot=inc_plot, NBODY=NBODY, multi_fig=multi_fig, dt=drift_analysis_type, T_lim=T_lim)
    
    return [wdforce, wdforce_phase, arq8d]

def added_mass_pot_damping(plota=0, dof_plot=[1,2,3,4,5,6], multi_fig=False, T_lim = [0, 25], param_out=[]):
    
    if not param_out:
        param_out = output_params()

    ULEN = param_out[0][1]
    NBODY = param_out[0][5]
    arq1 = np.loadtxt('force.1')
    
    # print('added_mass_pot_damping: ULEN = {:.1f}'.format(ULEN))
    # Unique with no sort
    _, idx = np.unique(arq1[:, 0], return_index=True)
    per = np.array([arq1[index, 0] for index in sorted(idx)])
    # searching the evaluated dof's
    dof = arq1[arq1[:,0]==per[0],1:3]
    
    dof_aux = [1,2,3]
    aux1 = []
    for n in range(NBODY-1):
        #print(n)
        for z in dof_aux:
            aux1.append(z+6*(n+1))
    
    [dof_aux.append(n) for n in aux1]
    
    #dim = np.ones(arq1.shape)
    dim = []
    for ii in arq1:
        if ii[1]==ii[2]:
            if ii[1] in dof_aux:
                k=3
            else:
                k=4
        else:
            k=5
        dim.append([1,1,1,1.025*(ULEN**k),(2*np.pi/ii[0])*1.025*(ULEN**k)])
        
    arq1d = arq1*dim
    added_mass=[]
    pot_damp=[]
    
    dof1 = np.unique(arq1[:,1:3],axis=0)
    
    #dof1 = [[1,1],[1,3],[1,5],[2,2],[2,4],[2,6],[3,1],[3,3],[3,5],[4,2],[4,4],[4,6],[5,1],[5,3],[5,5],[6,2],[6,4],[6,6]]
    
    # Added Mass
    # 
    # Plane motion ->  lower frequency
    # - At least 1 dof is surge, sway or yaw. Ex: [1,1],[1,3]
    # Out-plane motion -> mean between the mean and the max added mass
    # - Both dof are out-plane. Ex: [3,3], [3,4]
    
    aux_mad_matrix = []
    aux_pdamp_matrix = []
    for x in dof1:      
        aux = [] #massa adicional
        aux2 = [] #amorteciment
        for jj in per:
            aux.append(arq1d[(arq1d[:, 1] == x[0]) & (arq1d[:, 2] == x[1]) & (arq1d[:, 0] == jj),3])
            aux2.append(arq1d[(arq1d[:, 1] == x[0]) & (arq1d[:, 2] == x[1]) & (arq1d[:, 0] == jj),4])
        added_mass.append(aux)
        pot_damp.append(aux2)
        # print(x)
        if (x == [3,3]).all() or (x == [3,5]).all() or (x == [4,4]).all() or (x == [5,3]).all() or (x == [5,5]).all():
            aux_mad_matrix.append(np.mean([np.mean(aux),np.max(aux)]))
            aux_pdamp_matrix.append(np.mean([np.mean(aux2),np.max(aux2)]))
        else:
            pos_min_freq = np.argmax(per)
            aux_mad_matrix.append(np.array(aux[pos_min_freq]))
            aux_pdamp_matrix.append(np.array(aux2[pos_min_freq]))
    
    aux_mad_matrix = np.array(aux_mad_matrix)
    aux_pdamp_matrix = np.array(aux_pdamp_matrix)
    added_mass_matrix = np.zeros((6*NBODY, 6*NBODY))
    pot_damp_matrix = np.zeros((6*NBODY, 6*NBODY))
    
    cont = 0
    for x in dof1:
        pos1 = int(x[0]-1)
        pos2 = int(x[1]-1)
        added_mass_matrix[pos1,pos2] = aux_mad_matrix[cont]
        pot_damp_matrix[pos1,pos2] = aux_pdamp_matrix[cont]
        cont+=1
    
    added_mass = np.transpose(added_mass)
    added_mass = added_mass[0]
    pot_damp = np.transpose(pot_damp)
    pot_damp = pot_damp[0]

    # plots
    if plota==1:      
        plot_curves('a', arq1d, per, dof_plot, [], NBODY=NBODY, multi_fig=multi_fig, T_lim=T_lim)
        plot_curves('b', arq1d, per, dof_plot, [], NBODY=NBODY,multi_fig=multi_fig, T_lim=T_lim)
    
    print('')
    print(' * Added Mass and Potential Damping')
    return [added_mass, pot_damp, dof1, arq1d, added_mass_matrix, pot_damp_matrix]

def dynamic_params(param_out, mad):
    # function to evaluate the dynamic parameter as Natural Periods, Viscous Damping coefs
    print(' * Evaluating Dynamic Parameters')

    # [params, axis, vol, cb, cg, rest_coef, nome_out, GMt, GMl, M, Bvisc, C, Cext] = param_out
    # [added_mass, pot_damp, dof1, arq1d, added_mass_matrix, pot_damp_matrix] = mad

    NBODY = param_out[0][5]
    M = param_out[9]
    Bvisc = param_out[10]
    C = param_out[11]
    Cext = param_out[12]
    MA = mad[4]
    Bpot = mad[5]

    Tn = []
    Bc = []
    ca = []
    cv = []
    for jj in range(NBODY):
        tn_aux = []
        bc_aux = []
        cv_aux = []
        ca_aux = []
        print('BODY ' + str(jj+1))
        for ii in range(6):
            if (C[jj][ii][ii] + Cext[jj][ii][ii]) != 0:
                tn_aux.append(2 * np.pi * np.sqrt( (M[jj][ii][ii] + MA[ii][ii]) / (C[jj][ii][ii] + Cext[jj][ii][ii]) ))
                bc_aux.append(2 * np.sqrt( (M[jj][ii][ii] + MA[ii][ii]) * (C[jj][ii][ii] + Cext[jj][ii][ii]) ))
                cv_aux.append(Bvisc[jj][ii][ii] / bc_aux[ii])
                ca_aux.append((Bvisc[jj][ii][ii] + Bpot[ii][ii]) / bc_aux[ii])
            else:
                tn_aux.append(0.0)
                bc_aux.append(0.0)
                cv_aux.append(0.0)
                ca_aux.append(0.0)
            print(' M_' + str(ii+1) + ' = {:.2f}'.format(M[jj][ii][ii]) + \
                ' MA_' + str(ii+1) + ' = {:.2f}'.format(MA[ii][ii]) + \
                ' C_' + str(ii+1) + ' = {:.2f}'.format((C[jj][ii][ii] + Cext[jj][ii][ii])) + \
                ' Tn_' + str(ii+1) + ' = {:.2f}'.format(tn_aux[ii]) + \
                ' Bpot_' + str(ii+1) + ' = {:.2f}'.format(Bpot[ii][ii]) + \
                ' Bvisc_' + str(ii+1) + ' = {:.2f}'.format(Bvisc[jj][ii][ii]) + \
                ' cv_' + str(ii+1) + ' = {:.2f}'.format(cv_aux[ii]) + \
                ' ca_' + str(ii+1) + ' = {:.2f}'.format(ca_aux[ii])  )

                # RAO = F / [-(M+A)*w^2 + (B+Bext)*w + (K+Kext)]
    # 1- Evaluate the critical damping
    # 2- Evaluate the matrices M, A, B, Bext, K and Kext

def point_rao(points):
    # function to evaluate the rao in specific points    
    # Entry:
    # points = [[x1, y1, z1], [x2, y2, z2], ...]
    
    points = np.array(points)
    
    out_rao = raos(0)
    
#    per = out_rao[2]
#    inc = out_rao[3]
#    dof = out_rao[4]
    rao_c = out_rao[6]
    
    rao_p_i = []
    rao_p_j = []
    rao_p_k = []
    for pt in points:
        rao_p_i.append([ rao_c[0] - pt[1] * rao_c[5] + pt[2] * rao_c[4]])
        rao_p_j.append([ rao_c[1] - pt[2] * rao_c[3] + pt[0] * rao_c[5]])              
        rao_p_k.append([ rao_c[2] - pt[0] * rao_c[4] + pt[1] * rao_c[3]])
    
    return [rao_p_i, rao_p_j, rao_p_k]

#debuggers
#op = output_params()
#raos(1,remove_per=[18,19.3])
#raos(1)
#wave_forces(1,[1],[0])
#drift_forces_momentum(1)
#added_mass_pot_damping(plota=0)
#point_rao([[100, 20, 30], [-100, 10, 0]])
#verify_NBODY()
