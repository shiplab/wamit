import numpy as np
import matplotlib.pyplot as plt
import re
import math

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
    
def output_params(): 
    # Reading running output parameters
    nome_out = 'force.out'
    arq_out_aux = open(nome_out,'r')
    arq_out = arq_out_aux.readlines()
    arq_out_aux.close()
    padrao = "[+-]?\\d+(?:\\.\\d+)?(?:[eE][+-]?\\d+)?"
        
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
        if ('XBODY =' in x)==True:
            [xbody,ybody,zbody,phibody] = [float(i) for i in re.findall(padrao,x)]
        if ('Volumes (VOLX,VOLY,VOLZ):' in x)==True:
            [xvol,yvol,zvol] = [float(i) for i in re.findall(padrao,x)]
        if ('Center of Buoyancy (Xb,Yb,Zb):' in x)==True:
            [xb,yb,zb] = [float(i) for i in re.findall(padrao,x)]
        if ('C(3,3),C(3,4),C(3,5):' in x)==True:
            [a0,a1,a2,a3,a4,a5,c33,c34,c35] = [float(i) for i in re.findall(padrao,x)]
        if ('C(4,4),C(4,5),C(4,6):' in x)==True:
            [a0,a1,a2,a3,a4,a5,c44,c45,c46] = [float(i) for i in re.findall(padrao,x)]
        if ('C(5,5),C(5,6):' in x)==True:
            [a0,a1,a2,a3,c55,c56] = [float(i) for i in re.findall(padrao,x)]
        if ('Center of Gravity  (Xg,Yg,Zg):' in x)==True:
            [xg,yg,zg] = [float(i) for i in re.findall(padrao,x)]
                     
    params = [g,ulen,rho,water_depth,water_depth_aux]
    axis = [xbody,ybody,zbody,phibody]
    vol = [xvol,yvol,zvol]
    cb = [xb,yb,zb]
    cg = [xg,yg,zg]
    rest_coef = [c33*g*rho*(ulen**2),c34*g*rho*(ulen**3),c35*g*rho*(ulen**3),c44*g*rho*(ulen**4),c45*g*rho*(ulen**4),c46*g*rho*(ulen**4),c55*g*rho*(ulen**4),c56*g*rho*(ulen**4)]
    
    dof_rest_coef =[[3,3],[3,4],[3,5],[4,4],[4,5],[4,6],[5,5],[5,6]]
    
    C = np.zeros((6,6))
    cont=0
    for x in dof_rest_coef:
        C[x[0]-1,x[1]-1] = rest_coef[cont]
        cont+=1
        
    mad = added_mass_pot_damping()
    mad = mad[4]
    
    frc_out = read_frc() #substitute by read_mmx after
    M = frc_out[0]
    Mass = M[0,0]
    Cext = frc_out[2]
    
    for i in [3,4]:
        if i == 3:
            GMt = (C[i,i]+Cext[i,i])/(Mass*g)
        elif i == 4:
            GMl = (C[i,i]+Cext[i,i])/(Mass*g)        
        
    
    print('g = {:.2f} m/s^2'.format(g))
    print('ULEN = {:.2f}'.format(ulen))
    if water_depth_aux == 1:
        print('Water Depth = ' + water_depth)
    else:
        print('Water Depth = {:.2f} m'.format(water_depth))
    print('Vols = ' + '[' + ', '.join(["{:.2f}".format(v) for v in vol]) + '] m^3')
    print('Mass = {:.2f} t'.format(Mass))
    print('CoB = ' + '[' + ', '.join(["{:.2f}".format(v) for v in cb]) + '] m')
    print('CoG = ' + '[' + ', '.join(["{:.2f}".format(v) for v in cg]) + '] m')
    print('Wamit Axis = ' + '[' + ', '.join(["{:.2f}".format(v) for v in axis]) + '] m')      
    print('GMt = {:.2f} m'.format(GMt))
    print('GMl = {:.2f} m'.format(GMl))
    
    return [params,axis,vol,cb,cg,rest_coef,nome_out]
    
def read_frc():
    padrao = "[+-]?\\d+(?:\\.\\d+)?(?:[eE][+-]?\\d+)?"
    # Reading Force file
    nome_frc = 'ship.frc'
    arq_frc_aux = open(nome_frc,'r')
    arq_frc = arq_frc_aux.readlines()
    arq_frc_aux.close()
    
    mass=[]
    for x in arq_frc[5:11]:
        mass.append([float(i) for i in re.findall(padrao,x)])
        
    damp=[]
    for x in arq_frc[12:12+6]:
        damp.append([float(i) for i in re.findall(padrao,x)])

    rest_coef_ext=[]
    for x in arq_frc[12+7:12+7+6]:
        rest_coef_ext.append([float(i) for i in re.findall(padrao,x)])
        
    mass = np.array(mass)
    damp = np.array(damp)
    rest_coef_ext = np.array(rest_coef_ext)
    
    return [mass,damp,rest_coef_ext]

def plot_curves(tipo,arq_n_d,per,dof_plot,inc_plot):
    cont = 0
    curve = []
    ax = []
    #Visualization parameter
    t_inf = 0
    t_sup = 25
    #Name parameters
    names_dof = ['Surge', 'Sway', 'Heave',
                 'Roll', 'Pitch', 'Yaw']
    
    deg_multiplier = 1
    
    if tipo == 'rao':
        col_inc=1
        col_dof=2
        col_plot=3
        name = 'Response Amplitude Operator'
        units_dof = ['[m/m]', '[m/m]', '[m/m]',
                     '[deg/m]', '[deg/m]', '[deg/m]']
        deg_multiplier = 180 / np.pi
    elif tipo == 'mdf':
        col_inc=1
        col_dof=3
        col_plot=4
        name = 'Mean drift forces'
        units_dof = ['[kN]', '[kN]', '[kN]',
                     '[kN.m]', '[kN.m]', '[kN.m]']
    elif tipo == 'wf':
        col_inc=1
        col_dof=2
        col_plot=3
        name = 'Wave forces'
        units_dof = ['[kN]', '[kN]', '[kN]',
                     '[kN.m]', '[kN.m]', '[kN.m]']
    
    dof_plot_len = len(dof_plot)
    
    if dof_plot_len<=3:
        n_col = 1
        n_lin = dof_plot_len
    else:
        n_col = 2
        n_lin = math.ceil(dof_plot_len/2)
        
    f1 = plt.figure(None, (5.7*n_col, 2.5*n_lin))
    f1.canvas.set_window_title(name)
    for ii in dof_plot:
        aux = []
        leg = []

        for jj in inc_plot:
            aux.append(arq_n_d[(arq_n_d[:, col_dof] == ii) & (arq_n_d[:, col_inc] == jj), col_plot])
            leg.append('Inc = ' + str(jj) + 'deg')

        aux = np.array(aux)
        curve.append(aux)
        ax = plt.subplot(n_lin, n_col, cont+1)
        plt.grid(axis='both')
        idx = (per >= t_inf) & (per <= t_sup)

        # Defines a de-multiplier for plots that is necessary change units from rad to deg
        if (ii == 4) | (ii == 5) | (ii == 6):
            plt.plot(per, curve[cont].transpose() * deg_multiplier)
            y_inf = curve[cont].min() * deg_multiplier
            y_sup = curve[cont].max() * deg_multiplier
        else:
            plt.plot(per, curve[cont].transpose())
            y_inf = curve[cont][:, idx].min()
            y_sup = curve[cont][:, idx].max()
        plt.xlim([t_inf, t_sup])

        if y_inf != y_sup:
            wd = .05*abs(y_sup - y_inf)
            plt.ylim([y_inf-wd, y_sup+wd])

        if cont == 0:
            plt.legend(leg, loc='best')

        plt.ylabel(names_dof[ii-1] + ' ' + units_dof[ii-1])
        plt.xlabel('T [s]')
        cont+=1
    
    plt.tight_layout()
    plt.show()


def raos(plota=0, dof_plot=[1,2,3,4,5,6], inc_plot=[0,45,90,135,180]):
    # from matplotlib.ticker import FormatStrFormatter
    # from scipy import interpolate
    param_out = output_params()
    # Inputs (after must be imported from a configuration file)
    arq4 = np.loadtxt('force.4')
    ULEN = param_out[0][1]
    NBODY = 1   # implement function to read the number of bodies
    
    # Column:  0-Period, 1-Incidence angle, 2-DOF, 3-Amp, 4-Phase, 5-Real, 6-Imag.
    # OPTN.4:    PER    BETA    I    Mod(ξi)    Pha(ξi)    Re(ξi)    Im(ξi)

    # Unique with no sort
    per_a, idx = np.unique(arq4[:, 0], return_index=True)
    per = np.array([arq4[index, 0] for index in sorted(idx)])

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
        plot_curves('rao',arq4d,per,dof_plot,inc_plot)
    
    return [rao,rao_phase,per,inc,dof,arq4d,rao_c]


def wave_forces(plota=0,dof_plot=[1,2,3,4,5,6],inc_plot=[0,45,90,135,180]):
    param_out = output_params()
    arq2 = np.loadtxt('force.3')
    ULEN = param_out[0][1]
    # Unique with no sort
    per_a, idx = np.unique(arq2[:, 0], return_index=True)
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
        plot_curves('wf',arq2d,per,dof_plot,inc_plot)
        
    return [wforce,wforce_phase,arq2d]

def drift_forces_momentum(plota=0,dof_plot=[1,2,6],inc_plot=[0,45,90,135,180]):
    
    param_out = output_params()
    arq8 = np.loadtxt('force.8')
    ULEN = param_out[0][1]
    t_inf = 5
    t_sup = 25
    inc_plot = np.arange(0, 181, 45)
    # Unique with no sort
    per_a, idx = np.unique(arq8[:, 0], return_index=True)
    per = np.array([arq8[index, 0] for index in sorted(idx)])

    inc = np.unique(arq8[:, 1])
    dof = np.unique(arq8[:, 3])

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
    dof = [1,2,3,4,5,6]
    for ii in dof:#[dof[i] for i in [0,1,5]]:
        aux = []
        aux2 = []
        for jj in per:
            if ii==3 or ii==4 or ii==5:
                aux.append(np.zeros(len(inc)))
                aux2.append(np.zeros(len(inc)))
            else:
                aux.append(arq8d[(arq8d[:, 3] == ii) & (arq8d[:, 0] == jj), 4])
                aux2.append(arq8d[(arq8d[:, 3] == ii) & (arq8d[:, 0] == jj), 5])
                
        wdforce.append(aux)
        wdforce_phase.append(aux2)
        
    # plots
    if plota==1:      
        plot_curves('mdf',arq8d,per,dof_plot,inc_plot)
    
    return [wdforce,wdforce_phase,arq8d]

def added_mass_pot_damping(plota=0):
    arq1 = np.loadtxt('force.1')
    ULEN = read_ULEN()
#    print('added_mass_pot_damping: ULEN = {:.1f}'.format(ULEN))
    # Unique with no sort
    per_a, idx = np.unique(arq1[:, 0], return_index=True)
    per = np.array([arq1[index, 0] for index in sorted(idx)])
    
    dim = np.ones(arq1.shape)
    
    dim = []
    
    for ii in arq1:
        if ii[1]==ii[2]:
            if ii[1]==1 or ii[1]==2 or ii[1]==3:
                k=3
            else:
                k=4
        else:
            k=5
   
        dim.append([1,1,1,1.025*(ULEN**k),(2*np.pi/ii[0])*1.025*(ULEN**k)])
        
    arq1d = arq1*dim
    added_mass=[]
    pot_damp=[]
    dof1 = [[1,1],[1,3],[1,5],[2,2],[2,4],[2,6],[3,1],[3,3],[3,5],[4,2],[4,4],[4,6],[5,1],[5,3],[5,5],[6,2],[6,4],[6,6]]
    
    # Added Mass
    # 
    # Plane motion ->  lower frequency
    # Out-plane motion -> mean between the mean and the max added mass
    
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
        
        if x == [3,3] or x == [3,5] or x == [4,4] or x == [5,3] or x == [5,5]:
            aux_mad_matrix.append(np.mean([np.mean(aux),np.max(aux)]))
            aux_pdamp_matrix.append(np.mean([np.mean(aux2),np.max(aux2)]))
        else:
            pos_min_freq = np.argmax(per)
            aux_mad_matrix.append(np.array(aux[pos_min_freq]))
            aux_pdamp_matrix.append(np.array(aux2[pos_min_freq]))
    
    aux_mad_matrix = np.array(aux_mad_matrix)
    aux_pdamp_matrix = np.array(aux_pdamp_matrix)
    added_mass_matrix = np.zeros((6,6))
    pot_damp_matrix = np.zeros((6,6))
    
    cont = 0
    for x in dof1:
        pos1 = x[0]-1
        pos2 = x[1]-1
        added_mass_matrix[pos1,pos2] = aux_mad_matrix[cont]
        pot_damp_matrix[pos1,pos2] = aux_pdamp_matrix[cont]
        cont+=1
    
    added_mass = np.transpose(added_mass)
    added_mass = added_mass[0]
    pot_damp = np.transpose(pot_damp)
    pot_damp = pot_damp[0]
    
    return [added_mass,pot_damp,dof1,arq1d,added_mass_matrix,pot_damp_matrix]

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
#raos(1)
#wave_forces(1)
#drift_forces_momentum(1)
#added_mass_pot_damping(plota=0)
#point_rao([[100, 20, 30], [-100, 10, 0]])