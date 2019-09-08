import matplotlib.pyplot as plt
import numpy as np
import math

def plot_curves(tipo, arq_n_d, per, dof_plot, inc_plot, multi_fig=False, dt = 'm', T_lim = [0, 30]):
    cont = 0
    curve = []
    # ax = []
    #Visualization parameter
    t_inf = T_lim[0]
    t_sup = T_lim[1]
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
        name = 'Mean drift forces: '
        pos_name = {'m': 'Momentum', 'p': 'Pressure', 'c': 'Control Surface'}
        name = name + pos_name[dt]
        units_dof = ['[kN]', '[kN]', '[kN]',
                     '[kN.m]', '[kN.m]', '[kN.m]']
    elif tipo == 'wf':
        col_inc=1
        col_dof=2
        col_plot=3
        name = 'Wave forces'
        units_dof = ['[kN]', '[kN]', '[kN]',
                     '[kN.m]', '[kN.m]', '[kN.m]']
    elif tipo == 'a':
        col_inc=2
        col_dof=1
        col_plot=3
        name = 'Added Mass'
        units_dof = ['[t]', '[t]', '[t]',
                     '[t.m]', '[t.m]', '[t.m]']
    elif tipo == 'b':
        col_inc=2
        col_dof=1
        col_plot=4
        name = 'Potential Damping'
        units_dof = ['[t/s]', '[t/s]', '[t/s]',
                     '[t.m/s]', '[t.m/s]', '[t.m/s]']
    
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

        if len(inc_plot) == 0:
            it_aux = [ii]
        else:
            it_aux = inc_plot

        for jj in it_aux:
            aux.append(arq_n_d[(arq_n_d[:, col_dof] == ii) & (arq_n_d[:, col_inc] == jj), col_plot])
            leg.append('Inc = ' + str(jj) + 'deg')

        aux = np.array(aux)
        curve.append(aux)
        ax = plt.subplot(n_lin, n_col, cont+1)
        plt.grid(axis='both')
        idx = (per >= t_inf) & (per <= t_sup)

        # Defines a deg-multiplier for plots that is necessary change units from rad to deg
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
            if len(inc_plot) > 0:
                plt.legend(leg, loc='best')

        plt.ylabel(names_dof[ii-1] + ' ' + units_dof[ii-1])
        plt.xlabel('T [s]')
        cont+=1
    
    plt.tight_layout()
    if multi_fig == False:
        plt.show()