import os.path
import matplotlib.pyplot as plt

if os.path.isfile('input_aw.txt'):
    input = open('input_aw.txt','r')
else:
    import input_aw
    input = open('input_aw.txt','r')

contents = input.read()

contents = contents.split('\n')

for line in contents:
    exec(line)
    
import analysis_wamit

if added_mass == True:
    a = analysis_wamit.added_mass_pot_damping(plota=1, dof_plot=dof_plot, multi_fig=True, T_lim=T_lim)

if rao == True:
    r = analysis_wamit.raos(plota=1, dof_plot=dof_plot, inc_plot=inc_plot, multi_fig=True, T_lim=T_lim)

if drift == True:
    dof_aux=[]
    for dof in dof_plot:
        if dof != 3 and dof!=4 and dof!=5:
            dof_aux.append(dof)
    for dt in drift_analysis:
        d = analysis_wamit.drift_forces(plota=1, drift_analysis_type=dt, dof_plot=dof_aux, inc_plot=inc_plot, multi_fig=True, T_lim=T_lim)

if save_fig == True:
    for i in plt.get_fignums():
        f = plt.figure(i)
        name = f.canvas.get_window_title()
        name = name.replace(' ','_')
        name = name.replace(':','')
        plt.savefig('fig_' + name + '.svg')

plt.show()