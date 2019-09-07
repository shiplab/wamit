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
    a = analysis_wamit.added_mass_pot_damping(1, dof_plot, multi_fig=True)

if rao == True:
    r = analysis_wamit.raos(1, dof_plot, inc_plot, multi_fig=True)

if drift == True:
    dof_aux=[]
    for dof in dof_plot:
        if dof != 3 and dof!=4 and dof!=5:
            dof_aux.append(dof)

    d = analysis_wamit.drift_forces_momentum(1, dof_aux, inc_plot, multi_fig=True)

if save_fig == True:
    for i in plt.get_fignums():
        f = plt.figure(i)
        name = f.canvas.get_window_title()
        name = name.replace(' ','_')
        plt.savefig('fig_' + name + '.pdf')

plt.show()