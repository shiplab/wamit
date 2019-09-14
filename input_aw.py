# this function generates the file input_aw.txt in the wamit case folder. So user can adjust options and execute the complete analysis

input = open('input_aw.txt','w')

s='''# Angles to plot
inc_plot = list(range(0, 181, 45))

# Degrees of freedom to plot
dof_plot = [1,2,3,4,5,6]

# Limit of period axis
T_lim = [0,30]

# Desired Analysis
added_mass = True
wave_forces = True
rao = True
drift = True
drift_analysis = ['m'] # 'm' - momentum, 'p' - pressure, 'c' - control surface

# WNF Mean Drift
wnf_drift = ['m']

# Figures
save_fig = True
show_fig = False
'''
input.write(s)

input.close()

print(' ')
print('Analysis Wamit input file "input_aw.txt" generated.')
print(' ')