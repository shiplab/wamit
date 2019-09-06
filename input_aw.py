# this function generates the file input_aw.txt in the wamit case folder. So user can adjust options and execute the complete analysis

input = open('input_aw.txt','w')

s='''# Angles to plot
inc_plot = [0,45,90]
# Degrees of freedom to plot
dof_plot = [1,2,3,4,5,6]

# Desired Analysis
added_mass = True
rao = True
drift = True

# Save Figures
save_fig = False
'''
input.write(s)

input.close()

print('Analysis Wamit input file "input_aw.txt" generated.')