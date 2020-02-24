import analysis_wamit

p_out = analysis_wamit.output_params()

a = analysis_wamit.added_mass_pot_damping(param_out=p_out)

dp = analysis_wamit.dynamic_params(param_out=p_out, mad=a)