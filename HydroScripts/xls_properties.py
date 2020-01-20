import analysis_wamit
import gdf_class

p_out = analysis_wamit.output_params()

print(p_out)

ship = gdf_class.GDF('ship.gdf')

ship.props(T=[], KG=[])

print(ship)