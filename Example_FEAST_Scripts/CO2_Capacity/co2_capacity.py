import json
import feasst

from FEASST_saturation_adsorption import Build_saturation_calc_json, estimate_max_loading

# Compute an isotherm of CO2 in a MOF to check for saturation

# Settings for the current run
Sorbent = 'ZOYKAB_clean'
Calc_args = Build_saturation_calc_json(200., Sorbent, 'CO2_Trappe')

# Output Watermark
print('FEASST Version: ', feasst.version())
print()
print('Simulation Arguments:')
print(Calc_args)
print()

# Run the isotherm loop
short_isotherm, Calc_args, mc_runner, isotherm_dict, ewald_args = estimate_max_loading(threshold=0.005,**Calc_args)
Calc_args["Saturation_isotherm"] = isotherm_dict

# Outpu
output_file = "CO2_sat_200K_"+Sorbent+".json"
with open(output_file, "w") as f:
    json.dump(Calc_args, f, indent=4)
print('Results stored in file: ', output_file)
