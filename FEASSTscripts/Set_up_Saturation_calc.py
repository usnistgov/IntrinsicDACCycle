import json
import numpy as np


from monte_carlo_saturation_adsorption import Build_saturation_calc_json, estimate_max_loading


Calc_args = Build_saturation_calc_json(300., 'ABAVIJ_clean', 'CO2_Trappe')


# Converge to 0.1%
short_isotherm, mc_runner = estimate_max_loading(threshold=0.001,**Calc_args)
for line in short_isotherm:
    print(line)
print('Estimated Maximum Loading', short_isotherm[-1][1], short_isotherm[-1][2])


Calc_args["Saturation_isotherm"] = short_isotherm

with open("/users/asm6/CO2_sat_300K_ABAVIJ.json", "w") as f:
    json.dump(Calc_args, f, indent=4)