import json
import feasst

from Set_up_Kh_calc import Build_calc_json
from Kh_utilities import atomistic_system_generator, henry_AlwaysReject

# Compute Kh and extrapolation coefficients

# Settings for the current run
Sorbent = 'ZOYKAB_clean'
Calc_args = Build_calc_json(300., Sorbent, 'N2_Trappe')
Calc_args['seed'] = 1200512782

# Output Watermark
print('FEASST Version: ', feasst.version())
print()
print('Simulation Arguments:')
print(Calc_args)
print()

# Actual Calculation
results_300 = henry_AlwaysReject(**Calc_args)

# Output
output_file = "results_300_N2_" + Sorbent.replace('_clean','') + ".json"
# make file string-compatible
del results_300['system_generator_class']
with open(output_file, "w") as f:
    json.dump(results_300, f, indent=4)
print('Results stored in file: ', output_file)
