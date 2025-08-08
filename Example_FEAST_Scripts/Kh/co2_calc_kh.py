import json
import feasst

from Set_up_Kh_calc import Build_calc_json
from Kh_utilities import atomistic_system_generator, henry_AlwaysReject

# Compute Kh and extrapolation coefficients

# Settings for the current run
Sorbent = 'ZOYKAB_clean'
Calc_args = Build_calc_json(350., Sorbent, 'CO2_Trappe')
Calc_args['seed'] = 1200527482

# Output Watermark
print('FEASST Version: ', feasst.version())
print()
print('Simulation Arguments:')
print(Calc_args)
print()

# Actual Calculation
results_350 = henry_AlwaysReject(**Calc_args)

# Output
output_file = Sorbent + ".results.json"
# make file string-compatible
del results_350['system_generator_class']
with open(output_file, "w") as f:
    json.dump(results_350, f, indent=4)
print('Results stored in file: ', output_file)
