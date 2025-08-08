import json
import feasst
import numpy as np

from Set_up_Atmo_check import Build_GCMC_json, run_gcmc

# Run a single GCMC simulation of N2 in a MOF at 1atm

# Settings for the current run
Sorbent = 'ZOYKAB_clean'
Calc_args = Build_GCMC_json(300, 1.01325, Sorbent, ['N2_Trappe'], [1.])
Calc_args['seed'] = '92532459'

# Ouptut Watermark
print('FEASST Version: ', feasst.version())
print()
print('Simulation Arguments:')
print(Calc_args)
print()

# Actual Calculation
results_300 = run_gcmc(Calc_args)

# Output
output_file = f'results_300K_1atmN2_{Sorbent}.json'
with open(output_file, "w") as f:
    json.dump(results_300, f, indent=4)
print('Results stored in file: ', output_file)
