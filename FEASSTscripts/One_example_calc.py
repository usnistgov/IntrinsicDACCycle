# -*- coding: utf-8 -*-
"""
Created on Thu May 12 14:38:53 2022

@author: Austin McDannald
"""

import json
from Set_up_Kh_calc import Build_calc_json
#import sys
#sys.path.insert(1, '/home/asm6/software/FEAST/porous_network_utils')
#from porous_network_utils import atomistic_system_generator, henry_AlwaysReject
#sys.path.insert(0, '/home/asm6/software/FEAST/porous_network_utils/porous_network_utils/')
#import atomistic_system_generator, henry_AlwaysReject
from porous_networks_utils.porous_network_utils import atomistic_system_generator, henry_AlwaysReject


Calc_args = Build_calc_json(350., 'ABAVIJ_clean','CO2_Trappe')

results_350 = henry_AlwaysReject(**Calc_args)
del results_350['system_generator_class']

with open('/users/asm6/results_350_ABAVIJ.json', 'w') as f:
#with open('results_350_ABAVIJ.json', 'w') as f:
    json.dump(results_350, f, indent=4)