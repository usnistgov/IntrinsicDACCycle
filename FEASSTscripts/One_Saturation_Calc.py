import glob
import json
import sys

from monte_carlo_saturation_adsorption import Build_saturation_calc_json, estimate_max_loading


def Run_a_saturation_calc(sorbent):

    #chop of the file extenision and directory, leaving only the sorbent name
    sorbent = sorbent.replace(".json", "")
    sorbent = sorbent.replace("/users/asm6/DAC_data/CSD_FEASST_Materials/Materials/", "")

    #Check if that sorbent has already been done
    list_of_completed = glob.glob('/tmp/asm6/Saturation/*.json')
    #stip off the file extension and directory from the strings leaving just the sorbent name.
    list_of_completed = [sub.replace(".json", "") for sub in list_of_completed]
    list_of_completed= [sub.replace("/tmp/asm6/Saturation/CO2_sat_200K_", "") for sub in list_of_completed]

    test = sorbent in list_of_completed

    if test == True:
        pass
    else:

        #build the calculation arguments
        Calc_args = Build_saturation_calc_json(200., sorbent, 'CO2_Trappe')

        short_isotherm, Calc_args, mc_runner, isotherm_dict, ewald_args = estimate_max_loading(threshold=0.005,**Calc_args)
        
        Calc_args["ewald_args"] = ewald_args
        Calc_args["isotherm"] = isotherm_dict


        with open(f"/tmp/asm6/Saturation/CO2_sat_200K_{sorbent}.json", "w") as f:
            json.dump(Calc_args, f, indent=4)

print("Running", sys.argv[1])

sorbent = sys.argv[1]
Run_a_saturation_calc(sorbent)