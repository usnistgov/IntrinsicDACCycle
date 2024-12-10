import multiprocessing
import tqdm
import glob
import json

from monte_carlo_saturation_adsorption import Build_saturation_calc_json, estimate_max_loading


list_of_sorbents = glob.glob('/users/asm6/DAC_data/CSD_FEASST_Materials/Materials/*.json')

#stip off the file extension and directory from the strings leaving just the sorbent name.
list_of_sorbents = [sub.replace(".json", "") for sub in list_of_sorbents]
list_of_sorbents = [sub.replace("/users/asm6/DAC_data/CSD_FEASST_Materials/Materials/", "") for sub in list_of_sorbents]



def Run_a_saturation_calc(sorbent):

    #build the calculation arguments
    Calc_args = Build_saturation_calc_json(200., sorbent, 'CO2_Trappe')

    short_isotherm, Calc_args, mc_runner, isotherm_dict, ewald_args = estimate_max_loading(threshold=0.005,**Calc_args)
    
    Calc_args["ewald_args"] = ewald_args
    Calc_args["isotherm"] = isotherm_dict


    with open(f"/users/asm6/DAC_data/Saturation/CO2_sat_200K_{sorbent}.json", "w") as f:
        json.dump(Calc_args, f, indent=4)


if __name__ == "__main__":
    p = multiprocessing.Pool(64)
    for _ in tqdm.tqdm(p.imap_unordered(Run_a_saturation_calc, list_of_sorbents), total = len(list_of_sorbents)):
        pass