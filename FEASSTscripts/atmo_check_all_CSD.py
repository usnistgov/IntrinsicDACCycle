# -*- coding: utf-8 -*-
"""
Created on Fri May 13 16:31:26 2022

@author: Austin McDannald
"""
import os
cpu_s = os.cpu_count()

import multiprocessing
import tqdm
import glob
import json

from Set_up_Atmo_check import Build_GCMC_json, run_gcmc


#Get all the files names of the materials
sorbent_mat_files = glob.glob('/wrk/asm6/CSD_data/CSD_FEASST_Materials/Materials/*.json')
sorbent_mat_names = []
sorbent_mat_names_parsed = []
for mat in sorbent_mat_files:
    mat1 = mat.replace("/wrk/asm6/CSD_data/CSD_FEASST_Materials/Materials/","")
    mat2 = mat1.replace(".json", "")
    sorbent_mat_names.append(mat2)
    mat3 = mat2.replace("_clean", "")
    mat4 = mat3.replace("_h", "_clean_h")
    sorbent_mat_names_parsed.append(mat4)

#Get the list of all the sorbents that have the Atmosphere Check
list_atmo_files = glob.glob('/wrk/asm6/CSD_data/Atmosphere_check/results_300K_1atmN2*')
list_atmo_sorbents = []
for l in list_atmo_files:
    l1 = l.replace("/wrk/asm6/CSD_data/Atmosphere_check/results_300K_1atmN2_", "")
    l2 = l1.replace(".json", "")
    list_atmo_sorbents.append(l2)

##Get the sorbents that have completed Kh calculations:
#sorbent_Kh_results_files = glob.glob('/wrk/asm6/CSD_data/Results/results_300_N2*')
#sorbent_Kh_results = []
#for s in sorbent_Kh_results_files:
#    s1 = s[42:]
#    s2 = s1[:-5]
#    sorbent_Kh_results.append(s2)

##Get the list of sorbents that have the Atmosphere Check
#list_atmo_files = glob.glob('/wrk/asm6/CSD_data/Atmosphere_check/results_300K_1atmN2*')
##Get just the names of the sorbents
#list_atmo_sorbents = []
#for l in list_atmo_files:
#    l1 = l[56:] #strip off the directory part of the string
#    l2 = l1[:-5] #strip off the file type part of the string
#    list_atmo_sorbents.append(l2) #append to the list of sorbents

#remaining_sorbents = list(set(sorbent_Kh_results) - set(list_atmo_sorbents))

remaining_sorbents = list(set(sorbent_mat_names_parsed) - set(list_atmo_sorbents))



def Run_a_sorbent_calc(Sorbent):
    #Double check that the calc hasn't been done yet
    #Get the list of all the sorbents that have the Atmosphere Check
    list_atmo_files = glob.glob('/wrk/asm6/CSD_data/Atmosphere_check/results_300K_1atmN2*')
    list_atmo_sorbents = []
    for l in list_atmo_files:
        l1 = l.replace("/wrk/asm6/CSD_data/Atmosphere_check/results_300K_1atmN2_", "")
        l2 = l1.replace(".json", "")
        list_atmo_sorbents.append(l2)

    #Get the file name including any "_clean" suffix
    sorbent_mat_files = glob.glob('/wrk/asm6/CSD_data/CSD_FEASST_Materials/Materials/*.json')
    sorbent_mat_names = []
    sorbent_mat_names_parsed = []
    for mat in sorbent_mat_files:
        mat1 = mat.replace("/wrk/asm6/CSD_data/CSD_FEASST_Materials/Materials/","")
        mat2 = mat1.replace(".json", "")
        sorbent_mat_names.append(mat2)
        mat3 = mat2.replace("_clean", "")
        mat4 = mat3.replace("_h", "_clean_h")
        sorbent_mat_names_parsed.append(mat4)

    Sorbent_file_name = [elements for elements in sorbent_mat_names if(Sorbent in elements)][0]


    if Sorbent not in list_atmo_sorbents:
        Calc_args = Build_GCMC_json(300, 1.01325, Sorbent_file_name, ['N2_Trappe'], [1.])
    
    
        results_300 = run_gcmc(Calc_args)
    
        with open(f'/wrk/asm6/CSD_data/Atmosphere_check/results_300K_1atmN2_{Sorbent}.json', 'w') as f:
            json.dump(results_300, f, indent=4)
        

# Run_a_sorbent_calc(remaining_sorbents[0])

if __name__ == '__main__':
    
    
    p = multiprocessing.Pool(40)
    for _ in tqdm.tqdm(p.imap_unordered(Run_a_sorbent_calc, remaining_sorbents), total = len(remaining_sorbents)):
        pass
    