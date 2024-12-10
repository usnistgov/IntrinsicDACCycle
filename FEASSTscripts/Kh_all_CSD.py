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

from Set_up_Kh_calc import Build_calc_json
from porous_networks_utils.porous_network_utils import atomistic_system_generator, henry_AlwaysReject

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


#Get the list of all the sorbents that have the Henry Constant Calculation
sorbent_results_files = glob.glob('/wrk/asm6/CSD_data/Results/results_300_N2*')
sorbent_results = []
for l in sorbent_results_files:
    l1 = l.replace("/wrk/asm6/CSD_data/Results/results_300_N2_", "")
    l2 = l1.replace(".json", "")
    sorbent_results.append(l2)



##Read all the sorbent files
##list_of_sorbent_files = glob.glob('/wrk/asm6/CSD_data/CSD_FEASST_Materials/Materials/*_clean.json')
#list_of_sorbent_files = glob.glob('/wrk/asm6/CSD_data/CSD_FEASST_Materials/Materials/*.json') #all the .json's

##Get just the names of the sorbents
#list_of_sorbents = []
#for l in list_of_sorbent_files:
#    l1 = l[50:] #strip off the directory part of the string
#    l2 = l1.replace(".json", "") #strip off the file type part of the string
#    if "_clean_h" in l2:
#    	l3 = l2
#    else:
#    	l3 = l2.replace("_clean", "") #strip off the "clean" tag part of the string
#    list_of_sorbents.append(l3) #append to the list of sorbents

##Get the sorbents that have been completed:
#sorbent_results_files = glob.glob('/wrk/asm6/CSD_data/Results/results_300_N2*')
#sorbent_results = []
#for s in sorbent_results_files:
#    s1 = s[42:]
#    s2 = s1[:-5]
#    sorbent_results.append(s2)

remaining_sorbents = list(set(sorbent_mat_names_parsed) - set(sorbent_results))




def Run_a_sorbent_calc(Sorbent):
    #Double check that the calc hasn't been done yet
    #Get the list of all the sorbents that have the Henry Constant Calculation
    sorbent_results_files = glob.glob('/wrk/asm6/CSD_data/Results/results_300_N2*')
    sorbent_results = []
    for l in sorbent_results_files:
        l1 = l.replace("/wrk/asm6/CSD_data/Results/results_300_N2_", "")
        l2 = l1.replace(".json", "")
        sorbent_results.append(l2)

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

    if Sorbent not in sorbent_results:

        Calc_args = Build_calc_json(300., Sorbent_file_name,'N2_Trappe')
    
        results_300 = henry_AlwaysReject(**Calc_args)
        del results_300['system_generator_class']
    
        with open(f'/wrk/asm6/CSD_data/Results/results_300_N2_{Sorbent}.json', 'w') as f:
            json.dump(results_300, f, indent=4)

if __name__ == '__main__':
    
    
    p = multiprocessing.Pool(40)
    for _ in tqdm.tqdm(p.imap_unordered(Run_a_sorbent_calc, remaining_sorbents), total = len(remaining_sorbents)):
        pass
    