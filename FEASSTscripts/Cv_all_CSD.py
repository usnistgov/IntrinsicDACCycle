# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 11:35:48 2022

@author: Austin McDannald
"""

import multiprocessing
import tqdm
import glob
import json

from xgboost_predict_cv import cif_to_cv

list_of_sorbents = glob.glob('/wrk/asm6/CSD_data/CSD_FEASST_Materials/CIF_files/*.cif')




temperatures=[250.00,
              275.00,
              300.00,
              325.00,
              350.00,
              375.00,
              400.00]

models = "/wrk/asm6/CSD_data/Cv_predictions/ensemble_models_smallML_120_10"

def Run_a_Cv_calc(sorbent):
    cif_to_cv(sorbent, models, temperatures)

# for sorb in list_of_sorbents:
#     Run_a_Cv_calc(sorb)


# Run_a_Cv_calc(list_of_sorbents[0])
if __name__ == '__main__':
        
    p = multiprocessing.Pool(20)
    for _ in tqdm.tqdm(p.imap_unordered(Run_a_Cv_calc, list_of_sorbents), total = len(list_of_sorbents)):
        pass