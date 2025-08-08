# -*- coding: utf-8 -*-
"""
Created on Thu May 12 14:26:44 2022

@author: Austin McDannald
"""

import sys
import numpy as np
import os
import json
import copy

import feasst as fst

from porous_network_utils import atomistic_system_generator, henry_AlwaysReject

# CODATA Physical Constants
physical_constants = fst.CODATA2018()
kB = physical_constants.boltzmann_constant()  #J/K
Na = physical_constants.avogadro_constant()  #1/mol
amu = 1.66053906660e-27 #kg

def uncertainty(results,prop_type='Kh'):
    """
    Computes standard uncertainty in MC-integration-sampled property
    results = dictionary from HenryAlwaysReject
    prop_type = either Kh (default) or Qst
    returns the standard uncertainty (already divided by the degrees of freedom)
    """
    
    if prop_type == 'Kh':
        var = results['Kcoeff_var'][0]
    elif prop_type == 'Qst':
        var = ((-results['Kcoeff'][1]/(results['Kcoeff'][0]**2))**2) * results['Kcoeff_var'][0]
        var += ((1./results['Kcoeff'][0])**2) * results['Kcoeff_var'][1]
        var += ((-results['Kcoeff'][1]/(results['Kcoeff'][0]**2))) * ((1./results['Kcoeff'][0])) * 2. * results['Kcoeff_covar'][0][1]
    
    return np.sqrt(var/ float(results['trials']))

def extrap_uncertainty(results,dbeta,prop_type='Kh'):
    """
    Computes standard uncertainty in MC-integration-sampled property
    results = dictionary from HenryAlwaysReject
    dbeta = (beta - beta_0) in same units as results
    prop_type = either Kh (default) or Qst
    returns the standard uncertainty (already divided by the degrees of freedom)
    """
    
    num_coeffs = len(results['Kcoeff'])
    #num_coeffs = 5
    if prop_type == 'Kh':
        var = results['Kcoeff_var'][0]
        # jth sensitivity coefficient is (dbeta**j)
        # Assume Zero Covariance
        #var = sum([(dbeta**(2*j))*u2j for j, u2j in enumerate(results['Kcoeff_var'][0:num_coeffs]) ])
        # Propagate Covariance
        var = (
            sum(map(
                sum,[[(dbeta**(j+l))*results['Kcoeff_covar'][j][l] for l in range(num_coeffs)] for j in range(num_coeffs)]
            ))
        ) 
    elif prop_type == 'Qst':
        # pre-calculate the sensitivity coefficients
        if dbeta == 0.:
            sens_coeff = [0.]*num_coeffs
            sens_coeff[0] = -results['Kcoeff'][1]/(results['Kcoeff'][0]**2)
            sens_coeff[1] = 1./results['Kcoeff'][0]
        else:
            sens_coeff = [
                float(j)*(dbeta**(j-1)) / sum([results['Kcoeff'][j]*(dbeta**j) for j in range(0,num_coeffs)])
                - sum([float(j)*results['Kcoeff'][j]*(dbeta**(j-1)) for j in range(1,num_coeffs)]) * (dbeta**j)/
                sum([results['Kcoeff'][j]*(dbeta**j) for j in range(0,num_coeffs)])**2
                for j in range(num_coeffs)
            ]
        # Assume Zero Covariance
        #var = sum([
        #    (sens_coeff[j]**2)*u2j for j, u2j in enumerate(results['Kcoeff_var'][0:num_coeffs])
        #])
        # Propagate Covariance
        var = (
            sum(map(
                sum,[[sens_coeff[l]*sens_coeff[j]*results['Kcoeff_covar'][j][l] for l in range(num_coeffs)] for j in range(num_coeffs)]
            ))
        ) 
        
    return np.sqrt(var/ float(results['trials']))


def Build_calc_json(Tref, Sorbent, Sorbate):
    #Tref is the refernce temperature in K
    #Sorbent is a string of the name of the Sorbent
    #Sorbate is a string of the name of the Sorbate 
    
#    adsorbent = f'/users/asm6/DAC_data/CSD_FEASST_Materials/Materials/data.{Sorbent}'
    adsorbate = f'data.{Sorbate}'
    adsorbent = f'data.{Sorbent}'
    
    
    # Get necessary info out of the material data file
#    material_data_json = f'/users/asm6/DAC_data/CSD_FEASST_Materials/Materials/{Sorbent}.json'
    material_data_json = f'{Sorbent}.json'
    material_data = json.load(open(material_data_json,mode='r'))
    
    box = [ 
        material_data['cell_matrix'][0][0],
        material_data['cell_matrix'][1][1],
        material_data['cell_matrix'][2][2]
    ]
    tilts = [
        material_data['cell_matrix'][0][1],
        material_data['cell_matrix'][0][2],
        material_data['cell_matrix'][1][2]
    ]
    
    
    # Build the dictionary
    
    # "Fast" calculation
    args_AlwaysReject = {
        'beta': (1. / Tref) * (1. / kB) * (1. / Na) * 1000.,
        'trials': int(1.e6),  #will be overwritten if auto_converge=True
        'seed': 1200512782,
        'debug': True, #verbose output
        'scale_factor': 10., #helps avoid numerical overflow
        'ncoeffs': 21,
        'system_generator_class': atomistic_system_generator,
        'system_definition': {
            'adsorbate': {
                'particle': adsorbate
            },
            'adsorbent': {
                'particle': adsorbent
            },
            'forcefield_type': 'LennardJones_plus_Ewald',
            'tail_type': 'FS',
            'Ewald_Parameters': {
                'tolerance': 1.e-5,
                'method': 'DL_POLY'
            },
            'box': box,
            'tilts': tilts,
        },
    
        'auto_converge': True,
        'convergence_crit': {
            'trials_per_test': 10000,  # trials between tests
            'min_trials': 1000000,       # minimum trials
            'rel_unc_0': 0.01,             # maximum relative standard error in Kcoeff[0]
            'key': 10,                     # index of coefficient to apply convergence
            'rel_unc_key': 0.05,           # maximum relative standard error in Kcoeff[key]
            'max_trials': 50000000        # maximum trials (overrides 'trials' key)
        }
    
    }
    return args_AlwaysReject
