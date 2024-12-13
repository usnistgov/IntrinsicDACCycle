# -*- coding: utf-8 -*-
"""
Created on Thu May 12 14:26:44 2022

@author: Austin McDannald
"""

import sys
import os
import json
import argparse

import numpy as np

import feasst as fst

physical_constants = fst.CODATA2018()
kB = physical_constants.boltzmann_constant()  #J/K
Na = physical_constants.avogadro_constant()  #1/mol


def Build_GCMC_json(temperature, target_p, Sorbent, Sorbates, y):
    #Temperature is the refernce temperature in K
    #target_p is the total pressure in bar
    #Sorbent is a string of the name of the Sorbent
    #Sorbate is a list of strings of the name of the Sorbate gasses 
    #y is the list of mole fractions of each gas species
    if len(Sorbates) == 1:
        y = [1.]
    
    #convert temperature to KJ/mol
    temperature *= kB * Na / 1000.  # <- kJ/mol
    # Estimate betamu vector from Ideal Gas EOS
    target_p *= 1.e5 * 1.e-30 * Na / 1000.  # units: kJ/mol/ang^3
    betamu = [np.log(yi * target_p / temperature) for yi in y]
    
    adsorbent = f'/wrk/asm6/CSD_data/CSD_FEASST_Materials/Materials/data.{Sorbent}'
    
    adsorptives = []
    for Sorbate in Sorbates:
        adsorbate = f'data.{Sorbate}'
        adsorptives.append(adsorbate)
    
    # Get necessary info out of the material data file
    material_data_json = f'/wrk/asm6/CSD_data/CSD_FEASST_Materials/Materials/{Sorbent}.json'
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
    args_GCMC = {
        'temperature': temperature, #in KJ/mol
        'target_p': target_p, #in KJ/mol/ang^3
        'betamu': betamu,
        'beta': str(1. / temperature),
        'steps_per': int(1e4),
        'seed': str(92532459),  # use a consistent seed when debugging)}
        'trials_EQ': int(1.e6), #Equilibration Burn-in trials
        'trials_PROD': int(5.e6), #Production Run (data collected)
        # 'trials_EQ': int(1.e4), #Equilibration Burn-in trials
        # 'trials_PROD': int(2.e4), #Production Run (data collected)       
        'box': box,
        'tilts': tilts,
        'adsorptives': adsorptives,
        'adsorbent': adsorbent,
        'ff_type': 'ewald',
        'tail_type': 'FS', 
        }
        
    return args_GCMC

def DL_POLY_ewald(rcut, boxsize, tolerance):
    """ Parameters from the DL_POLY Algorithm
    https://doi.org/10.1080/002689798167881
    """
    # DL_POLY Algorithm
    eps = min(tolerance, 0.5)
    xi = np.sqrt(np.abs(np.log(eps * rcut)))
    alpha = np.sqrt(np.abs(np.log(eps * rcut * xi))) / rcut
    chi = np.sqrt(-np.log(eps * rcut * ((2. * xi * alpha)**2)))
    kmax = [int(0.25 + boxsize[i] * alpha * chi / np.pi) for i in range(3)]
    return alpha, kmax



def run_gcmc(GCMC_args):
    """
    Set up and run a very simple GCMC simulation
    """
    trials_EQ = GCMC_args['trials_EQ']
    trials_PROD = GCMC_args['trials_PROD']
    steps_per = GCMC_args['steps_per']
    ff_type = GCMC_args['ff_type']
    tail_type = GCMC_args['tail_type']
    
    thermo_params = {'beta': GCMC_args['beta']}
    temperature = GCMC_args['temperature']
    betamu = GCMC_args['betamu']
    
    adsorptives = GCMC_args['adsorptives']
    adsorbent = GCMC_args['adsorbent']
    box = GCMC_args['box']
    tilts = GCMC_args['tilts']
    
    mc = fst.MakeMonteCarlo()
    mc.set(fst.MakeRandomMT19937(fst.args({'seed': GCMC_args['seed']})))
    
    config_args = dict()
    index = 0

    for part in adsorptives:
        config_args['particle_type' + str(index)] = part
        thermo_params['chemical_potential' + str(index)] = str(betamu[index] *
                                                               temperature)
        index += 1

    # Add the MOF to the config args
    
    config_args['particle_type' +
                str(len(adsorptives))] = adsorbent  # remember, count from zero

    # Instantiate the Configuration Object
    config = fst.MakeConfiguration(
        fst.args(
            dict(
                {
                    'side_length0': str(box[0]),
                    'side_length1': str(box[1]),
                    'side_length2': str(box[2]),
                    'xy': str(tilts[0]),
                    'xz': str(tilts[1]),
                    'yz': str(tilts[2]),
                }, **config_args)))

    # Add 1 MOF particle to the system
    config.add_particle_of_type(len(adsorptives))  # remember, count from zero

    # Add the Configuration Object
    mc.add(config)

    # Assemble Forcefield

    if ff_type == 'ewald':
        #-----------------------------------
        # Ewald plus LJ
        #
        # Collect or Set the Ewald Parameters
        rcut = max([
            config.model_params().select('cutoff').value(i)
            for i in range(config.model_params().select('cutoff').size())
        ])
        print(rcut)
        alpha, kmax = DL_POLY_ewald(rcut, box, 1.e-5)
        print('kmax', kmax)
        ewald_args = {
            'alpha': str(alpha),
            'kxmax': str(kmax[0]),
            'kymax': str(kmax[1]),
            'kzmax': str(kmax[2])
        }
        # Ewald Fourier-space energy
        mc.add(fst.MakePotential(fst.MakeEwald(fst.args(ewald_args))))

        # Combine the LJ and Realspace Coulomb in a single loop
        TwoBodyFactory = fst.MakeModelTwoBodyFactory()
        if tail_type == 'CS':
            TwoBodyFactory.add(fst.MakeLennardJonesCutShift())
        elif tail_type == 'FS':
            TwoBodyFactory.add(fst.MakeLennardJonesForceShift())
        elif tail_type == 'LRC':
            TwoBodyFactory.add(fst.MakeLennardJones())
        else:
            raise Exception('Unknown Tail Correction:', tail_type)
        TwoBodyFactory.add(fst.MakeChargeScreened())
        mc.add(fst.MakePotential(TwoBodyFactory))

        # # If requested, add Long-range Corrections
        # if tail_type == 'LRC':
        #     self.add(fst.MakePotential(fst.MakeLongRangeCorrections()))

        # Correct to remove spurious intramolecular Coulomb
        #   VisitModelBond -> means perform this energy calculation on bonded interactions
        mc.add(
            fst.MakePotential(fst.MakeChargeScreenedIntra(),
                              fst.MakeVisitModelBond()))

        # Ewald Fourier-space self-energy
        mc.add(fst.MakePotential(fst.MakeChargeSelf()))
        #--------------------
    elif ff_type == 'LJ':
        #--------------------
        # Plain Lennard-Jones Version
        mc.add(fst.MakePotential(fst.MakeLennardJones()))
        #--------------------
    else:
        raise Exception('ERROR: unknown force field', ff_type)

    # Metropolis Strategy
    mc.set(fst.MakeThermoParams(thermo_params))
    mc.set(fst.MakeMetropolis())  # need this to set criteria

    # Monte Carlo Trial Moves
    for particle_type in range(mc.configuration().num_particle_types() - 1):
        # TrialTranslate is universal
        mc.add(
            fst.MakeTrialTranslate(
                fst.args({
                    'particle_type': str(particle_type),
                    'weight': '1.',
                    'tunable_param': '1.'
                })))
        # Decide between DCCB and vanilla moves
        if mc.configuration().particle_type(particle_type).num_sites() > 1:
            mc.add(
                fst.MakeTrialRotate(
                    fst.args({
                        'particle_type': str(particle_type),
                        'weight': '1.',
                        'tunable_param': '1.'
                    })))
        mc.add(
            fst.MakeTrialTransfer(
                fst.args({
                    'particle_type': str(particle_type),
                    'weight': '4'
                })))

    # Housekeeping Operations
    mc.add(
        fst.MakeCheckEnergy(
            fst.args({
                'steps_per': str(steps_per),
                'tolerance': '0.0001'
            })))
    mc.add(
        fst.MakeTune(
            fst.args({
                'steps_per': str(steps_per),
                #"stop_after_phase": "0"
            })))

    # Analysis Routines
    #  Basic Sim log
    mc.add(
        fst.MakeLog(
            fst.args({
                'steps_per': str(steps_per),
                'file_name': 'clones' + '_log.txt',
                #"file_name_append_phase": "True",
                'clear_file': 'True'
            })))

    # Equilibration Run
    mc.attempt(trials_EQ)

    # Assign Analyzers for Production Run
    #  Energy Tracker
    mc.add(
        fst.MakeEnergy(
            fst.args({
                'file_name': 'en' + '.txt',
                #"file_name_append_phase": "True",
                #"start_after_phase": "0",
                'steps_per_write': str(steps_per),
                'steps_per_update': '1',
                'multistate': 'True'
            })))
    #  Particle Counter
    if mc.configuration().num_particle_types() > 1:
        for particle in range(mc.configuration().num_particle_types() - 1):
            mc.add(
                fst.MakeNumParticles(
                    fst.args({
                        'particle_type':
                        str(particle),
                        'file_name':
                        'num' + str(particle) + '_' + '.txt',
                        #"file_name_append_phase": "True",
                        #"start_after_phase": "0",
                        'steps_per_write':
                        str(steps_per),
                        'steps_per_update':
                        '1',
                        'multistate':
                        'True'
                    })))



    mc.attempt(trials_PROD)
    
    #Extract the results from the MC object (SWiG objects)
    #Energy
    energy_analyzer = mc.analyze(1) 
    E_average = energy_analyzer.analyze(0).accumulator().average()
    E_std = energy_analyzer.analyze(0).accumulator().stdev()
    E_bstd = energy_analyzer.analyze(0).accumulator().block_stdev()
    
    GCMC_args['Energy_average'] = E_average
    GCMC_args['Energy_std'] = E_std
    GCMC_args['Energy_block_std'] = E_bstd
    
    #Number of molecules of each adsorptive:
    for particle in range(mc.configuration().num_particle_types() - 1):
    # subtract 1 because the the "last" particle is the MOF
        num_analyzer = mc.analyze(2+particle)
        num_average = num_analyzer.analyze(0).accumulator().average() #average number of particles
        num_std = num_analyzer.analyze(0).accumulator().stdev()
        num_bstd = num_analyzer.analyze(0).accumulator().block_stdev()
        
        GCMC_args[f'molecules_adsorptives{particle}_average'] = num_average
        GCMC_args[f'molecules_adsorptives{particle}_std'] = num_std
        GCMC_args[f'molecules_adsorptives{particle}_block_std'] = num_bstd
    
    return GCMC_args
    
# run_gcmc(0)




