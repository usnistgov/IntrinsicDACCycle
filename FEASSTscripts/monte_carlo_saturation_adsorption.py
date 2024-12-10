import feasst as fst
import os
import numpy as np
import pandas as pd
import json

physical_constants = fst.CODATA2018()
kB = physical_constants.boltzmann_constant()  #J/K
Na = physical_constants.avogadro_constant()  #1/mol


def mc_setup(**kwargs):
    """
    Function that builds the MC Object and Adds Essential Functionality
    """
    
    #Unpack the arguments
    args = dict(kwargs)
    
    seed = args['seed']
    beta = 1. / args['temperature']
    box = args['box']
    tilts = args.get('tilts', [0., 0., 0.])
    cutoff = args['cutoff']
    
    # TMMC Run
    mc = fst.MonteCarlo()
    mc.set(fst.MakeRandomMT19937(fst.args({'seed': str(seed)})))
    
    # System
    adsorbate = args['adsorbate']['particle']
    adsorbent = args['adsorbent']['particle']
    
    # File existence check
    if not os.path.exists(adsorbate):
        raise Exception('Adsorbate particle file missing: ', adsorbate)
    if not os.path.exists(adsorbent):
        raise Exception('Adsorbent particle file missing: ', adsorbent)
    
    #physical_constants CODATA2018
    
    
    # System and Config Objects
    system = fst.System()
    config = fst.Configuration(
        fst.args({
            'side_length0': str(box[0]),
            'side_length1': str(box[1]),
            'side_length2': str(box[2]),
            'xy': str(tilts[0]),
            'xz': str(tilts[1]),
            'yz': str(tilts[2]),
            'particle_type0': adsorbate,
            'particle_type1': adsorbent,
            'add_particles_of_type1': str(1),
            'cutoff': str(cutoff)
        }))
    system.add(config)

    # Potential Groups
    if args['forcefield'] == 'LennardJones':
        system.add(
            fst.MakePotential(fst.MakeLennardJones())
        )
    elif args['forcefield'] == 'LennardJones_plus_Ewald':
        #Because we're looking for saturation, the repulsive r^-12 term will dominate. These settings are computationally convienent for that situation. 
        ewald_args = {
            'alpha': str(0.2), 
            'kxmax': str(5), 'kymax': str(5), 'kzmax': str(5)
        } 
        system.add(fst.MakePotential(fst.MakeEwald(fst.args(ewald_args))))
        TwoBodyFactory = fst.MakeModelTwoBodyFactory()
        TwoBodyFactory.add(fst.MakeLennardJones())
        TwoBodyFactory.add(fst.MakeChargeScreened())
        system.add(fst.MakePotential(TwoBodyFactory))
        system.add(
            fst.MakePotential(fst.MakeChargeScreenedIntra(),
                              fst.MakeVisitModelBond())
        )
        system.add(fst.MakePotential(fst.MakeChargeSelf()))
    else:
        raise Exception('Unknown forcefield', forcefield)
    
    mc.set(system)
    mc.set(fst.MakeMetropolis())
    
    # MC Moves
    mc.add(
        fst.MakeTrialTranslate(
            fst.args({
                'weight': '3.',
                'tunable_param': '2.',
                'particle_type': '0'
            })))
    mc.add(
            fst.MakeTrialRotate(
                fst.args({
                    'weight': '2.',
                    'tunable_param': '2.',
                    'particle_type': '0'
                })))
    mc.add(
        fst.MakeTrialTransfer(
            fst.args({
                'weight': '5.',
                'particle_type': '0'
            })))

    # Tune adjustable moves on the fly
    # mc.add(fst.MakeTune(fst.args({'steps_per': str(10000)})))
    mc.add(fst.MakeTune())
    
    # Basic Logging
    sorbent_name = args['adsorbent']['particle']
    sorbent_name = sorbent_name.replace("/users/asm6/DAC_data/CSD_FEASST_Materials/Materials/data.", "")
    mc.add(
        fst.MakeLog(
            fst.args({
                # 'steps_per': str(1.e4),
                # 'trials_per': str(1.e6),
                'trials_per': str(1.e8),
                'output_file': f'/tmp/asm6/Saturation/logfiles/{sorbent_name}_logfile.txt',
                'clear_file': 'true'
            })))
    
    mc.add(
        fst.MakeCheckEnergy(
            fst.args({
                # 'steps_per': str(1.e4),
                # 'trials_per': str(1.e4),
                # 'trials_per_write': str(1e6),
                # 'trials_per_update': str(1e6),
                'trials_per_write': str(1e8),
                'trials_per_update': str(1e8),
                'tolerance': str(1.e-6)
            })))
        


    return mc, ewald_args


def isotherm_step(mc,**kwargs):
    """
    Perform one step of MC to estimate the loading for a specified fugacity
    """
    
    #Unpack the arguments
    args = dict(kwargs)
    beta = 1. / args['temperature']
    fugacity = args['fugacity']
    eq_trials = args['eq_trials']
    prod_trials = args['prod_trials']
    step = args['step_number']

    mu_step = np.log(beta * fugacity)/beta
    print('Running Isotherm Step at mu = ', mu_step, step)

    init_loading = mc.configuration().num_particles() -1 # minus 1 for the adsorbent
    print('Initial Loading', init_loading)
    
    # Set the Thermodynamic Parameters
    mc.set(
        fst.MakeThermoParams(
            fst.args({
                'beta': str(beta),
                'chemical_potential0': str(mu_step)
            })))
    
    # Equilibration (burn-in)
    mc.attempt(eq_trials)
    print('Finished Equilibration / Burn-in; Activating <N> Accumulator')
    
    # Particle Counter
    sorbent_name = args['adsorbent']['particle']
    sorbent_name = sorbent_name.replace("/users/asm6/DAC_data/CSD_FEASST_Materials/Materials/data.", "")
    mc.add(
        fst.MakeNumParticles(
            fst.args({
                'particle_type': str(0),
                'output_file': f'/tmp/asm6/Saturation/histogram_files/{sorbent_name}_histogram_step'+str(step)+'.txt'
            })))

    # Production Run
    mc.attempt(prod_trials)
    print('Finished Production Run')

    final_loading = mc.configuration().num_particles() -1 # minus 1 for the adsorbent
    print('Final Loading', final_loading)
    
    # Save and return the average particle count
    analyzer = mc.analyze(mc.num_analyzers() - 1)
    Navg = analyzer.accumulator().average()
    Nstdev = analyzer.accumulator().block_stdev()

    # Remove the Particle Counter (to prepare for the next step)
    #  MAJOR NOTE: in older FEASST, there is no way to remove an Analyzer
    #mc.run(fst.MakeRemoveStepper(
    #    fst.args({
    #        'name': 'NumParticles'
    #    })
    #))
    
    return Navg, Nstdev, mc, mu_step, final_loading


def estimate_max_loading(threshold=0.01,**kwargs):
    """
    Perform a loop over fugacities to estimate the isotherm.
    Stop loop when the fractional step change is smaller than the specified threshold
    """
    
    # Unpack Arguments
    args = dict(kwargs)

    # Build the MC Object
    mc_runner, ewald_args = mc_setup(**args)
    print('Object Setup Complete')

    # List for storing Isotherm Points
    isotherm = []
    isotherm_dict = {}

    #  First Isotherm Step
    args['fugacity'] = 0.001
    args['step_number'] = 1
    loading, stdev, mc_runner, mu, particles = isotherm_step(mc_runner, **args)
    isotherm.append([args['fugacity'], loading, stdev])

    isotherm_dict[1] = {"step_number":args["step_number"],
                        "fugacity": args["fugacity"],
                        "mu": mu,
                        "loading": loading,
                        "stdev": stdev,
                        "final_particles": particles}

    # Loop through fugacities
    finished = False
    while not finished:
        args['fugacity'] += 0.001
        args['step_number'] += 1
        loading, stdev, mc_runner, mu, particles = isotherm_step(mc_runner, **args)

        # Check for Convergence of the loading
        last_loading = isotherm[-1][1]
        print('Current Loading', loading)
        print('Last Loading', last_loading)
        if last_loading != 0.0:
            print('Change since last', np.abs((loading-last_loading)/last_loading))
            if np.abs((loading-last_loading)/last_loading) < threshold:
                finished = True
        print(args['fugacity'], loading, stdev)
        isotherm.append([args['fugacity'], loading, stdev])
        attempted_trials = mc_runner.get_trial_factory().num_attempts()
        succesful_trials = mc_runner.get_trial_factory().num_success()
        isotherm_dict[args["step_number"]] = {"step_number":args["step_number"],
                                                "fugacity": args["fugacity"],
                                                "mu": mu,
                                                "loading": loading,
                                                "stdev": stdev,
                                                "particles": particles,
                                                "cumulative_attempted_trials": attempted_trials,
                                                "cumulative_succesful_trials": succesful_trials}
        print()
        #finished = True # break
        if args['step_number'] > 100:
            finished = True
            isotherm_dict["terminated_early"] = True

    print('finished Isotherm')
    #for line in isotherm:
    #    print(line)

    return isotherm, args, mc_runner, isotherm_dict, ewald_args


def Build_saturation_calc_json(Tref, Sorbent, Sorbate):
    #Tref is the refernce temperature in K
    #Sorbent is a string of the name of the Sorbent
    #Sorbate is a string of the name of the Sorbate 
    
    adsorbent = f'/users/asm6/DAC_data/CSD_FEASST_Materials/Materials/data.{Sorbent}'
    adsorbate = f'data.{Sorbate}'
    
    # Get necessary info out of the material data file
    material_data_json = f'/users/asm6/DAC_data/CSD_FEASST_Materials/Materials/{Sorbent}.json'
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
    args = {
        'seed': 36897,
        'temperature': Tref * kB * Na / 1000., # kJ/mol
        'cutoff': 12.0, #Angstroms

        'box': box,
        'tilts': tilts,

        'adsorbate': {
            'particle': adsorbate
        },
        'adsorbent': {
            'particle': adsorbent
        },
        'forcefield': 'LennardJones_plus_Ewald'
    }

    # Set the number of MC Trials for each isotherm step
    args['eq_trials'] = 5000
    args['prod_trials'] = 5000

    return args