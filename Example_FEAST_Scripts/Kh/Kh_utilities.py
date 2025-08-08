# -*- coding: utf-8 -*-
"""Build FEASST System for further operations"""
# pylint: disable-msg=invalid-name   #because I use snake_case
# pylint: disable-msg=too-many-locals

import os
import math
import copy
import json
import numpy as np
import feasst as fst

from set_ewald import set_ewald


class atomistic_system_generator(fst.System):
    """Class, derived from feasst.System() to generate the atomistic confinement system"""
    def __init__(self, **kwargs):
        """Add the Config() object to the derived system object"""
        super().__init__()

        args = dict(kwargs)
        adsorbate = args['adsorbate'][
            'particle']  #NOTE: using FEASST terminology; particle = molecule
        adsorbent = args['adsorbent'][
            'particle']  #NOTE: using FEASST terminology; particle = molecule
        box = args['box']
        tilts = args.get('tilts', [0., 0., 0.])
        adsorbate_initial_fill = args['adsorbate_initial_fill']

        # File existence check
        if not os.path.exists(adsorbate):
            raise FileNotFoundError('Adsorbate particle file missing: ',
                                    adsorbate)
        if not os.path.exists(adsorbent):
            raise FileNotFoundError('Adsorbent particle file missing: ',
                                    adsorbent)

        # FEASST configuration object
        config = fst.Configuration(
            fst.args({
                'side_length0': str(box[0]),
                'side_length1': str(box[1]),
                'side_length2': str(box[2]),
                'xy': str(tilts[0]),
                'xz': str(tilts[1]),
                'yz': str(tilts[2]),
                'particle_type0': adsorbate,
                'particle_type1': adsorbent
            }))
        num0 = adsorbate_initial_fill
        #num0 = 1  #Henry Constant MC Integration using select particle
        #num0 = 0  #Henry Constant MC Integration using AlwaysReject
        #num0 = 0  #GC Simulation
        num1 = 1  #Adsorbent solid is "one" particle
        for _ in range(num0):
            config.add_particle_of_type(0)
        for _ in range(num1):
            config.add_particle_of_type(1)
        self.add(config)

        # Add Potential Groups
        self.set_potential_group(args)

    def set_potential_group(self, args):
        """Generalized function to specify the potential forcefield"""
        # pylint: disable-msg=too-many-branches
        forcefield = args['forcefield_type']
        tail_type = args['tail_type']

        # Add the specified forcefield
        if forcefield == 'LennardJones':
            if tail_type == 'CS':
                self.add(fst.MakePotential(fst.MakeLennardJonesCutShift()))
            elif tail_type == 'FS':
                self.add(fst.MakePotential(fst.MakeLennardJonesForceShift()))
            elif tail_type == 'LRC':
                self.add(fst.MakePotential(fst.MakeLennardJones()))
                self.add(fst.MakePotential(fst.MakeLongRangeCorrections()))
            elif tail_type == 'CUT':
                self.add(fst.MakePotential(fst.MakeLennardJones()))
            else:
                raise AttributeError('Unknown Tail Correction:', tail_type)
        elif forcefield == 'LennardJones_plus_Ewald':
            # Collect or Set the Ewald Parameters
            ewald_args = self.ewald_arguments(args)
            args['Ewald_Parameters']['set'] = ewald_args

            # Ewald Fourier-space energy
            self.add(fst.MakePotential(fst.MakeEwald(fst.args(ewald_args))))

            # Combine the LJ and Realspace Coulomb in a single loop
            TwoBodyFactory = fst.MakeModelTwoBodyFactory()
            if tail_type == 'CS':
                TwoBodyFactory.add(fst.MakeLennardJonesCutShift())
            elif tail_type == 'FS':
                TwoBodyFactory.add(fst.MakeLennardJonesForceShift())
            elif tail_type == 'LRC':
                TwoBodyFactory.add(fst.MakeLennardJones())
            elif tail_type == 'CUT':
                TwoBodyFactory.add(fst.MakeLennardJones())
            else:
                raise AttributeError('Unknown Tail Correction:', tail_type)
            TwoBodyFactory.add(fst.MakeChargeScreened())

            # Cell Lists are based on the maximum cutoff
            rcut = max([  #pylint: disable=consider-using-generator
                self.configuration().model_params().select('cutoff').value(i)
                for i in range(self.configuration().model_params().select(
                    'cutoff').size())
            ])
            if self.configuration().domain().is_tilted():
                # override cell list for non-cuboid domain
                num_cells = [1, 1, 1]
            else:
                num_cells = [
                    math.trunc(self.configuration().domain().side_length(i) /
                               rcut) for i in range(3)
                ]
            if min(num_cells) >= 3 and max(num_cells) > 3:
                print('Activating cell list for full potential')
                self.add(
                    fst.MakePotential(
                        TwoBodyFactory,
                        fst.MakeVisitModelCell(
                            fst.args({'min_length': str(rcut)}))))
            else:
                print('Disabling cell list for full potential')
                self.add(fst.MakePotential(TwoBodyFactory))

            # If requested, add Long-range Corrections
            if tail_type == 'LRC':
                self.add(fst.MakePotential(fst.MakeLongRangeCorrections()))

            # Correct to remove spurious intramolecular Coulomb
            #   VisitModelBond -> means perform this energy calculation on bonded interactions
            self.add(
                fst.MakePotential(fst.MakeChargeScreenedIntra(),
                                  fst.MakeVisitModelBond()))

            # Ewald Fourier-space self-energy
            self.add(fst.MakePotential(fst.MakeChargeSelf()))

            # Null the energy of the N=0 macrostate by subtracting the energy of the empty MOF
            # print(self.energy(),
            #       'current energy, in generator class, before nulling')
            # empty_energy = self.energy()
            # print('Nulling Background Energy')
            # self.add(
            #     fst.MakePotential(
            #         fst.MakeBackground(
            #             fst.args({'constant': str(-empty_energy)}))))
            # print(self.energy(),
            #       'current energy, in generator class, after nulling')

        else:
            raise AttributeError('ERROR: Unknown forcefield type: ' +
                                 forcefield)

    def ewald_arguments(self, args):
        """
        Parse the system_definition dictionary and either set
        or calculate the Ewald parameters
        """

        if 'tolerance' in args['Ewald_Parameters']:
            # Set parameters based on tolerance and a selected method
            method = args['Ewald_Parameters']['method']
            tolerance = args['Ewald_Parameters']['tolerance']
            alpha, kmax = set_ewald(method,
                                    self,
                                    tolerance=tolerance,
                                    minmax='min')
            ewald_args = {
                'alpha': str(alpha),
                'kxmax': str(kmax[0]),
                'kymax': str(kmax[1]),
                'kzmax': str(kmax[2])
            }
        else:
            # Fall back to direct specification
            if 'alpha' in args['Ewald_Parameters']:
                # Default is direct specification of Ewald damping parameter (alpha)
                alpha = args['Ewald_Parameters']['alpha']
            else:
                # Alternate option is input of alpha*min(Lx,Ly,Lz)
                alpha = args['Ewald_Parameters']['alpha_ratio']
                box = [
                    self.configuration().domain().side_length(i)
                    for i in range(3)
                ]
                alpha /= min(box)

            if 'kmax' in args['Ewald_Parameters']:
                # Default is direct specification of the number of k-vectors
                kmax = args['Ewald_Parameters']['kmax']
                ewald_args = {
                    'alpha': str(alpha),
                    'kxmax': str(kmax[0]),
                    'kymax': str(kmax[1]),
                    'kzmax': str(kmax[2])
                }
            else:
                # Alternate option is input of max(k^2)
                ewald_args = {
                    'alpha': str(alpha),
                    'kmax_squared':
                    str(args['Ewald_Parameters']['kmax_squared'])
                }

        return ewald_args


def henry_AlwaysReject(**kwargs):
    """Compute the temperature-extrapolation coefficients for the Henry Constant"""
    #Unpack arguments
    args = dict(kwargs)
    seed = args['seed']
    beta = args['beta']
    trials = args['trials']
    scale_factor = args['scale_factor']
    if scale_factor > 1.:
        scale_factor = scale_factor**(-1.)
    ncoeffs = args['ncoeffs']
    debug = args['debug']

    # Auto convergence parameters
    auto_converge = args.get('auto_converge', False)  # fall-back to false
    if auto_converge:
        convergence_crit = args['convergence_crit']
        trials = convergence_crit[
            'max_trials']  # overwrite trials with max_trials
        converge_key = convergence_crit[
            'key']  # which coefficient for convergence evaluation

    #Assemble the FEASST system object
    system_generator = args['system_generator_class']
    generator_args = args['system_definition']
    generator_args[
        'adsorbate_initial_fill'] = 0  #Note the difference here. "Select" approach vs. AlwaysReject
    # temporary parameter, until Kh and TMMC codes can use same system object
    system = system_generator(**generator_args)
    #write the ensemble parameters to disk
    ensemble_dict = {'beta': beta}
    if 'box' in args['system_definition']:
        ensemble_dict['box']: args['system_definition']['box']
    if 'forcefield_type' in args['system_definition']:
        if args['system_definition'][
                'forcefield_type'] == 'LennardJones_plus_Ewald':
            ensemble_dict['Ewald_Parameters'] = args['system_definition'][
                'Ewald_Parameters']
    with open('ensemble.json', mode='w') as f:
        json.dump(ensemble_dict, f, sort_keys=True, indent=4)

    #Monte Carlo Integration
    mc = fst.MonteCarlo()
    mc.set(fst.MakeRandomMT19937(fst.args({'seed': str(seed)})))
    mc.set(system)
    mc.set(
        fst.MakeThermoParams(
            fst.args({
                'beta': str(beta),
                'chemical_potential0': '1'
            })))
    mc.set(fst.MakeAlwaysReject())
    #MC Integration is done via a series of test insertions:
    # 1. Attempt an insertion at a random position, with random orientation
    # 2. Calculate the deltaE for the test insertion
    # 3. "Reject" the insertion to preserve the bare adsorbent state
    trial_move = fst.MakeTrialAdd(
        fst.args({
            'particle_type': '0',
            'new_only': 'true'
        }))
    mc.add(trial_move)
    #mc.add(fst.MakeLogAndMovie(fst.args({"steps_per": str(int(1e4)), "file_name": "tutorial_1"})))
    henry_index = mc.num_analyzers()  # pylint: disable-msg=unused-variable
    mc.add(fst.MakeHenryCoefficient())

    # Bare adsorbent energy
    # (system is initialized in an empty state)
    energy_bare = mc.system().energy()
    #print('bare energy', energy_bare)

    # Accumulators, Running Averages, Pre-calculated functions
    counter = 0
    Kcoeff = [0.] * ncoeffs
    Kcoeff_covar = [[0. for l in range(ncoeffs)] for j in range(ncoeffs)]
    factorial = [math.factorial(j) for j in range(ncoeffs)]
    values = [0.] * ncoeffs

    for _ in range(trials):
        mc.attempt(1)
        #deltaU = trial_move.accept().energy_new()
        deltaU = trial_move.accept().energy_new() - energy_bare

        #print(deltaU)
        #raise Exception("BREAK")
        counter += 1

        # Compute all of the coefficients
        #  Current values are computed first to simplify covariance accumulation
        if math.isnan(deltaU) or deltaU > 1.e19:
            values = [0.] * ncoeffs
        else:
            for j in range(ncoeffs):
                values[j] = np.exp(-beta * deltaU) * (
                    deltaU * scale_factor)**j * ((-1.)**j) / factorial[j]

        # Update running averages
        for j in range(ncoeffs):
            Kcoeff[j] += (values[j] - Kcoeff[j]) / float(counter)
            for l in range(ncoeffs):
                Kcoeff_covar[j][l] += (values[j] * values[l] -
                                       Kcoeff_covar[j][l]) / float(counter)

        # Test Auto-convergence
        if auto_converge is True and counter >= convergence_crit['min_trials']:
            if counter % convergence_crit['trials_per_test'] == 0:
                # criteria applied to K0(beta)
                var0 = float(counter) * (Kcoeff_covar[0][0] -
                                         Kcoeff[0]**2) / float(counter - 1)
                rel_unc_0 = np.sqrt(var0 / float(counter)) / Kcoeff[0]
                # criteria applied to Kj(beta) for the Jth coefficient
                var_key = float(counter) * (
                    Kcoeff_covar[converge_key][converge_key] -
                    Kcoeff[converge_key]**2) / float(counter - 1)
                rel_unc_key = np.sqrt(
                    var_key / float(counter)) / Kcoeff[converge_key]
                print('attempting auto converge', rel_unc_0, rel_unc_key)
                if rel_unc_0 < convergence_crit[
                        'rel_unc_0'] and rel_unc_key < convergence_crit[
                            'rel_unc_key']:
                    print('breaking loop', counter)
                    print(rel_unc_0, convergence_crit['rel_unc_0'])
                    print(rel_unc_key, convergence_crit['rel_unc_key'],
                          converge_key)
                    break

        # Debugging output
        if debug and counter % 10000 == 0:
            # Compute the noise-to-signal ratio
            var = float(counter) * (Kcoeff_covar[0][0] -
                                    Kcoeff[0]**2) / float(counter - 1)
            ratio = np.sqrt(var / float(counter)) / Kcoeff[0]

            var = float(counter) * (Kcoeff_covar[-1][-1] -
                                    Kcoeff[-1]**2) / float(counter - 1)
            ratio_max = np.sqrt(var / float(counter)) / Kcoeff[-1]
            #print(counter, Kcoeff[0], ratio, ratio_max)
            print(counter, Kcoeff[0], Kcoeff[1])

    # Record the actual number of trials
    actual_trials = counter

    # Re-scale the running averages
    Kcoeff = [Kcoeff[j] / (scale_factor**j) for j in range(ncoeffs)]
    Kcoeff_covar = [[(Kcoeff_covar[j][l] /
                      (scale_factor**(j + l)) - Kcoeff[j] * Kcoeff[l]) *
                     float(trials) / float(trials - 1) for l in range(ncoeffs)]
                    for j in range(ncoeffs)]
    Kcoeff_var = [Kcoeff_covar[j][j] for j in range(ncoeffs)]

    if debug:
        print('done', trials)
        for i in range(ncoeffs):
            print(Kcoeff[i], np.sqrt(Kcoeff_var[i] / float(trials)))

    # Package the results in a dictionary
    results = copy.deepcopy(args)
    results['Kcoeff'] = Kcoeff
    results['Kcoeff_var'] = Kcoeff_var
    results['Kcoeff_covar'] = Kcoeff_covar
    results['trials'] = actual_trials  # store actual number of trials
    results['units'] = 'kJ/mol'

    return results
