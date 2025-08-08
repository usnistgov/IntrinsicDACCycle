# -*- coding: utf-8 -*-
"""Functions for setting the Ewald Parameters from Heuristics"""
import numpy as np
import feasst as fst
# pylint: disable-msg=invalid-name   #because I use snake_case


def DL_POLY_ewald(rcut, box, tolerance):
    """ Parameters from the DL_POLY Algorithm
    https://doi.org/10.1080/002689798167881
    """
    # DL_POLY Algorithm
    eps = min(tolerance, 0.5)
    xi = np.sqrt(np.abs(np.log(eps * rcut)))
    alpha = np.sqrt(np.abs(np.log(eps * rcut * xi))) / rcut
    chi = np.sqrt(-np.log(eps * rcut * ((2. * xi * alpha)**2)))
    kmax = [int(0.25 + box[i] * alpha * chi / np.pi) for i in range(3)]
    return alpha, kmax


def LAMMPS_ewald(system, tolerance):
    """ Parameters from the LAMMPS Algorithm
    https://doi.org/10.1080/002689798167881
    """
    # LAMMPS algorithm
    ewald = fst.MakeEwald(fst.args({
        'tolerance': str(tolerance),
    }))
    system.add(fst.MakePotential(ewald))
    system.precompute()
    alpha = system.configuration().model_params().property('alpha')
    kmax = [
        ewald.num_kx() - 2,
        int((ewald.num_ky() - 3) / 2),
        int((ewald.num_kz() - 3) / 2)
    ]
    return alpha, kmax


def MCCCS_ewald(rcut, box, mcccs_type):
    """Parameters from either MCCCS-Towhee or MCCCS-MN"""
    if mcccs_type == 'TOWHEE':
        alpha = 5.60 / min(box)
        kmax = [10] * 3
    elif mcccs_type == 'MN':
        alpha = 6.4 / min(box)
        kmax = [int((alpha**2) * x * rcut / np.pi) + 1 for x in box]
    if rcut < min(box) / 2.:
        print('WARNING: Ewald Parameters are based on MCCCS heuristic that\n',
              'is invalid for current system')
        print(' MCCS Heuristic is designed for rcut = min(box)/2')
        print('   rcut:     ', rcut)
        print('   min(box): ', min(box))
    return alpha, kmax


def Cassandra_ewald(rcut, box, tolerance):
    """
    Parameters from Cassandra
    https://doi.org/10.1080/08927029408022180

    DWS NOTE 2021-05-24: This recipe generates an \alpha
    that is much larger than others; kmax is much larger
    as well. I have a suspicion that the logarithm in the
    calculation of \alpha should be base-10, not base-e.
    """

    alpha = np.sqrt(-np.log(tolerance)) / rcut
    h_ewald_cut = 2. * (-np.log(tolerance)) / rcut
    kmax = [int(h_ewald_cut * x / 2. / np.pi) + 1 for x in box]
    return alpha, kmax


def set_ewald(ewald_type, system_obj, tolerance=1.e-5, minmax='min'):
    """
    Wrapper function to capture input from main program and farm work to subsidiary functions
    """
    # Identify the box dimensions
    box = [
        system_obj.configuration().domain().side_length(i) for i in range(3)
    ]

    # Set the rcut for further operations, from particle files
    if minmax == 'min':
        rcut = min([
            system_obj.configuration().model_params().select('cutoff').value(i)
            for i in range(system_obj.configuration().model_params().select(
                'cutoff').size())
        ])
    elif minmax == 'max':
        rcut = max([
            system_obj.configuration().model_params().select('cutoff').value(i)
            for i in range(system_obj.configuration().model_params().select(
                'cutoff').size())
        ])
    else:
        raise Exception('minmax must be min or max')

    if ewald_type in ['DL_POLY', 'RASPA']:
        # DL_POLY Algorithm
        alpha, kmax = DL_POLY_ewald(rcut, box, tolerance)
    elif ewald_type == 'LAMMPS':
        # LAMMPS Algorithm
        alpha, kmax = LAMMPS_ewald(system_obj, tolerance)
    elif ewald_type == 'Cassandra':
        # Cassandra Algorithm
        alpha, kmax = Cassandra_ewald(rcut, box, tolerance)
    elif ewald_type in ['TOWHEE', 'MCCCS-TOWHEE']:
        # MCCCS-TOWHEE
        alpha, kmax = MCCCS_ewald(rcut, box, 'TOWHEE')
    elif ewald_type == 'MCCCS-MN':
        # MCCCS-MN
        alpha, kmax = MCCCS_ewald(rcut, box, 'MN')
    else:
        raise Exception('ERROR Unknown Ewald Type: ', ewald_type)

    return alpha, kmax
