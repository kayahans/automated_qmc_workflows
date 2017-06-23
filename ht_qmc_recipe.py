#! /usr/bin/env python

##################################################################
##  (c) Copyright 2017-  by Kayahan Saritas                     ##
##################################################################

from nexus import generate_physical_system
from nexus import generate_pwscf,run_project
from nexus import generate_pw2qmcpack
from nexus import generate_qmcpack
from nexus import loop,linear,vmc,dmc

import sys
import os
import numpy as np
from datetime import datetime

# Global variables
mag_percent             = 0.7
dmc_samples_per_valence = 80000
dmc_steps_per_block     = 1
vmc_stats3              = 15000
vmc_stats               = 5000
date = "{:%B_%d_%H_%M}".format(datetime.now())


def pause():
    programPause = raw_input("Press the <ENTER> key to continue...")

def error_handler(string):
    print_scr(str(string) + " \n")
    print_scr("=====================================================")
    sys.exit()

def clean_gen_info(gen_info):
    #Update here if there other key items
    del gen_info["valency"]
    return gen_info

def check_scell_list(scell_list):
    if any(not isinstance(x, int) for x in scell_list):
        error_handler("Error: Non-integer values in scell_list")
    scell_list = sorted(scell_list)
    if any(x < 1  for x in scell_list):
        error_handler("Error: Min supercell size in scell_list is smaller than 1")
    scell_set = list(set(scell_list))
    if scell_list != scell_set:
        error_handler("Error: Non-unique elements in scell_list")
    return scell_list

def generate_structure(structure,**kwargs):
    if 'kgrid' in kwargs:
        structure["kgrid"] = kwargs["kgrid"]
    if 'kshift' in kwargs:
        structure["kshift"] = kwargs["kshift"]
    if 'scell' in kwargs:
        structure["tiling"] = kwargs["scell"]
    #structure = ", ".join(["=".join([key, str(val)]) for key, val in structure.items()])
    return  generate_physical_system(**structure)
# end def generate_structure

def get_supercells(axes,scell_list):
    # Find the most spherical diagonal supercell
    result = dict()

    a = np.matrix(axes)
    b_keep = []

    for scellsize in scell_list:
        max_lattice = 100000
        sa_range = (range(-scellsize, scellsize+1, 1))
        sa_range = [ x for x in sa_range if x != 0]
        for sa in sa_range:
            sb_range = (range(-abs(scellsize / sa), abs(scellsize/sa) + 1))
            sb_range = [x for x in sb_range if x != 0]
            if (isinstance(item, int) for item in sb_range):
                for sb in sb_range:
                    sc_range = (range(-abs(scellsize / (sa*sb)), abs(scellsize / (sa*sb)) + 1))
                    sc_range = [x for x in sc_range if x != 0]
                    if (isinstance(item, int) for item in sc_range):
                        for sc in sc_range:
                            if abs(sa*sb*sc) == scellsize:
                                b = np.array([sa, sb, sc])

                                b_tmp = np.multiply(a, b)
                                b_tmp = np.multiply(b_tmp, b_tmp)
                                b_tmp = np.sum(b_tmp, axis=1)
                                b_tmp = np.sqrt(b_tmp)
                                max_b = np.max(b_tmp)

                                if max_b <= max_lattice:
                                    max_lattice = max_b
                                    b_keep = b

        result.update({abs(sa*sb*sc):b_keep.tolist()})
    return result
# end def find_supercell

def transition_metal(atoms):
    # Find if there are any transition metals in the atoms list, assign 0.7 starting mag. to it in PWSCF
    t_metals = set([ 'Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn',\
    'Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',\
    'La','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg',\
    'Ac','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn'])
    if list(t_metals & set(atoms)):
        print_scr('= Magnetization is applied on : {\'' + list(t_metals & set(atoms))[0] + "' : " + str(mag_percent) + "}")
        return {list(t_metals & set(atoms))[0] : mag_percent}
    else:
        print_scr('= Magnetization is applied on : {\'' + list(atoms)[0] + "' : " + str(mag_percent) + "}")
        return {list(atoms)[0] : mag_percent}
# end def transition_metal

def num_val_e(structure):
    #print structure
    gen_info    = structure["generation_info"]
    particles   = structure["particles"]
    atoms = gen_info["atoms"]

    tot_valence = 0

    for atom in atoms:
        atom_valence = particles[atom]["charge"]-particles[atom]["core_electrons"]
        tot_valence = atom_valence + tot_valence

    return tot_valence
# end num_val_e

def dmc_steps(tot_val,dmc_job):

    tot_dmc_samples = tot_val * dmc_samples_per_valence
    dmc_processes = dmc_job["processes"]
    steps = tot_dmc_samples / dmc_processes

    return (steps, dmc_processes, tot_dmc_samples)
# end dmc_steps

def find_ecut(ecut,dft_pps):

    max_ecut = 0

    #If a file is provided as the dictionary of Ecut for pps: 1st row pp, 2nd row Ecut (Ry)
    if type(ecut) is str:

        if not os.path.exists(str):
            error_handler("No such directory exists: " + ecut)

        pp_file = ecut
        pp_lib = {}

        print_scr("Ecut is read from " + pp_file)
        with open(pp_file) as f:
            for line in f:
                (key, val) = line.split()
                pp_lib[key] = val


        for pp in dft_pps:
            if pp in pp_lib:
                if pp_lib[pp] > max_ecut:
                    max_ecut = pp_lib[pp]
            else:
                error_handler(pp + " is not found in " + pp_file)

    # Simple integer is provided
    elif type(ecut) is int:
        max_ecut = ecut

    #
    elif type(ecut) is dict:
        for pp, pp_ecut in ecut.iteritems():
            if type(pp) is str and type(pp_ecut) is int:
                if pp_ecut > max_ecut:
                    max_ecut = pp_ecut

    else:
        error_handler("Unknown type in ecut parameter")

    return max_ecut
# end find_ecut

def print_scr(str):
    f=open("ht_qmc."+date+".log", "a")
    print str
    f.write(str+ "\n")
    f.close()


def ht_qmc(structure, directory, dft_job, p2q_job, vmc_opt_job, dmc_job, dft_pps, qmc_pps, dft_functional, **kwargs):

    defaults = {
        'scell_list': [1],
        'ecut': 200,
        'vmc_steps': 1,
        'vmc_dt': 0.1,
        'phonon': False
    }

    for var, val in kwargs.iteritems():
        kwargs[var] = val

    undef_params = []
    for key, value in defaults.iteritems():
        if key in kwargs:
            defaults[key]=kwargs[key]
        else:
            undef_params.append(key)

    scell_list, ecut, vmc_steps, vmc_dt, phonon = [defaults.get(k) for k in ['scell_list', 'ecut', 'vmc_steps', 'vmc_dt', 'phonon']]


    print_scr("=====================================================")
    print_scr("= HT-QMC script")
    print_scr('= J. Chem. Theory Comput., 2017, 13 (5), pp 1943-1951 ')
    print_scr("=====================================================")

    if undef_params:
        print_scr("= Using default parameters for " + str(undef_params))
        print_scr("=====================================================")
    # gen_info is the dict from the generation info of the physical system in the input script
    gen_info = dict()

    for key, value in structure["generation_info"].iteritems():
        gen_info.update({key: value})
    gen_info = clean_gen_info(gen_info)

    # Structure of the unit cell
    structure = generate_structure(gen_info)

    # get the supercell matrices from integer list and create supercell matrices
    scell_list = check_scell_list(scell_list)
    scell_matrices = get_supercells(structure["structure"]["axes"], scell_list)

    # Total number of valence electrons in the input cell
    tot_val = num_val_e(structure)

    # Create input parameters, assuming 1 walker per process

    (tot_dmc_steps, dmc_walkers, tot_dmc_samples) = dmc_steps(tot_val, dmc_job)

    dmc_blocks = tot_dmc_steps / dmc_steps_per_block

    ecut = find_ecut(ecut,dft_pps)

    if phonon is True:
        print_scr("= Phonon calculations are not implemented yet")
    print_scr("= Using supercells of size: " + ' '.join(map(str, scell_list)) + " !               ")
    print_scr("=====================================================")
    # starting magnetization of a transition metal atom or a random atom in the init cell
    print_scr("=                DFT Settings: ")
    print_scr("= Using " + dft_functional + " DFT functional with " + str(ecut) + " Ry cutoff!    ")
    start_mag = transition_metal(gen_info["atoms"])
    print_scr("=====================================================")
    print_scr("=                VMC Settings: ")
    print_scr("= Using " + str(vmc_steps) + " VMC steps and " + str(vmc_dt) + " VMC timestep")
    print_scr("=====================================================")
    print_scr("=                DMC Settings: ")
    print_scr("= Using " + str(tot_dmc_steps) + " steps for DMC statistics!")
    print_scr("= Using " + str(dmc_walkers) + " processses and total of " + str(tot_dmc_samples) + " DMC samples!")
    print_scr("= (80000 samples per valence electron )")
    print_scr("=====================================================")


    sims = []
    scf_conv = generate_pwscf(
        # nexus inputs
        identifier='scf_conv',
        path=directory + '/scf_conv',
        job=dft_job,
        pseudos=dft_pps,
        system=structure,
        input_type='scf',
        input_dft=dft_functional,
        start_mag=start_mag,
        ecut=ecut,
        conv_thr=1e-8,
        mixing_beta=0.7,
    )
    sims.append(scf_conv)

    ks_min = 1000
    scell_min = []
    opt_prev = []

    for ks, scell in scell_matrices.iteritems():
        if ks_min > ks:
            ks_min = ks
            scell_min = scell

    for ks, scell in scell_matrices.iteritems():

        real_cell = dict(kgrid=(2, 2, 2), kshift=(1, 1, 1), scell=scell)
        generic_real = generate_structure(gen_info, **real_cell)
        opt_cell = dict(kgrid=(1, 1, 1), kshift=(0, 0.5, 0), scell=scell)
        generic_opt = generate_structure(gen_info, **opt_cell)

        scf_opt = generate_pwscf(
            identifier='scf_opt',
            path=directory + '/radius.' + str(ks) + '/scf_opt',
            job=dft_job,
            pseudos=dft_pps,
            system=generic_opt,
            input_type='scf',
            input_dft=dft_functional,
            nosym=True,
            spin_polarized=True,
            ecut=ecut,
            conv_thr=1e-7,
            mixing_beta=0.7,
            dependencies=(scf_conv, 'charge_density')
        )

        sims.append(scf_opt)

        p2q_opt = generate_pw2qmcpack(
            identifier='p2q',
            path=directory + '/radius.' + str(ks) + '/scf_opt',
            job=p2q_job,
            write_psir=False,
            dependencies=(scf_opt, 'orbitals')
        )
        sims.append(p2q_opt)

        scf_qmc = generate_pwscf(
            identifier='scf_qmc',
            path=directory + '/radius.' + str(ks) + '/scf_qmc',
            job=dft_job,
            pseudos=dft_pps,
            system=generic_real,
            input_type='scf',
            input_dft=dft_functional,
            nosym=True,
            spin_polarized=True,
            ecut=ecut,
            conv_thr=1e-7,
            mixing_beta=0.7,
            dependencies=(scf_conv, 'charge_density')
        )
        sims.append(scf_qmc)

        p2q_qmc = generate_pw2qmcpack(
            identifier='p2q',
            path=directory + '/radius.' + str(ks) + '/scf_qmc',
            job=p2q_job,
            write_psir=False,
            dependencies=(scf_qmc, 'orbitals')
        )

        sims.append(p2q_qmc)

        #VMC calculations

        vmc_blocks_1 = vmc_stats / vmc_steps
        vmc_blocks_2 = vmc_stats3 / vmc_steps
        linopt = linear(
            energy              = 0.0,
            unreweightedvariance= 1.0,
            reweightedvariance  = 0.0,
            timestep            = vmc_dt,
            warmupsteps         = 5000,
            walkers             = 1,
            steps               = vmc_steps,
            blocks 	        = vmc_blocks_1,
            substeps            = 12,
            samples 		= vmc_stats,
	    nonlocalpp          = True,
            usebuffer           = True,
            minwalkers          = 0.5,
            usedrift            = True,
            minmethod           = 'quartic',
            gpu          	    = 'no'
        )
        linopt2 = linopt.copy()
        linopt2.blocks = vmc_blocks_2
	linopt2.samples = vmc_stats3

        if ks == ks_min:
            opt = generate_qmcpack(
                identifier='opt',
                path=directory + '/radius.' + str(ks) + '/opt',
                job=vmc_opt_job,
                input_type='basic',
                system=generic_opt,
                bconds='ppp',
                pseudos=qmc_pps,
                jastrows=[('J1', 'bspline', 8), ('J2', 'bspline', 8)],
                calculations=[loop(max=4, qmc=linopt), loop(max=4, qmc=linopt2)],
                dependencies=[(p2q_opt, 'orbitals')]
            )
            opt_prev = opt
            sims.append(opt)
        else:
            opt = generate_qmcpack(
                identifier='opt',
                path=directory + '/radius.' + str(ks) + '/opt',
                job=vmc_opt_job,
                input_type='basic',
                system=generic_opt,
                bconds='ppp',
                pseudos=qmc_pps,
                jastrows=[('J1', 'bspline', 8), ('J2', 'bspline', 8)],
                calculations=[loop(max=4, qmc=linopt), loop(max=4, qmc=linopt2)],
                dependencies=[(p2q_opt, 'orbitals'),(opt_prev, 'jastrow')]
            )
            sims.append(opt)

        # DMC calculations



        vmc_dmc = vmc(
            timestep=vmc_dt,
            warmupsteps=100,
            blocks=vmc_blocks_2,
            gpu='no',
            steps=vmc_steps,
            substeps=12,
            samplesperthread=1
        )
        dmc_dmc = dmc(
            gpu='no',
            timestep=0.01,
            warmupsteps=500,
            blocks=dmc_blocks,
            steps=dmc_steps_per_block,
            nonlocalmoves=True
        )

        qmc = generate_qmcpack(
            identifier='qmc',  # identifier/file prefix
            path=directory + '/radius.' + str(ks) + '/qmc',  # directory for dmc run
            job=dmc_job,
            pseudos=qmc_pps,
            system=generic_real,
            input_type='basic',
            jastrows=[],
            calculations=[vmc_dmc, dmc_dmc],
            dependencies=[(p2q_qmc, 'orbitals'), (opt, 'jastrow')]
        )

        sims.append(qmc)

    run_project(sims)

