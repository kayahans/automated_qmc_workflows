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
import numpy as np

def pause():
    programPause = raw_input("Press the <ENTER> key to continue...")

def error_handler(string):
    print str(string) + " \n"
    print "====================================================="
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
        print '= Magnetization is applied on : {\'' + list(t_metals & set(atoms))[0] + "' : 0.7}"
        return {list(t_metals & set(atoms))[0] : 0.7}
    else:
        print '= Magnetization is applied on : {\'' + list(atoms)[0] + "' : 0.7}         "
        return {list(atoms)[0] : 0.7}
# end def transition_metal


def ht_qmc(structure, directory, dft_job, p2q_job, vmc_opt_job, dmc_job, dft_pps, qmc_pps, dft_functional, scell_list, ecut, vmc_steps, vmc_dt,**kwargs):
    for var, val in kwargs.iteritems():
        kwargs[var] = val
    # end for
    print "====================================================="
    print "= HT-QMC script"
    print '= J. Chem. Theory Comput., 2017, 13 (5), pp 1943-1951 '
    print "====================================================="

    # gen_info is the dict from the generation info of the physical system in the input script
    gen_info = dict()

    for key, value in structure["generation_info"].iteritems():
        gen_info.update({key: value})
    gen_info = clean_gen_info(gen_info)

    scell_list = check_scell_list(scell_list)
    scell_matrices = get_supercells(structure["structure"]["axes"], scell_list)

    print "= Using supercells of size: " + ' '.join(map(str, scell_list)) + " !               "
    print "====================================================="
    # starting magnetization of a transition metal atom or a random atom in the init cell
    print "=                DFT Settings: "
    print "= Using " + dft_functional + " DFT functional with " + str(ecut) + " Ry cutoff!    "
    start_mag = transition_metal(gen_info["atoms"])
    print "====================================================="
    print "=                VMC Settings: "
    print "= Using " + str(vmc_steps) + " and " + str(vmc_dt)
    print "====================================================="
    structure = generate_structure(gen_info)
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

        vmc_blocks_1 = 50000 / vmc_steps
        vmc_blocks_2 = 150000 / vmc_steps

        linopt = linear(
            energy              = 0.0,
            unreweightedvariance= 1.0,
            reweightedvariance  = 0.0,
            timestep            = vmc_dt,
            warmupsteps         = 5000,
            walkers             = 256,
            steps               = vmc_steps,
            blocks 	            = vmc_blocks_1,
            substeps            = 12,
            nonlocalpp          = True,
            usebuffer           = True,
            minwalkers          = 128,
            usedrift            = True,
            minmethod           = 'quartic',
            gpu          	    = 'yes'
        )
        linopt2 = linopt.copy()
        linopt2.blocks = vmc_blocks_2
        linopt2.warmupsteps = 10000

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
            warmupsteps=10000,  
            blocks=vmc_blocks_2,
            gpu='no',
            steps=vmc_steps,
            substeps=3,
            samplesperthread=128
        )
        dmc_dmc = dmc(
            gpu='no',
            timestep=0.01,
            warmupsteps=500,
            blocks=500,
            steps=5,
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
