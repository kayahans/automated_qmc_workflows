#! /usr/bin/env python

from ht_qmc_recipe import ht_qmc
from nexus import settings,Job,run_project
from nexus import generate_physical_system

# general settings for nexus
settings(
    pseudo_dir    = '../pseudopotentials',
    runs          = '',
    results       = '',
    status_only   = 0,
    generate_only = 1,
    sleep         = 2,
    machine       = 'ws16',
    )

# locations of pwscf, pw2qmcpack and qmcpack executables
pwscf           = 'pw.x'
pw2qmcpack      = 'pw2qmcpack.x'
qmcpack		    = 'qmcpack'

generic = generate_physical_system(
    lattice   = 'hexagonal',      # hexagonal cell shape
    cell      = 'primitive',      # primitive cell
    centering = 'P',              # primitive basis centering
    constants = (2.462,5.0),     # a,c constants
    units     = 'A',              # in Angstrom
    atoms     = ('C','C'),        # C primitive atoms
    basis     = [[ 0  , 0  , 0],  # basis vectors
                 [2./3,1./3, 0]],
    tiling    = (1,1,1),          # tiling of primitive cell
    kgrid     = (8,8,1),          # Monkhorst-Pack grid
    kshift    = (.5,.5,.5),       # and shift
    C         = 4                 # C has 4 valence electrons
)

ht_qmc(
    structure       = generic,
    directory       = './graphene',
    dft_job         = Job(nodes=2, minutes=120, app=pwscf),
    p2q_job         = Job(cores=1, minutes=15, app=pw2qmcpack),
    vmc_opt_job     = Job(nodes=4, threads=16,minutes=120, threads=8, app=qmcpack),
    dmc_job         = Job(nodes=256, threads=16, minutes=120, threads=8, app=qmcpack),
    dft_pps         = ['C.BFD.upf'],
    qmc_pps	        = ['C.BFD.xml'],
    dft_functional  = 'lda',
    ecut            = 160,
    scell_list      = [1,   2],
    vmc_steps       = 200,
    vmc_dt          = 0.1
)
