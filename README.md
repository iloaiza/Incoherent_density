# Simulation of incoherent light excitations

Code written by Ignacio Loaiza G.
ignacio.loaiza@mail.utoronto.ca

Used for generating all the data and figures in
[1] - I. Loaiza, A. F. Izmaylov, and P. Brumer. Computational approaches to efficient generation of the stationary state for incoherent light excitation. (Submitted)
	

## Usage
Run using the Julia language. First include all libraries with 'include("include.jl")'.

If this is the first run and rhodopsin.h5 does not exist, generate it using fdebug(). Otherwise load using 'loading("rhodopsin")'. This loads rhodopsin model Hamiltonian and Franck-Condon wavefunction. fdebug.jl module also has examples for generating the LVC Hamiltonian and Franck-Condon wavefunction.


## General methodologies
Dynamical averaging: see methodologies in PLOT.jl

Lindbladian evolution: see methodologies in PLOT.jl and lindblad.jl

Lanczos shift-and-invert: see methodologies in shift.jl and parallel.jl


### General list of modules:

chebyshev.jl: Chebyshev interpolation procedure and Clenshaw rule for density matrices

config.jl: configuration constants

fdebug.jl: module for quick generation of Hamiltonians and wavefunctions

general_functions.jl: generic functions used throughout the code

heaviside.jl: generation of cis and trans projectors. See the ones used in PLOTS.jl for a working example with a given basis

include.jl: include module, includes all necessary packages and modules in correct order

krylov.jl: general functions involving krylov subspaces and Lanczos algorithm

lindblad.jl: module for Lindbladian decoherence

LVC_pot: generation of LVC system Hamiltonian

parallel.jl: parallel implementation of shift-and-invert Lanczos procedure for parallel analysis of random seeds

PARALLEL_RUN.jl: executable julia file, runs parallel.jl

PLOT.jl: functions used for generating figures on paper for dynamic averaging and Lindblad procedures

rhodopsin.jl: generation of rhodopsin system Hamitlonian

sanity.jl: some sanity check functions for debugging

saving.jl: saving and loading capabilities, uses HDF5 (to avoid generating everything anew on each simulation)

shift.jl: shift-and-invert Lanczos procedure implementation and plotting driver

tipt.jl: time-independent perturbation theory for corrected seed vector in Lanczos procedure