# 2026-IFF

This repository is a home to the code written for the International Flavors and Fragrances (IFF) project in the 2026 University of Delaware Spring Hackathon (CHEG867) class. This project has involved running molecular simulations of different malodor and fragrance molecules in order to predict how they interact with the surface of wet/damp hair (see https://sites.udel.edu/midas-nrt/graduate-trainee-timeline/courses-for-nrt-trainees/nrt-hackathon-course/ for additional information and final report). 
This repository contains the codes used to set up and analyze the simulations, as well as the resulting data collected from the simulations.

The repository is organized into the following groups of code, and each directory contains the additional required information needed to use the code:
+ **inputs**: directory containing necessary forcefield files, simulation mdp runfiles, and odorant-specific topology files
+ **gen_sys**, **ff_data**, **gromacs_helpers**: Python scripts used to write the data files (including initial bead locations, bonds, and angles, bonded and nonbonded interaction coefficients, and bead masses) that are read by GROMACS for the simulations completed in this project
+ **hair_model**: model of the base system comprising substrate, lipid layers, and water vapor, contained ~69,800 CG beads
+ **DMS/LIN/PHE/SYR_single**: models for studying Aim 1: adsorption of a single odorant molecule
+ **DMS/LIN/PHE/SYR_droplet**: models for studying Aim 2: adsorption of odorant clusters
+ **results**: directory containing energy & temperature plots, RDFs, MSDs as produced by GROMACS software for each odorant and system type
