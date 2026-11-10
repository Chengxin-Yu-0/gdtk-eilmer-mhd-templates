# gdtk-eilmer-mhd-templates

-- Chengxin Yu, updated 10 November 2026

This repository contains some template user-defined-functions (source terms and supervisor) for users running any magnetohydrodynamics (MHD) simulations using Eilmer5/lmr in the Gas Dynamics Toolkit (GDTK).

At the start of each step, the supervisor file reads a magnetic field from an external file and interpolates it to Eilmer's fluid grid, assigning each cell some B (B.x, B.y, and B.z if the simulation is 3D) values. The source term file then calculates the Lorentz force in each cell at each step, adding it to the momentum equations. The syntax for the supervisor file (steady vs transient, MPI vs non-MPI) is rather different, therefore a template for each is included. A source terms file with Raiser's interpolated conductivity model is also included. <<This repository is actively being updated so there will be more source terms templates to come, especially a template of the simplest case, which is a spherical magnet that can be modelled with the dipole equation.>>

Currently these templates are experimenter-focussed. Note that this is not meant to be a comprehensive catalog and users very likely still need to make some edits to these files depending on what they are trying to simulate, but hopefully this repository provides a starting point and allows any changes to be documented and communicated.

Most of the work developing these codes were done by previous PhD students (Aaron Cunningham, Sebastiaan van Oeveren, Alexis Lefevre, Daniel Smith) and Peter Jacobs.