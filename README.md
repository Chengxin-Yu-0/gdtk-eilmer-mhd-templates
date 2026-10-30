# gdtk-eilmer-mhd-templates

These are some template MHD user-defined function (udf) files that may be helpful for anyone trying to run an MHD simulation in eilmer gdtk. These files are more experimenter-focussed.

The supervisor file allow the user to load in a magnetic field from any source and ramp the field strength up incrementally (to reduce possibility numerical errors). The source term file then compute the Lorentz force using the magnetic field assigned by the supervisor file. Note that different functions must be called for steady vs transient solver. 

Majority of the work was done by Aaron Cunningham and Peter Jacobs
