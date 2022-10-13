READMEfiles for build_ops directory

noc_ops_build is a directory containing the code for building the un-adjusted
monthly operators from the output of the parent model IRF run.

multi_yr_build is a directory containing the code to build the seasonal set
of operators, with the SSH (sea surface height) adjustment as a climatological
average for each month.

ops_exam is a directory containing some diagnostic analysis for the operators.
The one-consistency test is here.

The .slurm files are job procedures that shortcut through building the noc
operators to build the multi-year climatologically averaged set of operators.
There is an assumption that this is not a new version of POP, allowing reuse
of pre-built structures.

The runname.mat and Ffile's are "tag" files that define the build case, used to 
carry forward the build information between build phases, and to later analysis
using deep-water tracers.

