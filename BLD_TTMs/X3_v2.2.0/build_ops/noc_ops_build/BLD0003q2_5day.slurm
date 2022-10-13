#!/bin/tcsh
#
# GreenPlanet side 2 batch script to build noc operators 
# 
#SBATCH --job-name=BLD0003q2
#SBATCH -p sib2.9,nes2.8,brd2.4,has2.5
#SBATCH -N 1
#SBATCH --mem=20G
#SBATCH --ntasks-per-node=16
#SBATCH -t 24:00:00
#SBATCH -o BLD0003q2.out.%J
#SBATCH -e BLD0003q2.out.%J
#SBATCH --mail-type=END 
#SBATCH --mail-user=abardin@uci.edu

setenv TMPDIR /DFS-L/SCRATCH/moore/abardin/temp
mkdir –p $TMPDIR
module purge
module load matlab
module list
matlab -nodesktop -nodisplay -nosplash <<EOF
  % output from IRF run
  F.ref_dir  = (['/DFS-L/SCRATCH/moore/abardin/archive/', ...
                    'hdev0003q2/ocn/hist/']);
  F.ref_file =      'hdev0003q2.pop.h.4051-01-05.nc';
  F.ref_file_base = 'hdev0003q2.pop.h.';

  % IRF tracer definition
  F.IRFdir  = '~/X3_v2.2.0/CESM_run/IRF_mods/';
  F.IRFfile = 'IRF_offline_transport_tracers_gx3v7_20150314.nc'; 

  % IRF stencils (aka masks) 
  F.irf_stencil_dir = '/DFS-L/DATA/moore/abardin/IRF_masks/CESM2.2.0_X3/';

  % transport tracer ops
  F.ops_loc =  '/DFS-L/SCRATCH/moore/abardin/X3v2.2.0_ops/yr4051_5day_3q2/';
  if ~exist(F.ops_loc,'dir')
     eval(['!mkdir ',F.ops_loc]);
  end
  F.noc_dir = ([F.ops_loc,'noc/']);  % unadjusted ops
  if ~exist(F.noc_dir,'dir')
     eval(['!mkdir ',F.noc_dir]);
  end

  F.nyears   = 10;
  F.zeroyr   = 4050;  
  F.firstyr  = F.zeroyr + 1;
  F.day5     = true
  F.num_periods = 73;

  % directory for in-process build files (temporary)
  F.tmp_ops_dir = (['/DFS-L/SCRATCH/moore/abardin/temp/']);
  F                     
  
  MET = BuildMET(F);
  eval(['save ',F.noc_dir,'MET.mat MET -v7.3']);
  for year = F.firstyr:(F.firstyr+F.nyears-1)
     for period = 1:F.num_periods
        Build_state_data (F, year, period);
        Build_monthly_ops(F, year, period);   % contains additional F.parameter
                                              % definitions specific to v. 2.2.0
     end % for each period
  end % for year
  
exit
EOF
exit 

