#!/bin/tcsh
#
# GreenPlanet batch script to build operators from IRF output
# 
#SBATCH --job-name=BLDops
#SBATCH -p sib2.9,nes2.8,brd2.4,has2.5
#SBATCH -N 1
#SBATCH --mem=20G
#SBATCH --ntasks-per-node=16
#SBATCH -t 24:00:00
#SBATCH -o phyN2150.out.%J
#SBATCH -e phyN2150.out.%J
#SBATCH --mail-type=END 
#SBATCH --mail-user=abardin@uci.edu

setenv TMPDIR /DFS-L/SCRATCH/moore/abardin/temp
module purge
module load matlab
module list
matlab -nodesktop -nodisplay -nosplash <<EOF
  % output from IRF run
  IRFcasename = 'phyIRF2150';    % name in history archive
  short_casename = 'phyN2150';   % 8 char
  F.nyears   = 10;               % number of years in average
  F.zeroyr   = 2149;             % year before first archived year
  F.firstyr  = F.zeroyr + 1;     % first year to begin average 
  F.firstcalyr = F.firstyr;      % used for name of transport ops 
  F.day5     = false;            % monthly if day5 = false
  F.num_periods = 12;            % number of periods in a year (12 or 73)
  F.IRFcasename = IRFcasename;
  % directory for in-process build files (temporary)
  F.tmp_ops_dir = (['/DFS-L/SCRATCH/moore/abardin/temp/']);

  F.ref_dir  = (['/DFS-L/SCRATCH/moore/abardin/archive/', ...
                    IRFcasename,'/ocn/hist/']);

  F.ref_file_base = ([IRFcasename,'.pop.h.']); 
  F.ref_file = sprintf('%s%04i-01.nc',F.ref_file_base,F.firstyr);

  % IRF tracer definition
  F.IRFdir  = '~/X3_v2.2.0/CESM_run/IRF_mods/';
  F.IRFfile = 'IRF_offline_transport_tracers_gx3v7_20150314.nc'; 

  % IRF stencils (aka masks) 
  F.irf_stencil_dir = '/DFS-L/DATA/moore/abardin/IRF_masks/CESM2.2.0_X3/';

  % transport tracer ops
  F.ops_loc = sprintf('%s%04i_%s/', ...
              '/DFS-L/SCRATCH/moore/abardin/X3v2.2.0_ops/yr', ...
              F.firstyr,short_casename);
  if ~exist(F.ops_loc,'dir')
     eval(['!mkdir ',F.ops_loc]);
  end
  F.noc_dir = ([F.ops_loc,'noc/']);  % unadjusted ops
  if ~exist(F.noc_dir,'dir')
     eval(['!mkdir ',F.noc_dir]);
  end

 F.ops_dir  = sprintf('%sX3_%02iyr_%04i/',F.ops_loc,F.nyears,F.firstcalyr); 
  if ~ exist(F.ops_dir,'dir')
     eval(['!mkdir ', F.ops_dir]);
  end

  F  

  cd noc_ops_build                   
  
  MET = BuildMET(F);
  eval(['save ',F.noc_dir,'MET.mat MET -v7.3']);
  for year = F.firstyr:(F.firstyr+F.nyears-1)
     for period = 1:F.num_periods
        % Build_state_data (F, year, period);
        Build_monthly_ops(F, year, period);   
     end % for each period
  end % for year
  
 cd ../multi_yr_build
 Build_multi_period(F); 

 cd ..

 Ffilename = (['Ffile_',short_casename,'.mat']);
 eval(['save ',Ffilename,' F ']);

 eval(['save runname.mat short_casename ']);
exit  
EOF

pushd ~/X3_v2.2.0/tracers
sbatch tag_deep_tracers.slurm
popd

exit 
