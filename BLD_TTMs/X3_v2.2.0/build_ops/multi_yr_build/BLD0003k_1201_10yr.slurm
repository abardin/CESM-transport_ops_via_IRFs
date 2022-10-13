#!/bin/tcsh
#
# GreenPlanet side 2 batch script to build multi-year operators
# 
#SBATCH --job-name=BLD0003k
#SBATCH -p sib2.9,nes2.8,brd2.4,has2.5
#SBATCH -N 1
#SBATCH --mem=20G
#SBATCH --ntasks-per-node=8
#SBATCH -t 24:00:00
#SBATCH -o BLD0003k.out.%J
#SBATCH -e BLD0003k.out.%J
#SBATCH --mail-type=END 
#SBATCH --mail-user=abardin@uci.edu

setenv TMPDIR /DFS-L/SCRATCH/moore/abardin/temp
mkdir –p $TMPDIR
module purge
module load matlab
module list
matlab -nodesktop -nodisplay -nosplash <<EOF

  % transport tracer ops
  F.ops_loc =  '/DFS-L/SCRATCH/moore/abardin/X3v2.2.0_ops/yr1201_3k/';
  F.noc_dir = ([F.ops_loc,'noc/']);  %monthly unadjusted ops
  if ~exist(F.noc_dir,'dir')
    disp(['no noc dir']);
    keyboard
  end

  F.nyears      = 10;
  F.zeroyr      = 1200;
  F.calyr_zero  = 1200;  
  F.firstyr     = F.zeroyr + 1;
  F.firstcalyr  = F.calyr_zero + 1; 
  F.ops_dir  = ([F.ops_loc,'X3_',int2str(F.nyears),'yr_',int2str(F.firstcalyr),'k/']); 
  if ~ exist(F.ops_dir,'dir')
     eval(['!mkdir ', F.ops_dir]);
  end
  F

  Build_multi_period(F); 
  Build_mo_climate_state(F);                  
  
  exit
EOF
exit 

