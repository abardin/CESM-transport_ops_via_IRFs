#!/bin/csh
#===============================================================================
#  This is a CCSM batch job script to run the ones consistency test
#===============================================================================

#SBATCH --job-name=ones_ops
#SBATCH --exclusive
#SBATCH -p brd2.4,has2.5,ilg2.3,m-c1.9,m-c2.2,nes2.8,sib2.9
#SBATCH -N 1
#SBATCH --mem 20gb
#SBATCH --ntasks-per-node=16
#SBATCH -t 60:00:00
#SBATCH -o physci_new.out.%J
#SBATCH -e physci_new.out.%J
#SBATCH --mail-type=END 
#SBATCH --mail-user=abardin@uci.edu

module purge
module load matlab
module list
 
setenv TMPDIR /DFS-L/SCRATCH/moore/abardin/temp

matlab -nodesktop -nodisplay -nosplash <<EOF
 F.casename = 'X3phyNew';
 F.disp_case = 'X3 phyNew';
 F.nyears = 10;            % number of years to run test
 F.firstyr = 2090;
 F.timestep_hrs = 3;
 F.num_periods = 12;       % number of periods per year
 F.ops_dir = '/DFS-L/SCRATCH/moore/abardin/X3v2.2.0_ops/yr1201_3k/X3_10yr_1201k/'; 
 F.out_base = (['/DFS-L/SCRATCH/moore/abardin/',F.casename,'/']);
 F.out_dir = ([F.out_base,'opsExam_nyrs',int2str(F.nyears),'/']);
 if ~ exist(F.out_base,'dir')
     eval(['!mkdir ', F.out_base]);
 end
 if ~ exist(F.out_dir,'dir')
     eval(['!mkdir ', F.out_dir]);
 end

 F
 Run_ones_test(F);
 exit
EOF
exit 
