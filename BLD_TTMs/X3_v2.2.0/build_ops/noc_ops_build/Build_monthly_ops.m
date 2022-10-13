function Build_monthly_ops(F,year,period)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The purpose of Build_ops is to build the advection/diffusion operators
% for the offline model from the history files resulting from running
% the POP ocean model with the "IRF tracers" turned on. This script 
% invokes functions which output files containing a sparce matrix with each
% operator(A (advection), H (horizontal diffusion, D (vertical diffusion), 
% and T (total transport matrix), for each month of the year.
%
% The first part of Build_ops contains the definitions of the files 
% and parameters that control the build process. 
% It invokes the following functions:
%
%  make_ops_from_IRF_output takes care of 
%    -- looping through the netcdf history files to collect the data;
%    -- invokes functions that get the specific data needed for
%       each operator;
%    -- combines the operators and saves them in files.
%
%  ** Note that the IRF-mask files are required to be pre-built to execute
%  this program, which builds the monthly operators. The IRF-mask files
%  are listed as part of the input files listed below.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define variables that change with the configuration

  % read year from file
%   fID = fopen('year_no.txt','r');
%   formatSpec = '%f';
%   myyr = fscanf(fID,formatSpec);
%   fclose(fID);

  % read the job id
%  fID = fopen(['JOB_ID_',int2str(period),'.txt'],'r');
%  formatSpec = '%f';
%  F.JOB_ID = fscanf(fID,formatSpec);

  F.file_yr = year;

  % expected dimensions in the POP2 history files:
  OGCM_nx = 100; OGCM_ny = 116; OGCM_nz = 60;
  
  % history file to use for time-unchanging geometry

  % directory for IRF_masks, calculated in build_IRF_masks.m
  % F.irf_stencil_dir = '/DFS-L/DATA/moore/abardin/IRF_masks/CESM2.1.3_X3/';
  F.g_irf_mask_file = 'g_irf_mask.mat';
  F.s_irf_mask_file = 's_irf_mask.mat';
  F.e_irf_mask_file = 'e_irf_mask.mat';

  % directory from which to get the history data from the IRF OGCM run
  F.g_dir  = F.ref_dir;
  F.g_file = F.ref_file_base;         % first part of file name
  
  F.s_dir   = F.ref_dir;
  F.s_file  = F.ref_file_base;        % first part of filename
  F.e_dir   = F.ref_dir;
  F.e_file  = F.ref_file_base;        % first part of filename
  F.VDC_dir = F.ref_dir;              % info for Vertical Diffusion is available from any 
                                      % POP source file
  F.VDCfile = F.ref_file_base;
%  F.IRFdir  % defined by input       % dir /file for IRF tracer definition
%  F.IRFfile     
 
  F.k_ntracers = 25*5;                 % number of tracers for g (standard)
  F.s_ntracers = 15;                   % number of tracers for s 
                                       % (ovfl source) runs
  F.e_ntracers = 16;                   % number of tracers for e 
                                       % (ovfl entrainment) runs
                                       
%% run parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  P.monthly = ~F.day5;         % whether files are monthly or daily
  P.file_yr = F.file_yr;
  P.nmo =  12;               % parameter for number of days % months to get data 
                             % for (normally = 12, 1 for BUILD_FAST)
  P.fmo =  12*(P.file_yr - F.zeroyr - 1) + 1;  % parameter for first month to get data for 
                             

  P.fday = 1;                % first day of run to process for daily IRF  output
  P.day_span = 5;            % span between days for daily IRF output
  P.lday = 365;               % last day of run to process for daily IRF output
                                                   
  P.fcase = 1;               % used in make_ops_from_IRF_output.m
                             % provided to allow for the remaking of
                             % only part of the operators while debugging.
                             % first case to run (IRF g,klev = 1)
                             % Note: cases 1-5 process the regular
                             % IRF output k-level 1-5, where k-level
                             % is a set of depth levels.  Cases 6 and 7
                             % process the source and entrainmet 
                             % over-the-sill flow specific regions.
                             % Case 8 processes the vertical diffusion
                             % parameters.
  P.lcase = 4;               % last case to run -- normally run 1:4
                             % to cover all the operators.
                             % (hint: it is sometimes useful to use fewer
                             % cases when debugging
                             
                             
      seq_num =  F.num_periods*(P.file_yr - F.zeroyr - 1) + period;                             
      if seq_num == 1                                                    
         P.collect_pulse = true;    % whether to collect pulse location data to check
      else                          % if one and only one pulse is set for each point
         P.collect_pulse = false;   % Note: this is turned off after the first 
      end                           % period by the software, since it
                                    % is the same for all months.
                          
      P.max_exp_err = 5e-11;     % define large vs. expected error level for
                             % operator mass-balance and divergence checks.
                             % See POP Reference Manual (2010), Section 4.5.3,
                             % in which an error tolerance in the convergence
                             % criteria is estimated at between 10^-12 and 
                             % 10^-13. The average monthly operators 
                             % contain row divergence as large as 3e-11
                             % below the surface layer
 
     P.days_per_mo = zeros(12,1);
     P.days_per_mo(1)  = 31;
     P.days_per_mo(2)  = 28;
     P.days_per_mo(3)  = 31;
     P.days_per_mo(4)  = 30;
     P.days_per_mo(5)  = 31;
     P.days_per_mo(6)  = 30;
     P.days_per_mo(7)  = 31;
     P.days_per_mo(8)  = 31;
     P.days_per_mo(9)  = 30;
     P.days_per_mo(10) = 31;
     P.days_per_mo(11) = 30;
     P.days_per_mo(12) = 31;
     P.days_per_yr     = sum(P.days_per_mo);
% F
% P
%% Control flow  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
  % load MET
  eval(['load ',F.noc_dir,'MET.mat MET']);
  
  % Check that the dimensions are the expected ones  
  [ny,nx,nzz] = size(MET.MASK);
  if ((nx == OGCM_nx) && (ny == OGCM_ny) && (nzz == OGCM_nz+1) )
      % continue
  else % unexpected dimensions
      disp ( 'check unexpected grid dimensions');
      keyboard
  end % if dimensions as expected
    
  % make_ops_from_IRF_output for given period
  %   Builds iocn-sized operators
  %   Does checks for the IRF pulses, and for operators not-on-land.
  % 
  %   Output to files 
  %       - numbering is per seq_num
  %         starting from the zeroyr parameter defined above
  %       - files are output to the F.noc_dir directory
  %

  make_ops_from_IRF_output(MET,F,P, year,period);
  
%  complete_fn = (['complete',int2str(period),'.txt']);
%  fID = fopen(complete_fn,'w');
%  nbytes = fprintf(fID,'%s','3');
% fclose(fID);
  disp(['completed building ops in ',F.noc_dir,' for year ',int2str(year),', period ',int2str(period)]);

return
end % function Build_monthly_ops

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  30 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
