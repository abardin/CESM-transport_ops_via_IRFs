function Build_multi_period(F)
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
%  buildMET, which gathers the geometric data from a netcdf history file.
%
%  make_A_periodic
%    -- computes the annual periodic adjustment for the surface layer,
%       and applies it to all the A and T operators.
%    -- computes the monthly SSH-change adjustment for the surface layer.
%    -- saves files containing the operators, organized by month-of-year, 
%       in the final directory.
%   
%   build_model_state_data, which builds a set of monthly files
%       containing SSH, TEMP (potential temprature), PD (potential density)
%       SALT (salinity), RHO (in-situ density), IFRAC (sea-ice fraction).
%       Only SSH is used for the transport; the other variables are
%       frequently needed for biogeochemical models.
%
%  ** Note that the IRF-mask files are required to be pre-built to execute
%  this program, which builds the monthly operators. The IRF-mask files
%  are listed as part of the input files listed below.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define variables that change with the configuration
%   % read year from file
%  fID = fopen('year_no.txt','r');
%  formatSpec = '%f';
%  myyr = fscanf(fID,formatSpec);
%  fclose(fID);

%% Control inputs:
  % F.nyears
  % F.firstyr
    % F.zeroyr
      % F.noc_dir
        % F.ops_dir
  
%  myfn = (['Ayear',int2str(myyr)]);  % acknowlege year read
%  fID = fopen(myfn, 'w');
%  fprintf(fID, formatSpec, myyr);
%  fclose(fID);
      
   
  % expected dimensions in the POP2 history files:
  OGCM_nx = 100; OGCM_ny = 116; OGCM_nz = 60;
                                       
%% run parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if F.day5 == true
      P.monthly = false;
  else
      P.monthly = true;         % whether files are monthly or daily
  end
  P.nperiods =  F.num_periods*F.nyears;     % parameter for number of days % months to get data 
                            % for
  P.fperiod =  F.num_periods*(F.firstyr - F.zeroyr - 1) + 1;    % parameter for first month to get data for

%  P.fday = 1;                % first day of run to process for daily IRF output
%  P.day_span = 1;            % span between days for daily IRF output
%  P.lday = 10;               % last day of run to process for daily IRF output
                                                   

                          
  P.max_exp_err = 9e-11;     % define large vs. expected error level for
                             % operator mass-balance and divergence checks.
                             % See POP Reference Manual (2010), Section 4.5.3,
                             % in which an error tolerance in the convergence
                             % criteria is estimated at between 10^-12 and 
                             % 10^-13. The average monthly operators 
                             % contain row divergence as large as 5e-11
                             % below the surface layer.
                             
 % P.first_month_in_annual = P.fmo;   % first month included  in the annual average
                                 % included months = first +
                                 % months_per_year -1.
                                
 % P.months_per_annual = 12*F.nyears;   % number of months in a multi-year avg
   
  P.fignum = 10;
 
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

%% Control flow  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
  % Build the geometric data that is standard
  % MET = buildMET(F.ref_dir, F.ref_file);

  % Read the geometric data file from another ops directory
  MET = get_MET(F.noc_dir);
  
  % Check that the dimensions are the expected ones  
  [ny,nx,nzz] = size(MET.MASK);
  if ((nx == OGCM_nx) && (ny == OGCM_ny) && (nzz == OGCM_nz+1) )
      % continue
  else % unexpected dimensions
      disp ( 'check unexpected grid dimensions');
      keyboard
  end % if dimensions as expected
  %Save in ops dir being set up
  if ~exist(F.ops_dir,'dir')
      eval(['!mkdir ',F.ops_dir])
  end
  eval(['save ',F.ops_dir,'MET.mat MET ']);
  disp (['saved  MET in ',F.ops_dir]);
 
  % make_annual_ops
  % makes the annual average operator from the periodic  operators
  make_annual_ops(MET,F,P); 
  
  % make_A_periodic
  %   Computes A2 periodic adjustment for the surface layer. 
  %   Applies it to A and T operators for months in the annual set
  %   Output to files, organized by calendar month in the final directory,  
  %   for months in the annual set, plus the annual average.
  %   The files contain A, H, D, T, dxidt.
  %   Makes mass_balance and row divergence checks.
  make_A_periodic(F,P,MET);
  
  disp(['completed building annual and periodic ops in ',F.ops_dir])
% end of BUILD_OPS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  30 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
