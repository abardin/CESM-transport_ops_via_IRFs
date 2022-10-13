function make_annual_ops(MET,F,P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The purpose of make_annual_ops is to build the annually-averaged
% operators from the monthly operators, which have already been built.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % display input parameters
  % MET                     % structure with geometric grid definitions
  F                       % structure with file location definitions
  P                       % structure with control parameters

% annual averages

  disp('making annual average operators');

  [A, H, D, T, num_days] = ...
         annual_avg_dmo(MET,F,P);
    
% save annual averages to file
   out_file_name = ([F.noc_dir,'MTM0.mat']);
   cmd = (['save ',out_file_name,' A H D T num_days -v7.3;']);  
   eval(cmd);
   disp([ out_file_name,' saved']);
         
  return
end % make_annual_ops

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A_an,H_an,D_an,T_an, days_per_year] = ...
         annual_avg_dmo(MET, F, P) 
%
% Annual_avg computes the annual average from the monthly averages
% of the iocn-sized matrix operators that are within 
% P.first_month_in_annual + P.months_per_year -1.
%
% INPUTS:
%   MET                      % geometric parameters for the grid
%   F.ops_dir                % directory containing the monthly operators, 
%                             
%   P                        % control paramters
%    first_month_in_annual   % first month to be included in the annual
%                            % average - provides for the annual average
%                            % to include months from neighboring years.
%    months_per_annual       % months in the multi-year average 
%                            % different resolution finer than the months of
%                            % a year. 
% OUTPUTS:
%                            % All operators and parameter coefficient
%                            % arrays are the same size as input.
%    A_an                    % A operator annually averaged
%    H_an                    % H operator annually averaged
%    D_an                    % D operator annually averaged
%    T_an                    % T operator annually averaged
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get dimensions 
  iocn = MET.iocn;
  iocn_matrix_dim = length(iocn);

% Initialize the totals matrixes
  col_zeros   = (zeros(iocn_matrix_dim,1));
  Adv_total   = spdiags(col_zeros, 0, iocn_matrix_dim, iocn_matrix_dim);
  Hdiff_total = spdiags(col_zeros, 0, iocn_matrix_dim, iocn_matrix_dim);
  Dv_total    = spdiags(col_zeros, 0, iocn_matrix_dim, iocn_matrix_dim);
  T_op_total  = spdiags(col_zeros, 0, iocn_matrix_dim, iocn_matrix_dim);
  
  days_total  = 0;

  fperiod_ann = P.fperiod;
  nperiods = P.nperiods;
  lperiod_ann = fperiod_ann + nperiods -1;
 

% collect the totals
  for seq_num = fperiod_ann:lperiod_ann             % for periods in the annual average
     % get ops for the month
     op_file = (['MTM',int2str(seq_num),'.mat']);
     filename = [F.noc_dir, op_file];
     % disp ('loading'); disp(filename);
     load (filename); %  A H D T num_days
     days_total = days_total + num_days;
     T_op_total = T_op_total + num_days.*T; 
     Adv_total = Adv_total + num_days.*A;
     Hdiff_total = Hdiff_total + num_days.*H;
     Dv_total = Dv_total + num_days.*D;
  end % for months in multi-annual average

% Have totals, calculate average
  days_per_year = days_total./F.nyears;
  if days_per_year ~= P.days_per_yr
    disp('days do not add up, in annual_avg_dmo');
    keyboard
  end
  T_an     = T_op_total  ./ days_total;
  D_an     = Dv_total    ./ days_total;
  H_an     = Hdiff_total ./ days_total;
  A_an     = Adv_total   ./ days_total; 
  return
end % function annual_avg
%%%%%%%%%%%%%%%%%%%% 30 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ck_landfall %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [good] = ck_landfall(P,MET,mo, A, H, D, T) 
    good = true;                        % innocent unless proven guilty
    [ny,nx,nzz] = size(MET.MASK);
    KMT = MET.KMT;
    ocn_mask = zeros(ny,nx,nzz);        % make land mask without deletions for marginal seas
    for k = 1:nzz
       ocn_mask(:,:,k) = (k <= KMT(:,:));
    end % all k's 
    land_mask = ocn_mask == 0;
     
  % check for no D values on land (or seafloor) 
    D_on_land = D * land_mask(:);
    ck_D_on_land = max(abs(D_on_land)); % should be zeros
    if ck_D_on_land ~= 0
        disp([ 'D on land not zero for month ',int2str(mo)]);
        keyboard
        good = false;
    end
    
  % check for no A values on land (or seafloor) 
    A_on_land = A * land_mask(:);
    ck_A_on_land = max(abs(A_on_land)); % should be zeros
    if ck_A_on_land ~= 0
        disp( 'A on land not zero')
        good = false;
        keyboard
    end  
 
  % check for no H values on land (or seafloor) 
    H_on_land = H * land_mask(:);
    ck_H_on_land = max(abs(H_on_land)); % should be zeros
    if ck_H_on_land ~= 0
        disp( 'H on land not zero')
        good = false;
        keyboard
    end 
    
  % check for no T values on land (or seafloor) 
  % check T last; if there are problems with the A,H,or D operators.
  % want those to show up
    T_on_land = T * land_mask(:);
    ck_T_on_land = max(abs(T_on_land)); % should be zeros
    if ck_T_on_land ~= 0
        disp( 'T on land not zero')
        good = false;
        keyboard
    end
    
   return
end % function ck_landfall
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  30 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
