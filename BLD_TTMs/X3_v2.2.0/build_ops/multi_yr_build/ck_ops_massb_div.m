function [P] = ck_ops_massb_div(F,P,MET,C_period) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The purpose of ck_ops is to provide mass-balance and row-sum  checks on
% the building process, and allow graphical examination of the data when
% problems are detected.
% Checks A, H, and D operators for the specified period.
% mo is to be the CALENDAR period.
% Logs errors that are out of range, or that the check was successful.
%
% INPUTS
%    F                       % structure with file locations
%      ops_dir               % directory for the final operators
%    P                       % structure with control paramters
%      fignum                % first figure number to use
%      max_exp_err           % maximum expected error value.
%                            % Error analysis data and plots are 
%                            % output if this parameter is exceeded.
%    MET                     % structure with geometric parameters
%    C_period                % period number of the year (extension for MTM)
% OUTPUTS
%    P.fignum                % figure number is the next figure to use
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  disp(['checking operators for period ', C_period ]);
  
  fn = ([F.ops_dir,'MTM', C_period,'.mat']); % contains A,H,D,T,dxidt
  
  load (fn);
  % apply dxidt to A; note: dxidt is zeros for annual operator
 
  [good,P] = check_D (D,P,MET,C_period); if good, disp('D check OK'); end;
  [good,P] = check_H (H,P,MET,C_period); if good, disp('H check OK'); end;
  [good,P] = check_A (A,dxidt,P,MET,C_period); if good, disp('A check OK'); end;
 
  return
end % function ck_ops_massb_div
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% check_H %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [good,P] = check_H(H,P,MET,C_period)
% check the horizontal diffusion operator
% 
% INPUTS
%    H                       % horizontal diffusion operator
%    P                       % structure with control paramters
%      fignum                % first figure number to use
%      max_exp_err           % maximum expected error value.
%                            % Error analysis data and plots are 
%                            % output if this parameter is exceeded.
%    MET                     % structure with geometric parameters
%    C_period                      % period number of the year
% OUTPUTS
%    P.fignum                % figure number is the next figure to use
%    good                    % true if the check was successful.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  disp(['checking H for period ',C_period])
  good = true;                    
  fignum = P.fignum;
  [ny,nx,nz] = size(MET.MASK);
  iocn = MET.iocn;
  dVt = MET.VOL(iocn);
   
  % Proportional mass (col) error calculation
  Sum_H = ((H') * dVt(:)) ./ dVt(:);  % should be zeros'
  max_err_H_col = max(abs(Sum_H(:)));
  disp(['H max mass col error = ',num2str(max_err_H_col)]);
  
  % keyboard % force error: set max_err_H_col above P.max_exp_err
  if (max_err_H_col >= P.max_exp_err) % only do extra calc if err is significant
      good = false;
  
      % Make a distribution for a bar chart
  %    op_name = 'H';
  %    item_name = 'column sum err';
  %    fignum = bar_chart(max_err_H_col,Sum_H,MET,fignum, op_name, ...
  %                item_name, C_period);
              
  %   % Show the locations of the worst error values 
  %    H_sums_3d = zeros(ny,nx,nz) +  NaN;
  %    H_sums_3d(iocn) = Sum_H;
      
  %    [fignum] = show_worst_errs(max_err_H_col,H,Sum_H,H_sums_3d, ...
  %                             MET,fignum, op_name, item_name, C_period);
    
  end % if error greater than specified

  % H row sum errors check
  H_m = H * (0*iocn+1); % [delta conc per second]  
  H_row_max = max(abs(H_m(:)));
  disp(['H row max = ',num2str(H_row_max)]);
    
   % keyboard % force error; reset H_row_max above P.max_exp_err
   if (H_row_max >= P.max_exp_err) % do other stuff if significant error 
      good = false;
   %    H_m_3d = zeros(ny,nx,nz) + NaN;
   %    H_m_3d(iocn) = H_m;
      
      % Make a distribution for a bar chart
   %   op_name = 'H';
   %  item_name = 'row sum err';
   %  fignum = bar_chart(H_row_max,H_m(:),MET,fignum, op_name, ...
   %               item_name, C_period);
      
       % show worst row sum errors
   %   fignum = ...x
   %      show_worst_rows(H_row_max, H, H_m_3d(:), H_m_3d, ...
   %                     MET,fignum, op_name, ...
   %                     item_name, C_period);
                    
    end % if max err > exp error level   
  P.fignum = fignum;
  return
end % function check_H

%% check_A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [good,P] = check_A(A,dxidt,P,MET,C_period)
% check the horizontal diffusion operator
% 
% INPUTS
%    A                       % advection operator
%    P                       % structure with control paramters
%      fignum                % first figure number to use
%      max_exp_err           % maximum expected error value.
%                            % Error analysis data and plots are 
%                            % output if this parameter is exceeded.
%    MET                     % structure with geometric parameters
%    C_period                      % period number of the year
% OUTPUTS
%    P.fignum                % figure number is the next figure to use
%    good                    % true if the check was successful.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  disp(['checking A for period ',C_period]);
  good = true;
  fignum = P.fignum;
  [ny,nx,nz] = size(MET.MASK);
  iocn = MET.iocn;
  dVt = MET.VOL(iocn);
 
  % check sums of columns 
  % Proportional mass (col) error calculation
  Sum_A = ((A') * dVt(:)) ./ dVt(:);  % should be zeros
  max_err_A_col = max(abs(Sum_A(:)));
  disp(['A max mass col error = ',num2str(max_err_A_col)]);
  
  % keyboard % force error: set max_err_A_col > P.max_exp_err
    if (max_err_A_col >= P.max_exp_err) % only do extra calc if err is significant
      good = false;
      
      % Make a distribution for a bar chart
%      op_name = 'A';
%      item_name = 'column sum err';
%     fignum = bar_chart(max_err_A_col,Sum_A,MET,fignum, op_name, ...
%                 item_name, C_period);
%              
%      % Make a set of color distribution chart for err col 
%     A_sums_3d = zeros(ny,nx,nz) +  NaN;
%      A_sums_3d(iocn) = Sum_A;
      
%      [fignum] = show_worst_errs(max_err_A_col,A,Sum_A,A_sums_3d, ...
%                               MET,fignum, op_name, item_name, C_period);
    
    end  % if error greater than specified

    % row sum errors 
    A_m = A * (0*iocn+1); % [delta conc per second]   
    A_row_max = max(abs(A_m(:)));
    disp(['A row max = ',num2str(A_row_max)]);

   % keyboard  % force execution: set A_row_max > P.max_exp_err
   if (A_row_max >= P.max_exp_err) % do other stuff if significant error 
      good = false;
      A_m_3d = zeros(ny,nx,nz) + NaN;
      A_m_3d(iocn) = A_m;
      
     % Make a distribution for a bar chart     
     op_name = 'A';
     item_name = 'row sum err';
     fignum = bar_chart(A_row_max,A_m(:),MET,fignum, op_name, ...
                 item_name, C_period);
     % show the worst row errors          
      fignum = ...
        show_worst_rows(A_row_max, A, A_m_3d(:), A_m_3d, ...
                       MET,fignum, op_name, ...
                       item_name, C_period);
%     
   end % if max err > exp error level
    
  P.fignum = fignum;
  return
end % function check_A

%% Check_D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [good,P] = check_D(D,P,MET,C_period)
% check out the vertical diffusion operator
% 
% INPUTS
%    D                       % vertical diffusion operator
%    P                       % structure with control paramters
%      fignum                % first figure number to use
%      max_exp_err           % maximum expected error value.
%                            % Error analysis data and plots are 
%                            % output if this parameter is exceeded.
%    MET                     % structure with geometric parameters
%    C_period                      % period number of the year
% OUTPUTS
%    P.fignum                % figure number is the next figure to use
%    good                    % true if the check was successful.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  disp(['checking D for period ',C_period])
  good = true;
  fignum = P.fignum;
  [ny,nx,nz] = size(MET.MASK);
  iocn = MET.iocn;
  dVt = MET.VOL(iocn);
  
  % check sums of columns 
  % Proportional mass (col) error calculation for DV
  Sum_D = ((D') * dVt(:)) ./ dVt(:);  % should be zeros
  max_err_D_col = max(abs(Sum_D(:)));
  disp(['D max mass col error = ',num2str(max_err_D_col)]);
  disp('D in units of 1/sec');
  
   if (max_err_D_col >= P.max_exp_err) % only do extra calc if err is significant
      good = false;
%      % Make a distribution for a bar chart
%      op_name = 'D';
%      item_name = 'column sum err';
%      fignum = bar_chart(max_err_D_col,Sum_D,MET,fignum, op_name, ...
%                  item_name, C_period);
%              
%      % Make a set of color distribution chart for err col D 
%      D_sums_3d = zeros(ny,nx,nz) +  NaN;
%      D_sums_3d(iocn) = Sum_D;
%   
%      [fignum] = show_worst_errs(max_err_D_col,D,Sum_D,D_sums_3d, ...
%                               MET,fignum, op_name, item_name, C_period);
%    
    end  % if error greater than specified
 
    % D row sum errors 
    D_m = D *(0*iocn+1); % [delta conc per second]
   
    D_row_max = max(abs(D_m(:)));
    disp(['D row max = ',num2str(D_row_max)]);
    % keyboard  % force execution
    if (D_row_max >= P.max_exp_err) % do other stuff if significant error       
       good = false;
%       D_m_3d = zeros(ny,nx,nz) + NaN;
%       D_m_3d(iocn) = D_m;
%       
%      % Make a distribution for a bar chart
%      op_name = 'D';
%      item_name = 'row sum err';
%      fignum = bar_chart(D_row_max,D_m(:),MET,fignum, op_name, ...
%                  item_name, C_period);
%              
%      % Make a set of color distribution chart for err col D   
%      [fignum] = show_worst_rows(D_row_max,D,D_m_3d(:),D_m_3d, ...
%                               MET,fignum, op_name, item_name, C_period);  
    end % if max err > exp error level
 
  P.fignum = fignum;
  return
end % function check_D

%% bar_chart function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [next_fignum] = ...
    bar_chart(tot_max_err,opsum_1D,MET,fignum, op_name, item_name, C_period)   
% makes a bar chart for the distribution of the sum errors
% 
% INPUTS:
%     total_max_err          % the maximum error found for the entire ocean
%     opsum_1D               % 1-D (vector) of the sum showing the error values
%     MET                    % geometric data structure,
%     fignum                 % first figure number to use,
%     op_name                % text name of the operator for titles on displays
%     item_name              %text name of the caluculation, such as 'row sums', 
%     C_period                     % integer period number
%
% OUTPUTS:
%     next_fignum            % next figure number to be used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    iocn = find(MET.MASK);
    [ny,nx,nz] = size(MET.MASK);
    matrix_dim = length(iocn);
    
    dist = zeros(25,1);
    x_scale = tot_max_err/24;
    x_dist = [0:x_scale:tot_max_err]';
    for m_i = 1:matrix_dim
        d_i = round(abs(opsum_1D(m_i)/x_scale)) + 1;
        dist(d_i) = dist(d_i) + 1;
    end
    figure (fignum)
    bar(x_dist,dist);
    title ([' Dist of ',op_name, ' ', item_name, ...
            ' for period ',C_period]);
   next_fignum = fignum+1;
   return
end %function barchart

%% show_worst_errs function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function next_fignum = ...
         show_worst_errs(total_max_err,op,opsum_1D,op_sums_3d, ...
                        MET,fignum, op_name, ...
                        item_name, C_period)
% shows graphical and text information about the largest errors in the 
% column sums.
%
% INPUTS:
%     total_max_err          % maximum error found for the entire ocean
%     op                     % operator being analyzed, full-sized
%     opsum_1D               % 1-D (vector) of the sum showing the error values
%     op_sums_3d             % 3-D version of the sum showing the error values
%     MET                    % geometric data structure
%     fignum                 % first figure number to use
%     op_name                % text name of the operator for titles on displays
%     item_name              % text name of the caluculation, such as 'col sums' 
%     C_period                     % integer period number
% OUTPUTS:
%     next_fignum            % next figure number to be used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  iocn = find(MET.MASK);
  [ny,nx,nz] = size(MET.MASK);
    
  % find the worst error values  
  %  for column totals  %%
  display_list = zeros(nz,1);        % to collect levels with worst values
  disp(['worst ',op_name,'  ',item_name]);
  for worst = 1:5
      [max_val,column_I]= max(abs(opsum_1D(:)));
      nzr   = find(op(:,column_I)); % find nz entries in this column
      n_nzr = length(nzr);  % number of nz entries
      fmtnum = '%+8.4e';
      disp(['for error value ',num2str(max_val),' op index ', ...
            int2str(column_I)]);
      disp(' j   i   k   op_index   op value');
      for en = 1:n_nzr
         % determine the 3-D location
         i_3D = iocn(nzr(en));
         kk(en) = floor((i_3D-1) / (ny*nx)) + 1;
         ii(en) = floor(((i_3D-1) - ((kk(en)-1)*ny*nx))/ny) + 1;
         jj(en) = (i_3D-1) - ((kk(en)-1)*ny*nx) - ((ii(en)-1)*ny) + 1;
         disp([int2str(jj(en)),'  ', ... 
               int2str(ii(en)),'  ', ... 
               int2str(kk(en)),'  ', ...
               int2str(nzr(en)),'  ', ... 
               num2str(op(nzr(en),column_I),fmtnum) ]);
         if column_I == nzr(en)        % if this is the error location
            display_list(kk(en)) = 1;  % mark for display
         end
      end % for all nz rows
      disp (' ');
      opsum_1D(column_I) = 0; % eliminate the worst one
  end % for worst n
    
  % output displays for worst levels
  for k = 1:nz 
      if display_list(k) > 0
        op_sums_2d = op_sums_3d(:,:,k);
        op_disp = op_sums_2d ./total_max_err;   % scale to worst error
        maxerr = max(abs(op_sums_2d(:)) );
        figure (fignum)       
        pcolor (op_disp); shading flat ;
        caxis([-1 +1]); colorbar;       
        title([op_name,' ', item_name,' k=',int2str(k),', ' ...
                'period ',C_period, ...
                '. COLOR SCALED to MAX ERR']);
        maxerr_k = max(abs(op_sums_2d(:)) );
        xlabel(['maxerr level ',int2str(k),' = ',num2str(maxerr_k)]);
        fignum = fignum+1;
      end % if this level had a maximum value
  end % for k levels
  next_fignum = fignum+1;
  return
end %function show_worst_errs

%% show_worst_rows function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function next_fignum = ...
         show_worst_rows(total_max_err,op,opsum_1D,op_sums_3d, ...
                        MET,fignum, op_name, ...
                        item_name, C_period)
% shows graphical and text information about the largest errors in the row 
% sums
%
% INPUTS:
%     total_max_err          % maximum error found for the entire ocean
%     op                     % operator being analyzed, full-sized
%     opsum_1D               % 1-D (vector) of the sum showing the error values
%     op_sums_3d             % 3-D version of the sum showing the error values
%     MET                    % geometric data structure
%     fignum                 % first figure number to use
%     op_name                % text name of the operator for titles on displays
%     item_name              % text name of the caluculation, such as 'col sums' 
%     C_period                     % integer period number
% OUTPUTS:
%     next_fignum            % next figure number to be used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  iocn = find(MET.MASK);
  [ny,nx,nz] = size(MET.MASK);
    
  % find the worst error values  
  %  for row totals  %%
  display_list = zeros(nz,1);        % to collect levels with worst values
  disp(['worst ',op_name,'  ',item_name]);
  for worst = 1:5
      [max_val,row_I]= max(abs(opsum_1D(:)));
      nzr   = find(op(row_I,:)); % find nz entries in this row
      n_nzr = length(nzr);  % number of nz entries
      disp(['for error value ',num2str(max_val),' op index ', ...
            int2str(row_I)]);
      disp(' j   i   k   op index   op value');
      fmtnum = '%+8.4e';
      for en = 1:n_nzr
         % determine the 3-D location
         i_3D = iocn(nzr(en));
         kk(en) = floor((i_3D-1) / (ny*nx)) + 1;
         ii(en) = floor(((i_3D-1) - ((kk(en)-1)*ny*nx))/ny) + 1;
         jj(en) = (i_3D-1) - ((kk(en)-1)*ny*nx) - ((ii(en)-1)*ny) + 1;
         disp([int2str(jj(en)), '  ', ... 
               int2str(ii(en)), '  ', ... 
               int2str(kk(en)), '  ', ...
               int2str(nzr(en)),'  ', ... 
               num2str(op(row_I,nzr(en)),fmtnum)]);
         if row_I == nzr(en)        % if this is the error location
            display_list(kk(en)) = 1;  % mark for display
         end
      end % for all nz rows
      disp (' ');
      opsum_1D(row_I) = 0; % eliminate the worst one
  end % for worst n
    
  % output displays for worst levels
  for k = 1:nz 
      if display_list(k) > 0
        op_sums_2d = op_sums_3d(:,:,k);
        op_disp = op_sums_2d ./total_max_err;   % scale to worst error
        
        figure (fignum)       
        pcolor (op_disp); shading flat ;
        caxis([-1 +1]); colorbar;       
        title([op_name,' ', item_name,' for k level = ',int2str(k),', ' ...
                'for period ',C_period, ...
                '. COLOR SCALED to MAX ERR']);
        maxerr_k = max(abs(op_sums_2d(:)) );
        xlabel(['maxerr level ',int2str(k),' = ',num2str(maxerr_k)]);
        fignum = fignum+1;
      end % if this level had a maximum value
  end % for k levels
  next_fignum = fignum+1;
  return
end %function show_worst_rows
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  30 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

