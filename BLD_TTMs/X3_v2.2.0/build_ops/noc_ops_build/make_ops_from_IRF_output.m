function make_ops_from_IRF_output(MET,F,P, year, period)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The purpose of make_ops is to build the advection/diffusion operators
% for the offline model.  The program extracts the advection and diffusion
% tendencies from the history files resulting from running
% the POP ocean model with the "IRF tracers" turned on. It outputs a
% file containing a sparce matrix with each operator(A (advection), 
% H (horizontal diffusion), D (vertical diffusion), and T (total transport
% matrix) for the selected period.
%
% The archived history files are organized such that one file contains the
%  data for 1 month for all the regular points (g-mask) + 
% ovfl_entrainment (e_mask) + ovfl_source(s-mask).
%
% Each regular pulse is one of 125 "tracers" defined, allowing for spacial 
% separation of +2 and -2 grid boxes in each direction, so that its 
% contribution to the neighboring cells can be captured. This results 
% in a 5 x 5 x 5 box around each cell.
%
% There are two overflow special handling types of tracers, 
% associated with sill overflow special processing: one
% associated with the "source" region, behind the sill, and one associated
% with the entrainment region, partway down the overflow side of the sill.
% The advection of tracers in each of these areas extends many grid cells 
% away in a "product" region, near the base of the sill.  The source and 
% entrainment regions have 33 and 20 tracers associated with them, 
% respectively.
%
% The data from the set of history files is picked off and reorganized,
% such that it is possible to obtain all the contributions to a particular
% grid cell.  This is organized into a sparce matrix, with the impulse 
% locations indicated by the column index.  The row index within a column
% corresonds to the response location.  Each IRF corresponds to one column
% in the matrix.   This is written out to a file, one per month. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ** Note that the IRF-mask files are required to be pre-built to execute
%  this program, which builds the raw monthly operators. The IRF-mask files
%  are listed as part of the input variables listed below.
%
%  make_ops_from_IRF_output is  an outer layer that takes care of 
%    -- looping through the source files to collect the data;
%    -- invokes functions that get the specific data needed for
%       each operator;
%    -- combines the operators and saves them in files.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  persistent g_irf_mask s_irf_mask e_irf_mask IRF_ncid IRFdir

  % input parameters:
  % MET                     % structure with geometric grid definitions
  % F                       % structure with file location definitions
  % P                       % structure with control parameters

  [ny,nx,nzz] = size(MET.MASK);

  % additional geometry   
  matrix_dim = nx*ny*nzz;  % sparce matrix dimensions

  index_m=[1:matrix_dim]'; % used in get_x_data to calculate the linear 
                             % index of a value in the matrix'
                            
  % Get the irf_masks that separate the impulses computed in parallel.
  if (isempty(g_irf_mask))
     filename = [F.irf_stencil_dir,F.g_irf_mask_file];
     g = eval('load(filename)');
     g_irf_mask = g.g_irf_mask;
     clear g;
     if F.s_ntracers > 0
        filename = [F.irf_stencil_dir,F.s_irf_mask_file];
        s = eval('load(filename)');
        s_irf_mask = s.s_irf_mask;
        clear s;

        filename = [F.irf_stencil_dir,F.e_irf_mask_file];
        e = eval('load(filename)');
        e_irf_mask = e.e_irf_mask;
        clear e;
     end % if overflows
  end % if isempty

  % set up the additional geometric variables in structure GEO_PARAM 
  GEO_PARAM.ny = ny; GEO_PARAM.nx = nx; GEO_PARAM.nzz = nzz;
  GEO_PARAM.matrix_dim  = matrix_dim;
  GEO_PARAM.index_m     = index_m;
  GEO_PARAM.g_irf_mask  = g_irf_mask;
  GEO_PARAM.s_irf_mask  = s_irf_mask;
  GEO_PARAM.e_irf_mask  = e_irf_mask;
%
%% Data extraction processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Build one file per PERIOD for each operator.

  % open the IRF pulse file - same for all periods
  if isempty(IRF_ncid)   
     IRF_ncid = netcdf.open([F.IRFdir,F.IRFfile],'NC_NOWRITE');
  end % if isempty
  FILE_PARAM.IRF_ncid = IRF_ncid;

   seq_num = (year - F.firstyr)*F.num_periods + period;                   
   file_yr = year;  % source history file year
    
    % set up the filename for this period
    if F.day5
       [cal_mo,cal_day] = get_index_5day(period);
       data_file = sprintf('%s%04i-%02i-%02i.nc',F.g_file,F.file_yr,cal_mo,cal_day);
    else
       data_file = sprintf('%s%04i-%02i.nc',F.g_file,file_yr,period);
    end
        
    disp ('getting data from '); disp(['  ',F.g_dir,data_file]);
    FILE_PARAM.ncid = netcdf.open([F.g_dir,data_file],'NC_NOWRITE');
   
    % Retrieve the advection and diffusion fields and parameters 
    % from various files of the OGCM run, 
    % using the appropriate mask for the data type,
    % and form into the matrix operators
 
    for IRF_case = P.fcase:P.lcase 
     switch IRF_case

      case 1  % k_level = all %%%%%%%%%%%%%%%%%%%%%% 
       Adv_g  = sparse(matrix_dim, matrix_dim);     % start with zero
       Hdif_g = sparse(matrix_dim, matrix_dim);     % start with zero
       for klev = 1:5   
        FILE_PARAM.klev = klev;
        %disp(['klev = ',int2str(FILE_PARAM.klev)]);
        FILE_PARAM.ntracers = F.k_ntracers/5; % number of tracers for this case       
        % extract the data from the history file, return the vectors
        % for the current klev for making the sparse matrix.
        [icol,irow,a_val,h_val,Pulse_kx] = ...
          get_g_data(GEO_PARAM, FILE_PARAM, F, P.collect_pulse);
            
        % Build the sparse matrixes for this klev for this period        
        Adv_gx = sparse(irow,icol,a_val,matrix_dim,matrix_dim);                 
        Hdif_gx = sparse(irow,icol,h_val,matrix_dim,matrix_dim);
        
        Adv_g = Adv_g + Adv_gx;
        Hdif_g = Hdif_g + Hdif_gx;
           
        % Save the Pulse matrix
        if P.collect_pulse               % if the impulse is to be collected   
            eval(['save ',F.tmp_ops_dir,'Pulse_k',int2str(klev), ...
                  '.mat',' Pulse_kx -v7.3;']);
            disp([F.tmp_ops_dir,'Pulse_k',int2str(klev),'.mat saved']);
        end % if collect pulse
       end % for each klev 
       
      case 2  % sill overflow source data %%%%%%%%
         if F.s_ntracers > 0
            FILE_PARAM.klev = 1; % to structure the data pick up  
            FILE_PARAM.ntracers = F.s_ntracers; % number of tracers for this case      
            % extract the data from the history file,
            % add to the vectors for the current klev.
   
            [icol,irow,a_val_s,h_val_s,Pulse_s] = ...
               get_s_data(GEO_PARAM, FILE_PARAM, F, P.collect_pulse);
      
            % Build the sparse matrixes for this klev for this period 
            Adv_s = sparse(irow,icol,a_val_s,matrix_dim,matrix_dim);      
            Hdif_s = sparse(irow,icol,h_val_s,matrix_dim,matrix_dim);
    
            % Save the Pulse matrix
            if P.collect_pulse               % if the impulse is to be collected   
              eval(['save ',F.tmp_ops_dir,'Pulse_s.mat Pulse_s -v7.3;']);
              disp([F.tmp_ops_dir,'Pulse_s.mat saved']);                  
            end % if collect pulse
         end % if overflows
       
      case 3  % sill overflow entrainment data %%%
         if F.s_ntracers > 0
            FILE_PARAM.klev = 1;
            % disp(['sill overflow entrainment data processing']);     
            FILE_PARAM.ntracers = F.e_ntracers; % number of tracers for this case 
      
            % extract the data from the history file,
            % add to the vectors for the current klev. 
            [icol,irow,a_val_e,h_val_e,Pulse_e] = ...
               get_e_data(GEO_PARAM, FILE_PARAM, F, P.collect_pulse);
            Adv_e = sparse(irow,icol,a_val_e,matrix_dim,matrix_dim);    
            Hdif_e = sparse(irow,icol,h_val_e,matrix_dim,matrix_dim);
  
            % Save the Pulse matrix
            if P.collect_pulse                 % if the impulse to be collected   
               eval(['save ',F.tmp_ops_dir,'Pulse_e.mat Pulse_e -v7.3;']);
               disp([F.tmp_ops_dir,'Pulse_e.mat saved']);
            end % if collect pulse
         end % if overlows
         
      case 4  % VDC data %%%%%%%%%%%%%%%%%%%%%%%%%
        % disp(['Vertical Diffusion data processing']);
       
        % extract the data from the history file,
        [KTF, KBF] = get_VDC_data(MET,FILE_PARAM,GEO_PARAM);                 
        % [KTF, KPP_NLK] = get_VDC_data(FILE_PARAM,GEO_PARAM);
        
        % (obsolete?) Save KBF matrix of VDC data
        % updated value might be used in par_tridiagonal
        % out_file_name = ([F.tmp_ops_dir,'KBF_full_mo_',C_mo,'.mat']);
        % eval(['save ',out_file_name,' KBF  -v7.3;']);  
        % disp([ out_file_name,' saved']);
      
        % Build the VDC operator matrix      
        [D] = Vdiff(MET,GEO_PARAM,KTF,KBF);
      
      otherwise % "there is no other hand" %%%%%%
        
        disp (['unexpected case number ',int2str(IRF_case)]);
        keyboard
     end  % switch on IRF case   

    end   % loop on IRF cases
        
%% Combine Adv matrices to make operators %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% After all cases have extracted the tracer data for this period,
%  we are ready to combine the sparce matrixes for the operators.
%  For adv, build one for the regular k_levels, and 1 each for source and 
%  entrainment, then add them together. 
%  For hdif, build one. Since the data was acquired in the same manner
%  as for adv, follow the same building approach.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Make a complete Adv sparse matrix including the overflow regions
    % disp('Building Adv sparse matrix');
    A = Adv_g;
    if F.s_ntracers > 0
       A = A + Adv_s + Adv_e;
    end
          
    % Make a complete Hdif sparse matrix including the overflow regions
    % disp('Building Hdiff sparse matrix');
    H = Hdif_g;
    if F.s_ntracers > 0
       H = H + Hdif_s + Hdif_e;   
    end

    % Put together the total Advection-diffusion operator
    % disp('Building T, total Advection-Diffusion operator');
    T = A + H + D;
    if F.day5
        num_days = 5;
    else
        num_days = P.days_per_mo(period);
    end

    if P.collect_pulse
       [good] = pulse_ck(P,MET,F); if good, disp( 'Pulse Check OK'); end;  
       P.collect_pulse = false;           % only need to do this for 1 pass
       
       % check that there were no land points with operator values   
       [good] = ck_landfall(F,P,MET,seq_num, A, H, D, T); 
       if good, disp('ops landfall OK'); end
    end % if collecting pulse
    
    % make iocn-sized operators, with marginal-seas eliminated mask
    % Note MATLAB is much faster doing calculations by column than by rows.
    iocn = MET.MASK(:) > 0;
    At = (A(:,iocn))';     % mask off 
    A = (At(:,iocn))';
  
    Ht = (H(:,iocn))';     % mask off 
    H = (Ht(:,iocn))';
  
    Dt = (D(:,iocn))';     % mask off 
    D = (Dt(:,iocn))';
 
    
    Tt = (T(:,iocn))';     % mask off 
    T = (Tt(:,iocn))';
 
    out_file_name = ([F.noc_dir,'MTM',int2str(seq_num),'.mat']);
    eval(['save ',out_file_name,' A H D T num_days  -v7.3;']);
    disp([ out_file_name,' saved']);
 
    % close history file for the period and remove from local storage
    netcdf.close(FILE_PARAM.ncid);
    
  return
end % make_ops_from_IRF_output

%% ck_landfall %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [good] = ck_landfall(F,P,MET,seq_num, A, H, D, T) 
    good = true;                        % innocent unless proven guilty
    [ny,nx,nzz] = size(MET.MASK);
    ms_MET = ms_buildMET(F); % get mask without deletions 
    ocn_mask = ms_MET.MASK;
    land_mask = ocn_mask == 0;
     
  % check for no D values on land (or seafloor) 
    D_on_land = D * land_mask(:);
    ck_D_on_land = max(abs(D_on_land)); % should be zeros
    if ck_D_on_land ~= 0
        disp([ 'D on land not zero for period ',int2str(seq_num)]);
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
