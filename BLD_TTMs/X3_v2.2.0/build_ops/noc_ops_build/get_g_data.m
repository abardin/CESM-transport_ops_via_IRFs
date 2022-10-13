function [icol,irow,a_val,h_val,Pulse] = ...
    get_g_data(GEO_PARAM, FILE_PARAM, F, collect_pulse)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Retrieves the data for the global IRF tracers from the history file.
% The data are organized into row, column, and value data, ready to
% build the sparce matrixes for the operators.
%
% INPUTS:
%   GEO_PARAM:
%      ny, nx, nzz           % 3-D dimensions for the full grid size
%      g_irf_mask            % mask for regular irf pulse and responses
%      index_m               % linear index of a value in the matrix
%   FILE_PARAM:
%      ncid                  % key identifying the opened history file
%      klev                  % current k-lev for the history file (always
%                              1)
%      ntracers              % the number of tracers in the file
% 
% OUTPUTS:
%      icol                  % column index vector
%      irow                  % row index vector
%      a_val                 % corresponding value in A (advection) matrix
%                            % for position given by icol, irow
%      h_val                 % corresponding value in H (horizonal 
%                            % diffusion) matrix   
%      Pulse                 % 3-D locationss of the pulse's (as 1's) 
%                            % associated with this data.
%      Note: The additional level of land has been added to the z-dimension
%            of the parent POP model.
%      Note: Output is in YXZ order (Z varies slowest) in the Matlab
%            sparse matrix transport operators.  In CESM everything is in
%            XYZ.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  persistent MASKms

  ny = GEO_PARAM.ny;  nx = GEO_PARAM.nx; nzz = GEO_PARAM.nzz; nz = nzz-1;
  g_irf_mask = GEO_PARAM.g_irf_mask;
  index_m    = GEO_PARAM.index_m;
  % 
  ncid      = FILE_PARAM.ncid;
  kr        = FILE_PARAM.klev;
  ntracers  = FILE_PARAM.ntracers;
  if ntracers ~= 25
      disp('unexpected number of tracers for g_data')
      keyboard
  end
  
  IRF_ncid  = FILE_PARAM.IRF_ncid;
  
  % get the land mask including marginal seas
  if isempty(MASKms)
     METms = ms_buildMET(F);
     MASKms = METms.MASK;
  end

  % Initiate the column, row, and value vectors.
  icol =  [ ];
  irow =  [ ];
  a_val = [ ]; % for adv
  h_val = [ ]; % for hdif
  
  % Initiate collections for processing / tracer check
  Pulse = zeros(ny,nx,nzz);
  
  % loop through for each tracer in the file
  % three indices are necessary to get the names of the tracers, but the
  % third one (kr = klev) is constant for a given output file

  for ir = 1:5 % loop through for each tracer in the file
    for jr = 1:5 % second indexing necessary to get the 
                            % names of the tracers
       
      % extract the IRF impulse locations **
     
      cmd = sprintf ...
          ('varid = netcdf.inqVarID(IRF_ncid,''IRF_%i_%i_%i'');',ir,jr,kr);
      eval(cmd);
      IRF_d(:,:,:) = netcdf.getVar(IRF_ncid,varid);  
     
      IRF_p =permute ( IRF_d(:,:,:),[2,1,3]);  % change the order
                                               % (i,j,k)to (j,i,k) and store in
                                               % double precision
      IRF_p = (IRF_p .* (IRF_p < 1.0e36));     % eliminate the no data indicator
      IRF = logical(zeros(ny,nx,nzz)); 
      IRF(:,:,1:nz) = logical(IRF_p(:,:,:));   % resize with added z layer
      
      % remove the pulse locations on land 
      IRF = IRF .* MASKms;

      % advection %%    
 
      cmd = sprintf ...
      ('varid = netcdf.inqVarID(ncid,''ADV_3D_IRF_%i_%i_%i'');',ir,jr,kr);
      eval(cmd);
      var_d(:,:,:) = netcdf.getVar(ncid,varid);
  
      var_p =double(permute ( var_d(:,:,:),[2,1,3]));  % change the order
                                          % (i,j,k)to (j,i,k) and store in
                                          % double precision
      var_p = var_p .* (var_p < 1.0e36);  % eliminate the nodata indicator
      var = zeros(ny,nx,nzz);
      var(:,:,1:nz) = var_p(:,:,:);       % resize with added z layer
    
      % horizontal diffusion %%
     
      cmd = sprintf ...
      (['varid = netcdf.inqVarID(ncid,', ...
               '''HDIF_EXPLICIT_3D_IRF_%i_%i_%i'');'],ir,jr,kr);
      eval(cmd);
      hdif_d(:,:,:) = netcdf.getVar(ncid,varid);
    
      hdif_p =double(permute ( hdif_d(:,:,:),[2,1,3]));  % change the order
                                            % (i,j,k)to (j,i,k) and store
                                            % in double precision
      hdif_p = hdif_p .* (hdif_p < 1.0e36); % eliminate nodata indicator
      hdif = zeros(ny,nx,nzz);
      hdif(:,:,1:nz) = hdif_p(:,:,:);       % resize with added z layer
     
      % Use the IRF pulse locations to locate the piece of the irf_mask
      % to associate the data with the pulse
                                  
      col = g_irf_mask*(IRF(:).*index_m(:));
      inz = find(col);
      % find nonzero values of col
      max_col = max(col(inz(:)));
      if max_col > GEO_PARAM.matrix_dim
          disp(['col index too big ',int2str(max_col) ])
          keyboard
      end
      % append the new values and corresponding row and column
      % indices to the list for adv and hdif
      icol  = [icol; col(inz(:)) ];         % impulse location index
      irow  = [irow; index_m(inz(:)) ];     % response location  index
      a_val = [a_val; var(inz(:)) ];        % response value for adv
      h_val = [h_val; hdif(inz(:)) ];       % response value for hdif
 
      % check out IRF pulse  %%
      if collect_pulse % if to check out the pulse locations
         % have a look
         % col_3D = reshape(col,ny,nx,nzz);
     
         % Locate each pulse, record pulse; record frequency of each 
         % related cell having a non-zero value; provides for checking
         % later
      
         for z_IRF = 1:nzz
           for x_IRF = 1:nx
             for y_IRF = 1:ny
               if (IRF(y_IRF,x_IRF,z_IRF) ~= 0)
                   Pulse(y_IRF,x_IRF,z_IRF) = Pulse(y_IRF,x_IRF,z_IRF) + 1;
                   % to check that there is 1 and only 1 pulse each loc in
                   % ocean              
               end % if IRF pulse
             end % for y_IRF = 1:ny
           end % for x_IRF = 1:nx
         end % for z_IRF = 1:nzz
      end % if pulse collection wanted
      % end checout IRF pulse  %%
           
    end % for jr
  end % for ir

  % icol_size = size(icol)     % show number of matrix points
  % irow_size = size(irow)
  return
end % function get_g_data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  30   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
