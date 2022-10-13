function [icol_e,irow_e,a_val_e,h_val_e,Pulse] = ...
          get_e_data(GEO_PARAM, FILE_PARAM, F, collect_pulse)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Retrieves the data for the over-the-sill  overflow (OVF) entrainment region
% IRF tracers from the history file for the A (advection) and H (horizontal
% diffusion operators.
% The output data are organized into row, column, and value data, ready to
% build the sparce matrixes for the operators.
%
% INPUTS:
%   GEO_PARAM:
%      ny, nx, nzz           % 3-D dimensions for the full grid size
%      s_irf_mask            % mask for OVF source regions irf pulse and 
%                            % responses
%      index_m               % linear index of a value in the matrix
%   FILE_PARAM:
%      ncid                  % key identifying the opened history file
%      klev                  % current k-lev for the history file
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
  ny = GEO_PARAM.ny; nx = GEO_PARAM.nx; nzz = GEO_PARAM.nzz; nz = nzz-1;
  e_irf_mask = GEO_PARAM.e_irf_mask;
  index_m    = GEO_PARAM.index_m;
  ntracers   = FILE_PARAM.ntracers;
  ncid       = FILE_PARAM.ncid;
  IRF_ncid = FILE_PARAM.IRF_ncid;
  
  % get the land mask including marginal seas
  if isempty(MASKms)
     METms = ms_buildMET(F);
     MASKms = METms.MASK;
  end
  
  % Initiate the collection for entrainment data 
  icol_e =     [ ];
  irow_e =     [ ];
  a_val_e =    [ ]; % for adv
  h_val_e =    [ ]; % for hdif
  
  % Initiate collections for processing / tracer check
  Pulse = zeros(ny,nx,nzz);

  for ir = 1:ntracers                  % for each tracer
    
    % extract the IRF impulse locations %%    
    eval( sprintf ...
          ('varid = netcdf.inqVarID(IRF_ncid,''IRF_OVF_ENT_%02i'');',ir) );    
    IRF_d(:,:,:) = netcdf.getVar(IRF_ncid,varid);     
    IRF_p =permute ( IRF_d(:,:,:),[2,1,3]);  % change the order
                                             % (i,j,k)to (j,i,k) 
    IRF_p = IRF_p .* (IRF_p < 1.0e36);       % eliminate no data indicator
    IRF = zeros(ny,nx,nzz);
    IRF(:,:,1:nz) = logical(IRF_p(:,:,:));   % resize with added z layer                                   
    
    % remove the pulse locations on land 
      IRF = IRF .* MASKms;
    % advection %%
   
    eval( sprintf ...
      ('varid = netcdf.inqVarID(ncid,''ADV_3D_IRF_OVF_ENT_%02i'');',ir) );    
    var_d(:,:,:) = netcdf.getVar(ncid,varid);
    var_p =double(permute ( var_d(:,:,:),[2,1,3]));  % change the order
                                          % (i,j,k)to (j,i,k) and store in
                                          % double precision
    var_p = var_p .* (var_p < 1.0e36);    % eliminate the nodata indicator
    a_var = zeros(ny,nx,nzz);
    a_var(:,:,1:nz) = var_p(:,:,:);       % resize with added z layer
                                              
    % horizontal diffusion %
         
    eval( sprintf ...
      (['varid = netcdf.inqVarID(ncid,', ...
        '''HDIF_EXPLICIT_3D_IRF_OVF_ENT_%02i'');'],ir) );
    hdif_d(:,:,:) = netcdf.getVar(ncid,varid);
    hdif_p =double(permute ( hdif_d(:,:,:),[2,1,3])); % change the order
                                          % (i,j,k)to (j,i,k) and store in
                                          % double precision
    hdif_p = hdif_p .* (hdif_p < 1.0e36); % eliminate the nodata indicator
    hdif = zeros(ny,nx,nzz);
    hdif(:,:,1:nz) = hdif_p(:,:,:);       % resize with added z layer
                                     
    % Use the IRF pulse locations to locate the piece of the stencil
    %  to associate the data with the pulse %%%%%%%%%%%%%%

    col = e_irf_mask*(IRF(:).*index_m(:)); % index number of impulse for 
                                           % entrainment
    inz = find(col); % find the nonzero values
    % append the new values and corresponding row and column
    % indices to the list
    icol_e  = [icol_e;  col(inz(:)) ];     % impulse location col index
    irow_e  = [irow_e;  index_m(inz(:)) ]; % response location  index
    a_val_e = [a_val_e; a_var(inz(:)) ];   % response value for adv
    h_val_e = [h_val_e; hdif(inz(:)) ];    % response for hdif
     
    % ckoout IRF pulse locations  %%
    if collect_pulse          
      % Locate each pulse, record pulse; record frequency of each related cell having
      %  nz value; and check that the col_index is the same as for the
      %  pulse.
      
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
      end % for z_IRF = 1:nz
    end % if collect_pulse
    % end checkout IRF  %%
      
  end  % for ntracers

  return
end % function get_e_data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  30  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


