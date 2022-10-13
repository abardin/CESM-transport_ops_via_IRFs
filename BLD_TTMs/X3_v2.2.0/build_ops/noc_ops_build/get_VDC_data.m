function [KTF, KBF] = get_VDC_data(MET, FILE_PARAM, GEO_PARAM)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Retrieves the data for the Vertical Diffusion Coefficients from 
% the history file.
% The separately recorded VDC_KPP and the correction VDC_GM are added
% together to form KTF, kappa vertical diffusion coefficient.
% The data are in a full matrix.
%
% INPUTS:
%   GEO_PARAM:
%      ny, nx, nzz           % 3-D dimensions for the full grid size
%   FILE_PARAM:
%      ncid                  % key identifying the opened history file
%   MET                       
%      MASK                  % land/sea mask
% 
% OUTPUTS:
%   KTF                      % Kappa at the top of the grid cell, m^2/s. 
%                            % 3-D array
%   KBF                      % kappa at at the bottom of the grid cell m^2/s
%
%   Note: The additional level of land has been added to the z-dimension
%            of the parent POP model.
%   Note: Output is in YXZ order (Z varies slowest) in the Matlab
%            sparse matrix transport operators.  In CESM everything is in
%            XYZ.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ncid = FILE_PARAM.ncid;
  ny = GEO_PARAM.ny; nx = GEO_PARAM.nx; nzz = GEO_PARAM.nzz; nz = nzz-1;
  KMT = MET.KMT;

  % extract the VDC_KPP data %%
  eval(sprintf ...
      ('varid = netcdf.inqVarID(ncid,''VDC_S'');') );
  VDC_S_d(:,:,:) = netcdf.getVar(ncid,varid);
   
  VDC_S_p =double(permute ( VDC_S_d(:,:,:),[2,1,3]));  % change the order
                                      % (i,j,k)to (j,i,k) and store in
                                      % double precision
  % VDC_S_screen = (VDC_S_p < 1.0e36);  % eliminate the no data indicator
  VDC_S_p = VDC_S_p .* (VDC_S_p < 1.0e36);      % Just the data, please
  VDC_S_p = (VDC_S_p ./ (100^2));       % from (cm^2)(s^-1) to (m^2)(s^-1) 
  VDC_S   = zeros(ny,nx,nzz);
  VDC_S(:,:,1:nz) = VDC_S_p(:,:,:);         % resize with added z layer

  % extract the GM data %%

  cmd = sprintf ...
      ('varid = netcdf.inqVarID(ncid,''VDC_GM'');');
  eval(cmd);
  GM_d(:,:,:) = netcdf.getVar(ncid,varid);
  GM_p =double(permute ( GM_d(:,:,:),[2,1,3]));  % change the order
                                      % (i,j,k)to (j,i,k) and store in
                                      % double precision
                                      % eliminate the nodata indicator
  GM_p = GM_p .* (GM_p < 1.0e36);     % just real data and zeros
  GM_p = (GM_p ./ (100^2));           % from (cm^2)(s^-1) to (m^2)(s^-1) 
  GM = zeros(ny,nx,nzz);
  GM(:,:,1:nz) = GM_p(:,:,:);         % resize with added z layer          
 
  % add the VDC_S and GM components
  KBF = VDC_S + GM;  % KBF is Kappa at the BOTTOM of the grid cell
  KBF = KBF .* MET.MASK;              % set sea floor cells to zero for sure
  % check that KBF at KMT == 0
  cell_cnt = 0;
  for j = 1:ny
      for i = 1:nx
          if KMT(j,i) ~= 0                  % if not land at surface
              if KBF(j,i,KMT(j,i)) ~= 0     % check that bottom flux = 0
                  disp( 'bottom cell flux KBF not = 0')
                  cell_cnt = cell_cnt + 1;
                  keyboard
              end % if bottom flux ~= 0
          end % if not land
      end % for i
  end % for j
  if cell_cnt == 0
      % disp('KBF bottom flux all zero')
  else
      disp (['KBF bottom flux not zero for ',int2str(cell_cnt), ...
             ' locations' ]);
  end 
              
  %% convert KBF to KTF
 
  % Calculate kappa at the TOP of the T grid box 
  KTF = zeros(ny,nx,nzz);
  KTF(:,:,2:nzz) = KBF(:,:,1:nzz-1);

  return
end % function get_VDC_data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  30  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


