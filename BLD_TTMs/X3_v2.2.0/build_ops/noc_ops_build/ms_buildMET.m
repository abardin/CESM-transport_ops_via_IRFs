function [MET] = ms_buildMET(F)
% returns a structure, MET, containing geometry for the configuration,f
% including the marginal seas 
% Inputs:  
%      source for the configuration: derived from
%      F.ref_dir
%      F.ref_file
%     
%  configuration dimensions: ny,nx,nz DEFINED INTERNALLY
%
%  marginal_mask is an (additional) mask for the marginal seas to
%      be excluded from the ocean configuration for the operators to be 
%      used.
%  Notes: 1) Dimensions output are in meters.
%         2) Variables are generally output in 3_D, even if this means they are
%            replicated. Exceptions are 
%              KMT, which is 2-D, and
%              iocn, which is a vector of 1-D indices into the ocean cells
%              of the 3-D grid geometry.
%         3) Refer to the POP Reference Manual, Chapter 3, for figures
%            and further text depicting the geometry definitions.
%            In particular, see figures 3.1, 3.2, and 3.3. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%  CONFIGURATION DEFINITIONS - to be modified for a different configuration:
  ny = 116; nx = 100; nz = 60; nzz = nz+1; % extra layer at the bottom

  MASK_MARGINAL = ones(ny,nx,nzz);       % include all
  %
  % read in the metric coefficients frpm the netcdf file
  % ncfile = sprintf('%s%04i-%02i.nc',F.ref_file_base,year,month);
  ncfile = ([F.ref_dir,F.ref_file]);
  ncid   = netcdf.open(ncfile,'NOWRITE');
  
  % use the KMT (index of the bottom wet cell (land at the surface = 0)
  % to create the land sea mask  
  varid = netcdf.inqVarID(ncid,'KMT'); % number of wet cells above seafloor
  KMT_d(:,:) = netcdf.getVar(ncid,varid);
  KMT = permute (KMT_d,[2,1]); % reorder to (ny,nx) 
  % mask for ocn points only, full mask per OGCM, includes marginal seas
  [ny,nx] = size(KMT); 
  nz = max(KMT(:)); nzz = nz+1;        % Offline variables have a seafloor;
                                       % k index increased by 1;
  MASK = zeros(ny,nx,nzz);
  for k = 1:nzz
    MASK(:,:,k) = (k <= KMT(:,:));
  end % all depths 
  
  % check expected dimensions
  [nyck,nxck,nzzck] = size(MASK);
  if (nyck ~=ny) ||(nxck ~=nx) || (nzzck ~=nzz)
      disp('dimensions not as expected in buildMET');
      keyboard
  end
  
  %  mask out the small and irrelevant closed basins
  MASK = MASK .* MASK_MARGINAL;
  % 
  iocn = find(MASK(:) > 0); % 1D index into 3D structure for ocean grid cells
  MASKsurf = MASK(:,:,1);
  isurfocn = find(MASKsurf(:) >0); % 1D index into surface ocean grid cells

  % get HTE, horizontal length of the East face of the T-cell   
  varid = netcdf.inqVarID(ncid,'HTE');
  HTE_d(:,:) = netcdf.getVar(ncid,varid,'double' );
  HTE_d = HTE_d ./ (100);      % convert from centimeters to m
  HTE = permute (HTE_d,[2,1]); % reorder to (ny,nx)
  HTE = HTE(:,:,ones(nzz,1));  % extend to 3_D, with added layer at the bottom
  %
  % get HTN, horizontal length of the North face of the T-cell 
  varid = netcdf.inqVarID(ncid,'HTN');
  HTN_d(:,:) = netcdf.getVar(ncid,varid,'double' );
  HTN_d = HTN_d ./ (100);      % convert from centimeters to m
  HTN = permute (HTN_d,[2,1]); % reorder to (ny,nx)
  HTN = HTN(:,:,ones(nzz,1));  % extend to 3_D, with added layer at the bottom
  %
  % get DXU, horizontal length of the North face of the U_cell, centered on
  % U-points
  varid = netcdf.inqVarID(ncid,'DXU');
  DXU_d(:,:) = netcdf.getVar(ncid,varid,'double' );
  DXU_d = DXU_d ./ (100);      % convert from centimeters to m
  DXU = permute (DXU_d,[2,1]); % reorder to (ny,nx)
  DXU = DXU(:,:,ones(nzz,1));  % extend to 3_D, with added layer at the bottom
  %  
  % get DXT, horizontal East-West length of the T-cell, centerered at
  % T_points (that is, mid-Tcell).
  varid = netcdf.inqVarID(ncid,'DXT');
  DXT_d(:,:) = netcdf.getVar(ncid,varid,'double' );
  DXT_d = DXT_d ./ (100);      % convert from centimeters to m
  DXT = permute (DXT_d,[2,1]); % reorder to (ny,nx)
  DXT = DXT(:,:,ones(nzz,1));  % extend to 3_D, with added layer at the bottom
  %  
  % get DYU, horizontal length of the East face of  the U-cell
  varid = netcdf.inqVarID(ncid,'DYU');
  DYU_d(:,:) = netcdf.getVar(ncid,varid,'double' );
  DYU_d = DYU_d ./ (100);      % convert from centimeters to m
  DYU = permute (DYU_d,[2,1]); % reorder to (ny,nx)
  DYU = DYU(:,:,ones(nzz,1));  % extend to 3_D, with added layer at the bottom
  %  
  % get DYT, horizontal North-South length of the T-cell, centered at 
  % T_points (that is, mid-Tcell).
  varid = netcdf.inqVarID(ncid,'DYT');
  DYT_d(:,:) = netcdf.getVar(ncid,varid,'double' );
  DYT_d = DYT_d ./ (100);      % convert from centimeters to m
  DYT = permute (DYT_d,[2,1]); % reorder to (ny,nx)
  DYT = DYT(:,:,ones(nzz,1));  % extend to 3_D, with added layer at the bottom 
  %  
  % get DZT, vertical thickness of the T_cell (same in all horizontal
  % directions)
  % This variable is the same as dz in the POP reference manual.
  varid = netcdf.inqVarID(ncid,'dz');
  DZT(1,1,:) = netcdf.getVar(ncid,varid,'double' ) ./100; % convert cm to m
  DZT = DZT(1,1,[1:nzz-1]);    % extend to 3_D
  DZT(1,1,nzz) = 200;          % small non-zero depth for added layer
  DZT = DZT(ones(ny,1),ones(nx,1),:); 
  %               
  % get DZU, vertical thickness of the U_cell 
  % (Note that DZU(k=1) = (1/2)*DZT(k=1)
  % This variable is the same as dzw in the POP reference manual.
  varid = netcdf.inqVarID(ncid,'dzw');
  DZU(1,1,:) = netcdf.getVar(ncid,varid,'double' ) ./100; % convert cm to m
  DZU = DZU(1,1,[1:nzz-1]);    % extend to 3_D
  DZU(1,1,nzz) = 200;          % small non-zero depth for added layer
  DZU = DZU(ones(ny,1),ones(nx,1),:);
  DZU(:,:,1) = DZU(:,:,1);
  %  
  % get ZT, depth of T_points (positive, from the surface, mid-T_cell)
  varid = netcdf.inqVarID(ncid,'z_t');
  ZT(1,1,:) = netcdf.getVar(ncid,varid,'double' ) ./100; % convert cm to m
  ZT = ZT(1,1,[1:nzz-1]);    % extend to 3_D
  ZT(1,1,nzz) = ZT(1,1,nzz-1) + DZU(1,1,nzz);    % small non-zero depth for added layer
  ZT = ZT(ones(ny,1),ones(nx,1),:);
  %
  % get ZU, depth of the top of T_cells (positive, from the surface)
  % This variable is the same as z_w in the POP reference manual.
  varid = netcdf.inqVarID(ncid,'z_w');
  ZU(1,1,:) = netcdf.getVar(ncid,varid,'double' ) ./100; % convert cm to m
  ZU = ZU(1,1,[1:nzz-1]);    % extend to 3_D
  ZU(1,1,nzz) = ZU(1,1,nzz-1) + DZT(1,1,nzz);  % small non-zero depth for added layer
  ZU = ZU(ones(ny,1),ones(nx,1),:);
  %  
  % lat and long, at the middle of the T_cell
  varid = netcdf.inqVarID(ncid,'TLAT');
  TLAT_d(:,:) = netcdf.getVar(ncid,varid,'double' );
  TLAT = permute (TLAT_d,[2,1]);  % reorder to (ny,nx)
  TLAT = TLAT(:,:,ones(nzz,1));   % extend to 3_D

  varid = netcdf.inqVarID(ncid,'TLONG');
  TLONG_d(:,:) = netcdf.getVar(ncid,varid,'double' );
  TLONG = permute (TLONG_d,[2,1]);  % reorder to (ny,nx)
  TLONG = TLONG(:,:,ones(nzz,1));
  %  
  % TAREA, area of the top of a T_cell
  varid = netcdf.inqVarID(ncid,'TAREA');
  TAREA_d(:,:) = netcdf.getVar(ncid,varid,'double' );
  TAREA = permute (TAREA_d,[2,1]);  % reorder to (ny,nx)
  TAREA = TAREA ./ (100^2);         % convert from sq centimeters to sq m
  TAREA = TAREA(:,:,ones(nzz,1));   % extend to 3_D
  %
  % ANGLE, angle grid makes with latitude line
  varid = netcdf.inqVarID(ncid,'ANGLE');
  ANGLE_d(:,:) = netcdf.getVar(ncid,varid,'double' );
  ANGLE = permute (TAREA_d,[2,1]);  % reorder to (ny,nx)
  ANGLE = ANGLE(:,:,ones(nzz,1));   % extend to 3_D
  %   
  % volume of the T_cell (m^3)
  VOL = TAREA.*DZT;
  %  
  % the area of the top, east and north faces of the T_cell (m^2)
  TFACE = TAREA;
  EFACE = HTE.*DZT;
  NFACE = HTN.*DZT;
  %  
  % earth radius
  varid = netcdf.inqVarID(ncid,'radius');
  a = netcdf.getVar(ncid,varid,'double' ); 
  a = a/100;                        % convert from centimeters to m
  %
  netcdf.close(ncid);
  
  % create the structure 
  MET.MASK = MASK;
  MET.HTE  = HTE;
  MET.HTN  = HTN;
  MET.DZT  = DZT;
  MET.DXU  = DXU;
  MET.DYU  = DYU;
  MET.DXT  = DXT;
  MET.DYT  = DYT;
  MET.DZU  = DZU;
  MET.TAREA = TAREA;
  MET.ANGLE = ANGLE;
  MET.VOL   = VOL;
  MET.TFACE = TFACE;
  MET.EFACE = EFACE;
  MET.NFACE = NFACE;
  MET.TLAT  = TLAT;
  MET.TLONG = TLONG;
  MET.ZT = ZT;
  MET.ZU = ZU;
  MET.a  = a;
  MET.iocn = iocn;
  MET.KMT  = KMT;  
  MET.grain = 1;                   % full sized MET
  return 
end % function ms_ buildMET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  30 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




