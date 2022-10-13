function Build_state_data(F, year, period)
% Gathers basic ocean state variables from the run .nc files, and puts them
% into .mat files per month in the F.noc_dir.
% Monthly average values are generated, available to be combined into 
% climatological averages..
% Data in each file:
%  SSH, TEMP (potential temprature), PD (potential density)
%  SALT (salinity), RHO (in-situ density), IFRAC (sea-ice fraction),
% and many more 
% The data are in dimensions of ny X nx X nz for 3_D data, and ny X nx for
% 2-D data.
%% Added 5/19 extra fields for oxygen
%% modified for v2.1.3, with variable surface height thickness
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F.ref_dir        % archive directory for CESM run
% F.ref_file_base  % first part of filename in CESM ocn/hist source
% F.zeroyr         % in CESM ocn/hist source terms
% F.firstyr        % in CESM ocn/hist source terms
% F.nyears         % number of years to process
% F.noc_dir        % output directory for state data files
                   % MET.mat is assumed present in F.noc_dir

  seq_num = (year - F.firstyr)*F.num_periods + period; % for output file
  out_file_base = 'model_state_data_';

  % set up dimensions from MET
  if ~exist([F.noc_dir,'MET.mat'],'file')
    disp(['no MET file']);
    keyboard
  else
     eval (['load ',F.noc_dir,'MET.mat']);
  end
 
  [ny,nx,nz] = size(MET.MASK);
  iocn = MET.iocn;
  surf_mask = MET.MASK(:,:,1);
  isurfocn = find(surf_mask ~= 0);

   % set up the filename for this period
    if F.day5
       [cal_mo,cal_day] = get_index_5day(period);
       data_file = sprintf('%s%04i-%02i-%02i.nc',F.ref_file_base,year,cal_mo,cal_day);
    else
       data_file = sprintf('%s%04i-%02i.nc',F.ref_file_base,F.file_yr,period); 
    end
      filename = [F.ref_dir,data_file];
      ncid = netcdf.open(filename,'NC_NOWRITE'); 
      disp ('getting monthly data from '); disp(filename);

      % get the SSH monthly average  
      varid = netcdf.inqVarID(ncid,'SSH');
      SSHd(:,:) = netcdf.getVar(ncid,varid,'double' );
      SSHd = SSHd ./ (100);            % convert from centimeters to meters
      SSHp = permute (SSHd,[2,1]);     % reorder to (ny,nx)
      SSH  = zeros(ny,nx);
      SSH(isurfocn) = SSHp(isurfocn);  % eliminate marginal ocean 
                                       % and land grid-cells
      maxSSH = max(SSH(isurfocn)); 
      minSSH = min(SSH(isurfocn));     % sanity check
      if (maxSSH > 5) || (minSSH < -5 )
          disp ('SSH out of range')
          keyboard
      end
      
      % get the temperature monthly average  
      varid = netcdf.inqVarID(ncid,'TEMP');  % potential temp., degrees C
      TEMP_d(:,:,:) = netcdf.getVar(ncid,varid,'double'); 
      TEMPp = permute(TEMP_d(:,:,:),[2,1,3]);  %reorder to (ny,nx)
      TEMP = zeros(ny,nx,nz);
      TEMP(iocn) = TEMPp(iocn);
      maxTEMP = max(TEMP(iocn));
      minTEMP = min(TEMP(iocn));
      if (maxTEMP > 40 || minTEMP < -10)
          disp('TEMP out of range')
          keyboard
      end
      
      %% get the salinity monthly averate
      varid = netcdf.inqVarID(ncid,'SALT'); % gram/kilogram
      SALT_d(:,:,:) = netcdf.getVar(ncid,varid,'double'); 
      SALTp = permute(SALT_d(:,:,:),[2,1,3]);  %reorder to (ny,nx)
      SALT = zeros(ny,nx,nz);
      SALT(iocn) = SALTp(iocn);
      maxSALT = max(SALT(iocn));
      minSALT = min(SALT(iocn));
      if (maxSALT > 0.1 || minSALT < 0)
          disp('SALT out of range')
          keyboard
      end
      
      %% get the density monthly average
      varid = netcdf.inqVarID(ncid,'RHO'); % gram/cm^3
      RHO_d(:,:,:) = netcdf.getVar(ncid,varid,'double'); 
      RHOp = permute(RHO_d(:,:,:),[2,1,3]);  %reorder to (ny,nx)
      RHOp = RHOp.*(100^(3) / 1000); % convert units to kgm/m^3
      RHO = zeros(ny,nx,nz);
      RHO(iocn) = RHOp(iocn);
      maxRHO = max(RHO(iocn));
      minRHO = min(RHO(iocn));
      if (maxRHO > 1100|| minRHO < 1000)
          disp('RHO out of range')
          disp (['minRHO = ',num2str(minRHO,4), ...
                ', maxRHO = ',num2str(maxRHO,4) ]);
      end
      
      %% get the potential density monthly average
      varid = netcdf.inqVarID(ncid,'PD'); % gram/cm^3
      PD_d(:,:,:) = netcdf.getVar(ncid,varid,'double'); 
      PDp = permute(PD_d(:,:,:),[2,1,3]);  %reorder to (ny,nx)
      PDp = PDp.*(100^(3) / 1000); % convert units to kgm/m^3
      PD = zeros(ny,nx,nz);
      PD(iocn) = PDp(iocn);
      maxPD = max(PD(iocn));
      minPD = min(PD(iocn));
      if (maxPD > 1100|| minPD < 1000)
          disp('PD out of range')
          disp (['minPD = ',num2str(minPD,4), ...
                ', maxPD = ',num2str(maxPD,4) ]);
      end
      
      %% get the sea-ice fraction monthly average  
      varid = netcdf.inqVarID(ncid,'IFRAC');
      IFRACd(:,:) = netcdf.getVar(ncid,varid,'double' );
      IFRACp = permute (IFRACd,[2,1]);     % reorder to (ny,nx)
      IFRAC  = zeros(ny,nx);
      IFRAC(isurfocn) = IFRACp(isurfocn);  % eliminate marginal ocean 
                                           % and land grid-cells
      maxIFRAC = max(IFRAC(isurfocn)); 
      minIFRAC = min(IFRAC(isurfocn));             % sanity check
      if (maxIFRAC > 1) || (minIFRAC < 0 )
          disp ('IFRAC out of range')
          keyboard
      end
  
%%  % NO CFC data 

 %% get the CFC 11 concentration monthly average fmol/m^3
%      varid = netcdf.inqVarID(ncid,'CFC11'); % 
%      CFC11d(:,:,:) = netcdf.getVar(ncid,varid,'double'); 
%      CFC11p = permute(CFC11d(:,:,:),[2,1,3]);  %reorder to (ny,nx)
%      CFC11 = zeros(ny,nx,nz);
%      CFC11(iocn) = CFC11p(iocn);
%      maxCFC11 = max(CFC11(iocn));
%      minCFC11 = min(CFC11(iocn));
%      disp(['max CFC11 = ',num2str(maxCFC11,4), ...
%            ', min CFC11 = ',num2str(minCFC11,4)]);
     % if (maxO2 > 1e3 || minO2 < 0)
     %	     disp('O2 out of range');  
     % end

     %% get the CFC_sea-ice fraction monthly average  
%       varid = netcdf.inqVarID(ncid,'CFC_IFRAC');
%       CFC_IFRACd(:,:) = netcdf.getVar(ncid,varid,'double' );
%       CFC_IFRACp = permute (CFC_IFRACd,[2,1]);     % reorder to (ny,nx)
%       CFC_IFRAC  = zeros(ny,nx);
%       CFC_IFRAC(isurfocn) =CFC_IFRACp(isurfocn);  % eliminate marginal ocean 
%                                                   % and land grid-cells
%       maxCIFRAC = max(CFC_IFRAC(isurfocn)); 
%       minCIFRAC = min(CFC_IFRAC(isurfocn));             % sanity check
%       if (maxCIFRAC > 1) || (minCIFRAC < 0 )
%           disp ('CFC_IFRAC out of range')
%           keyboard
%       end
%
  %% get the CFC XKW monthly average  
%       varid = netcdf.inqVarID(ncid,'CFC_XKW');
%       CFC_XKWd(:,:) = netcdf.getVar(ncid,varid,'double' );
%       CFC_XKWp = permute (CFC_XKWd,[2,1]);     % reorder to (ny,nx)
%       CFC_XKWp = CFC_XKWp .* 0.01;             % convert cm/s to m/s
%       CFC_XKW  = zeros(ny,nx);
%       CFC_XKW(isurfocn) = CFC_XKWp(isurfocn);  % eliminate marginal ocean 
%                                                % and land grid-cells
%       maxXKW = max(CFC_XKW(isurfocn)); 
%       minXKW = min(CFC_XKW(isurfocn));             % sanity check
%       disp(['max XKW = ',num2str(maxXKW,4),', min XKW = ',num2str(minXKW,4)]);
%
   %% get the CFC atmospheric pressure monthly average (atm)  
%       varid = netcdf.inqVarID(ncid,'CFC_ATM_PRESS');
%       CFC_APd(:,:) = netcdf.getVar(ncid,varid,'double' );
%       CFC_APp = permute (CFC_APd,[2,1]);     % reorder to (ny,nx)
%       CFC_ATM_PRESS  = zeros(ny,nx);
%       CFC_ATM_PRESS(isurfocn) = CFC_APp(isurfocn);  % eliminate marginal ocean% 
%                                                     % and land grid-cells
%       maxAP = max(CFC_ATM_PRESS(isurfocn)); 
%       minAP = min(CFC_ATM_PRESS(isurfocn));         % sanity check
%       disp(['max ATM PRESS = ',num2str(maxAP,4), ...
%             ', min ATM PRESS = ',num2str(minAP,4)]);
%
%   %% get the CFC 11 surface flux monthly average (fmol/m^2/s) 
%       varid = netcdf.inqVarID(ncid,'STF_CFC11');
%       CFC11_PVd(:,:) = netcdf.getVar(ncid,varid,'double' );
%       CFC11_PVp = permute(CFC11_PVd,[2,1]);      % reorder to (ny,nx)
%       CFC11_PVp = CFC11_PVp .* (0.01)^2;         % convert cm^2 to m
%       STF_CFC11 = zeros(ny,nx);
%       STF_CFC11(isurfocn) = CFC11_PVp(isurfocn); % eliminate marginal ocean 
%                                                  % and land grid-cells
%       maxSTF_CFC11 = max(STF_CFC11(isurfocn)); 
%       minSTF_CFC11 = min(STF_CFC11(isurfocn));     % sanity check
%       disp(['max STF CFC11 = ',num2str(maxSTF_CFC11,4), ...
%             ', min STF CFC11 = ',num2str(minSTF_CFC11,4)]);

   %% get the ECOSYS atmospheric pressure monthly average (atm)  
      varid = netcdf.inqVarID(ncid,'ECOSYS_ATM_PRESS');
      ECO_APd(:,:) = netcdf.getVar(ncid,varid,'double' );
      ECO_APp = permute (ECO_APd,[2,1]);     % reorder to (ny,nx)
      ECOSYS_ATM_PRESS  = zeros(ny,nx);
      ECOSYS_ATM_PRESS(isurfocn) = ECO_APp(isurfocn);  % eliminate marginal ocean 
                                                    % and land grid-cells
      maxAP = max(ECOSYS_ATM_PRESS(isurfocn)); 
      minAP = min(ECOSYS_ATM_PRESS(isurfocn));         % sanity check
      disp(['max ECO ATM PRESS = ',num2str(maxAP,4), ...
            ', min ECO ATM PRESS = ',num2str(minAP,4)]);
      % if (maxIFRAC > 1) || (minIFRAC < 0 )
      %     disp ('IFRAC out of range')
      %     keyboard
      % end

   %% get the ECOSYS ice fraction monthly average (atm)  
      varid = netcdf.inqVarID(ncid,'ECOSYS_IFRAC');
      ECO_IFRd(:,:) = netcdf.getVar(ncid,varid,'double' );
      ECO_IFRp = permute (ECO_IFRd,[2,1]);     % reorder to (ny,nx)
      ECOSYS_IFRAC  = zeros(ny,nx);
      ECOSYS_IFRAC(isurfocn) = ECO_IFRp(isurfocn);  % eliminate marginal ocean 
                                                    % and land grid-cells
      % maxIF = max(ECOSYS_IFRAC(isurfocn)); 
      % minIF = min(ECOSYS_IFRAC(isurfocn));         % sanity check
      % disp(['max ECO IFRAC = ',num2str(maxIF,4), ...
      %      ', min ECO IFRAC = ',num2str(minIF,4)]);
      % if (maxIFRAC > 1) || (minIFRAC < 0 )
      %     disp ('IFRAC out of range')
      %     keyboard
      % end

 %% get the ECOSYS XKW monthly average (atm)  
      varid = netcdf.inqVarID(ncid,'ECOSYS_XKW');
      ECO_XKWd(:,:) = netcdf.getVar(ncid,varid,'double' );
      ECO_XKWp = permute (ECO_XKWd,[2,1]);     % reorder to (ny,nx)
      ECO_XKWp = ECO_XKWp ./100; % convert cm/s to m/s
      ECOSYS_XKW  = zeros(ny,nx);
      ECOSYS_XKW(isurfocn) = ECO_XKWp(isurfocn);  % eliminate marginal ocean 
                                                    % and land grid-cells
     % maxIF = max(ECOSYS_XKW(isurfocn)); 
     % minIF = min(ECOSYS_XKW(isurfocn));         % sanity check
     % disp(['max ECO XKW = ',num2str(maxIF,4), ...
      %      ', min ECO XKW = ',num2str(minIF,4)]);
      % if (maxIFRAC > 1) || (minIFRAC < 0 )
      %     disp ('IFRAC out of range')
      %     keyboard
      % end
 %% get the IAGE monthly average, years
      varid = netcdf.inqVarID(ncid,'IAGE'); % 
      IAGEd(:,:,:) = netcdf.getVar(ncid,varid,'double'); 
      IAGEp = permute(IAGEd(:,:,:),[2,1,3]);  %reorder to (ny,nx)
      IAGE = zeros(ny,nx,nz);
      IAGE(iocn) = IAGEp(iocn);
      maxIAGE = max(IAGE(iocn));
      minIAGE = min(IAGE(iocn));
      disp(['max IAGE = ',num2str(maxIAGE,4), ...
            ', min O2 = ',num2str(minIAGE,4)]);
      if (maxIAGE > 1e6 || minIAGE < -1e3)
	     disp('IAGE out of range');  
      end

 %% get the O2 concentration monthly average mmol O2 /m^3
      varid = netcdf.inqVarID(ncid,'O2'); % 
      O2d(:,:,:) = netcdf.getVar(ncid,varid,'double'); 
      O2p = permute(O2d(:,:,:),[2,1,3]);  %reorder to (ny,nx)
      O2 = zeros(ny,nx,nz);
      O2(iocn) = O2p(iocn);
      maxO2 = max(O2(iocn));
      minO2 = min(O2(iocn));
      disp(['max O2 = ',num2str(maxO2,4), ...
            ', min O2 = ',num2str(minO2,4)]);
      if (maxO2 > 1e3 || minO2 < 0)
	     disp('O2 out of range');  
      end

 %% get O2 production rate mmol /m^3 /s
      varid = netcdf.inqVarID(ncid,'O2_PRODUCTION'); % mmol /m^3 /s
      O2_PRODd(:,:,:) = netcdf.getVar(ncid,varid,'double');
      [pnx,pny,pnz] = size(O2_PRODd);
      O2_PRODUCTION = zeros(ny,nx,nz);
      O2_PRODUCTION(:,:,1:pnz) = double(permute(O2_PRODd,[2 1 3]) );
      O2_PROD_iocn = O2_PRODUCTION(iocn);
      O2_PRODUCTION = zeros(ny,nx,nz);
      O2_PRODUCTION(iocn) = O2_PROD_iocn;     
      maxO2P = max(O2_PRODUCTION(iocn));
      minO2P = min(O2_PRODUCTION(iocn));
      disp(['max O2 production = ',num2str(maxO2P,4), ...
            ', min O2 production = ',num2str(minO2P,4)]);
      if (maxO2P > 1e-2|| minO2P < 0)
          disp('O2 PRODUCTION out of range')   
      end

 %% get O2 consumption rate mmol /m^3 /s
      varid = netcdf.inqVarID(ncid,'O2_CONSUMPTION'); % mmol /m^3 /s
      O2_CONSd(:,:,:) = netcdf.getVar(ncid,varid,'double');
      [pnx,pny,pnz] = size(O2_CONSd);
      O2_CONSUMPTION = zeros(ny,nx,nz);
      O2_CONSUMPTION(:,:,1:pnz) = double(permute(O2_CONSd,[2 1 3]) );
      O2_CONS_iocn = O2_CONSUMPTION(iocn);
      O2_CONSUMPTION = zeros(ny,nx,nz);
      O2_CONSUMPTION(iocn) = O2_CONS_iocn;     
      maxO2C = max(O2_CONSUMPTION(iocn));
      minO2C = min(O2_CONSUMPTION(iocn));
      disp(['max O2 consumtion = ',num2str(maxO2C,4), ...
            ', min O2 consumption = ',num2str(minO2C,4)]);
      if (maxO2C > 1e-3|| minO2C < 0)
          disp('O2 consumption out of range')   
      end

 %% get the O2 Schmidt number monthly average (atm)  
      varid = netcdf.inqVarID(ncid,'SCHMIDT_O2');
      O2_Schd(:,:) = netcdf.getVar(ncid,varid,'double' );
      O2_Schp = permute (O2_Schd,[2,1]);     % reorder to (ny,nx)
      SCHMIDT_O2  = zeros(ny,nx);
      SCHMIDT_O2(isurfocn) = O2_Schp(isurfocn);  % eliminate marginal ocean 
                                                    % and land grid-cells
     % maxIF = max(SCHMIDT_O2(isurfocn)); 
      %minIF = min(SCHMIDT_O2(isurfocn));         % sanity check
      %disp(['max SCHMIDT O2 = ',num2str(maxIF,4), ...
      %      ', min SCHMIDT O2 = ',num2str(minIF,4)]);
      % if (maxIFRAC > 1) || (minIFRAC < 0 )
      %     disp ('IFRAC out of range')
      %     keyboard
      % end

 %% get the O2 Saturation monthly average (atm)  
      varid = netcdf.inqVarID(ncid,'O2SAT');
      O2SATd(:,:) = netcdf.getVar(ncid,varid,'double' );
      O2SATp = permute (O2SATd,[2,1]);     % reorder to (ny,nx)
      O2SAT  = zeros(ny,nx);
      O2SAT(isurfocn) = O2SATp(isurfocn);  % eliminate marginal ocean 
                                                    % and land grid-cells
      %maxIF = max(O2SAT(isurfocn)); 
      %minIF = min(O2SAT(isurfocn));         % sanity check
      %disp(['max O2SAT = ',num2str(maxIF,4), ...
       %     ', min O2SAT = ',num2str(minIF,4)]);
      % if (maxIFRAC > 1) || (minIFRAC < 0 )
      %     disp ('IFRAC out of range')
      %     keyboard
      % end

 %% get DOP
      varid = netcdf.inqVarID(ncid,'DOP'); %   mmol/m^3
      DOPd(:,:,:) = netcdf.getVar(ncid,varid,'double');
      DOPp = double(permute(DOPd,[2 1 3]) );
      DOP = zeros(ny,nx,nz);
      DOP(iocn) = DOPp(iocn);     
      maxDOP = max(DOP(iocn));
      minDOP = min(DOP(iocn));
      disp(['max DOP  = ',num2str(maxDOP,4), ...
            ', min DOP = ',num2str(minDOP,4)]);
     % if (maxPOC_remin > 1e-2|| minPOC_remin < 0)
     %     disp('POC remin out of range')   
     % end

 %% get PO4
      varid = netcdf.inqVarID(ncid,'PO4'); %   mmol/m^3
      PO4d(:,:,:) = netcdf.getVar(ncid,varid,'double');
      PO4p = double(permute(PO4d,[2 1 3]) );
      PO4 = zeros(ny,nx,nz);
      PO4(iocn) = PO4p(iocn);     
      maxPO4 = max(PO4(iocn));
      minPO4 = min(PO4(iocn));
      disp(['max PO4  = ',num2str(maxPO4,4), ...
            ', min PO4 = ',num2str(minPO4,4)]);
     % if (maxPOC_remin > 1e-2|| minPOC_remin < 0)
     %     disp('POC remin out of range')   
     % end

      varid = netcdf.inqVarID(ncid,'NITRIF'); %  N mmol/m^3/s
      NITRIFd(:,:,:) = netcdf.getVar(ncid,varid,'double');
      NITRIFp = double(permute(NITRIFd,[2 1 3]) );
      NITRIF = zeros(ny,nx,nz);
      NITRIF(iocn) = NITRIFp(iocn);     
      maxNITRIF = max(NITRIF(iocn));
      minNITRIF = min(NITRIF(iocn));
      disp(['max NITRIF  = ',num2str(maxNITRIF,4), ...
            ', min NITRIF = ',num2str(minNITRIF,4)]);
     % if (maxNITRIF > 1e-2|| minNITRIF < 0)
     %     disp('NITRIF out of range')   
     % end


 %% get SedDenitrif  % nmolN/cm^2/s 2D => N mmol/m^3/s
      % N mmol/m^3/s = N nmol/cm^2/s * [1e-6 mmol/nmol] * [1e4 cm^2/m^2] ./ [dz(KMT) m];
      varid = netcdf.inqVarID(ncid,'SedDenitrif'); %
      SedDenitrifd(:,:) = netcdf.getVar(ncid,varid,'double');
      SedDenitrifp = double(permute(SedDenitrifd,[2 1]) );
      SedDenitrif = zeros(ny,nx,nz);
      % convert to 3D N mmol/m^3/s
      for k = 1:nz-1
	      mapKMT = find(MET.KMT == k);
          SedDeN_lev = zeros(ny,nx);
          dzlev = MET.DZT(:,:,k);
          SedDeN_lev(mapKMT) = (SedDenitrifp(mapKMT)) .* 1e-2 ./ dzlev(mapKMT);
          SedDenitrif(:,:,k) = SedDeN_lev;
      end % for each level
      maxSedDenitrif = max(SedDenitrif(iocn));
      minSedDenitrif = min(SedDenitrif(iocn));
      disp(['max SedDenitrif  = ',num2str(maxSedDenitrif,4), ...
            ', min SedDenitrif = ',num2str(minSedDenitrif,4)]);
     % if (maxSedDenitrif > 1e-2|| minSedDenitrif < 0)
     %     disp('SedDenitrif out of range')   
     % end 

 %% get alkalinity 
      varid = netcdf.inqVarID(ncid,'ALK'); %
      ALKd(:,:,:) = netcdf.getVar(ncid,varid,'double');
      ALKp = double(permute(ALKd,[2 1 3]) );
      ALK = zeros(ny,nx,nz);
      ALK(iocn) = ALKp(iocn);     
      maxALK = max(ALK(iocn));
      minALK = min(ALK(iocn));
      disp(['max ALK  = ',num2str(maxALK,4), ...
            ', min ALK = ',num2str(minALK,4)]);
      if (maxALK > 3e3|| minALK < 0)
          disp('ALK out of range')   
      end

 %% get the PH in the surface monthly average (atm)  
      varid = netcdf.inqVarID(ncid,'PH');
      PHd(:,:) = netcdf.getVar(ncid,varid,'double' );
      PHp = permute (PHd,[2,1]);     % reorder to (ny,nx)
      PH  = zeros(ny,nx);
      PH(isurfocn) = PHp(isurfocn);  % eliminate marginal ocean 
                                                    % and land grid-cells
      maxIF = max(PH(isurfocn)); 
      minIF = min(PH(isurfocn));         % sanity check
      disp(['max PH (surface) = ',num2str(maxIF,4), ...
            ', min PH = ',num2str(minIF,4)]);
      % if (maxIFRAC > 1) || (minIFRAC < 0 )
      %     disp ('IFRAC out of range')
      %     keyboard
      % end

 %% get the Schmidt number for CO2 monthly average (atm)  
      varid = netcdf.inqVarID(ncid,'SCHMIDT_CO2');
      CO2_Schd(:,:) = netcdf.getVar(ncid,varid,'double' );
      CO2_Schp = permute (CO2_Schd,[2,1]);     % reorder to (ny,nx)
      SCHMIDT_CO2  = zeros(ny,nx);
      SCHMIDT_CO2(isurfocn) = CO2_Schp(isurfocn);  % eliminate marginal ocean 
                                                    % and land grid-cells
      %maxIF = max(SCHMIDT_CO2(isurfocn)); 
      %minIF = min(SCHMIDT_CO2(isurfocn));         % sanity check
      %disp(['max Schmidt CO2 = ',num2str(maxIF,4), ...
       %     ', min Schmidt CO2 = ',num2str(minIF,4)]);
      % if (maxIFRAC > 1) || (minIFRAC < 0 )
      %     disp ('IFRAC out of range')
      %     keyboard
      % end

 %% get the ATM CO2 monthly average (atm)  
      varid = netcdf.inqVarID(ncid,'ATM_CO2');
      ATM_CO2d(:,:) = netcdf.getVar(ncid,varid,'double' );
      ATM_CO2p = permute (ATM_CO2d,[2,1]);     % reorder to (ny,nx)
      ATM_CO2  = zeros(ny,nx);
      ATM_CO2(isurfocn) = ATM_CO2p(isurfocn);  % eliminate marginal ocean 
                                                    % and land grid-cells
      maxIF = max(ATM_CO2(isurfocn)); 
      minIF = min(ATM_CO2(isurfocn));         % sanity check
      disp(['max ATM CO2 = ',num2str(maxIF,4), ...
            ', min ATM CO2 = ',num2str(minIF,4)]);
      % if (maxIFRAC > 1) || (minIFRAC < 0 )
      %     disp ('IFRAC out of range')
      %     keyboard
      % end

 %% get the transition layer thickness monthly average (atm)  
      varid = netcdf.inqVarID(ncid,'TLT');
      TLTd(:,:) = netcdf.getVar(ncid,varid,'double' );
      TLTp = permute (TLTd,[2,1]);     % reorder to (ny,nx)
      TLTp = TLTp ./100; % convert cm to m
      TLT  = zeros(ny,nx);
      TLT(isurfocn) = TLTp(isurfocn);  % eliminate marginal ocean 
                                                    % and land grid-cells
      maxIF = max(TLT(isurfocn)); 
      minIF = min(TLT(isurfocn));         % sanity check
      disp(['max TLT transition layer thickness = ',num2str(maxIF,4), ...
            ', min TLT = ',num2str(minIF,4)]);
      % if (maxIFRAC > 1) || (minIFRAC < 0 )
      %     disp ('IFRAC out of range')
      %     keyboard
      % end

 %% get the interior depth level monthly average   
      varid = netcdf.inqVarID(ncid,'INT_DEPTH');
      INT_DEPTHd(:,:) = netcdf.getVar(ncid,varid,'double' );
      INT_DEPTHp = permute (INT_DEPTHd,[2,1]);     % reorder to (ny,nx)
      INT_DEPTHp = INT_DEPTHp ./100; % convert cm to m
      INT_DEPTH  = zeros(ny,nx);
      INT_DEPTH(isurfocn) = INT_DEPTHp(isurfocn);  % eliminate marginal ocean 
                                                    % and land grid-cells
     % maxIF = max(INT_DEPTH(isurfocn)); 
     % minIF = min(INT_DEPTH(isurfocn));         % sanity check
      %disp(['max INT DEPTH = ',num2str(maxIF,4), ...
       %     ', min INT DEPTH = ',num2str(minIF,4)]);
      % if (maxIFRAC > 1) || (minIFRAC < 0 )
      %     disp ('IFRAC out of range')
      %     keyboard
      % end
 
   %% get the average velocity in the u-direction
      varid = netcdf.inqVarID(ncid,'UVEL'); % cm/s
      UVELd(:,:,:) = netcdf.getVar(ncid,varid,'double'); 
      UVELp = permute(UVELd(:,:,:),[2,1,3]);  %reorder to (ny,nx)
      UVELp = UVELp.*.01; % convert units to m/s
      UVEL = zeros(ny,nx,nz);
      UVEL(iocn) = UVELp(iocn);
     % maxUVEL = max(UVEL(iocn));
     % minUVEL = min(UVEL(iocn));
     % disp(['max UVEL = ',num2str(maxUVEL,4), ...
      %      ', min UVEL = ',num2str(minUVEL,4)]);
     % if (maxRHO > 1100|| minRHO < 1000)
     %     disp('RHO out of range')
     %     disp (['minRHO = ',num2str(minRHO,4), ...
     %           ', maxRHO = ',num2str(maxRHO,4) ]);
     % end
 
   %% get the average velocity in the v-direction
      varid = netcdf.inqVarID(ncid,'VVEL'); % cm/s
      VVELd(:,:,:) = netcdf.getVar(ncid,varid,'double'); 
      VVELp = permute(VVELd(:,:,:),[2,1,3]);  %reorder to (ny,nx)
      VVELp = VVELp.*.01; % convert units to m/s
      VVEL = zeros(ny,nx,nz);
      VVEL(iocn) = VVELp(iocn);
     % maxVVEL = max(VVEL(iocn));
     % minVVEL = min(VVEL(iocn));
     % disp(['max VVEL = ',num2str(maxVVEL,4), ...
      %      ', min VVEL = ',num2str(minVVEL,4)]);
     % if (maxRHO > 1100|| minRHO < 1000)
     %     disp('RHO out of range')
     %     disp (['minRHO = ',num2str(minRHO,4), ...
     %           ', maxRHO = ',num2str(maxRHO,4) ]);
     % end

   %% get the average velocity in the w-direction
      varid = netcdf.inqVarID(ncid,'WVEL'); % cm/s
      WVELd(:,:,:) = netcdf.getVar(ncid,varid,'double'); 
      WVELp = permute(WVELd(:,:,:),[2,1,3]);  %reorder to (ny,nx)
      WVELp = WVELp.*.01; % convert units to m/s
      WVEL = zeros(ny,nx,nz);
      WVEL(iocn) = WVELp(iocn);
      %maxWVEL = max(WVEL(iocn));
      %minWVEL = min(WVEL(iocn));
      %disp(['max WVEL = ',num2str(maxWVEL,4), ...
      %      ', min WVEL = ',num2str(minWVEL,4)]);
     % if (maxRHO > 1100|| minRHO < 1000)
     %     disp('RHO out of range')
     %     disp (['minRHO = ',num2str(minRHO,4), ...
     %           ', maxRHO = ',num2str(maxRHO,4) ]);
     % end

   %% get the average velocity squared in the u-direction
      varid = netcdf.inqVarID(ncid,'UVEL2'); % cm^2/s^2
      UVEL2d(:,:,:) = netcdf.getVar(ncid,varid,'double'); 
      UVEL2p = permute(UVEL2d(:,:,:),[2,1,3]);  %reorder to (ny,nx)
      UVEL2p = UVEL2p.*(.01)^2; % convert units to m^2/s^2
      UVEL2 = zeros(ny,nx,nz);
      UVEL2(iocn) = UVEL2p(iocn);
    %  maxUVEL2 = max(UVEL2(iocn));
     % minUVEL2 = min(UVEL2(iocn));
     % disp(['max UVEL2 = ',num2str(maxUVEL2,4), ...
     %       ', min UVEL2 = ',num2str(minUVEL2,4)]);
     % if (maxRHO > 1100|| minRHO < 1000)
     %     disp('RHO out of range')
     %     disp (['minRHO = ',num2str(minRHO,4), ...
     %           ', maxRHO = ',num2str(maxRHO,4) ]);
     % end

   %% get the average velocity squared in the v-direction
      varid = netcdf.inqVarID(ncid,'VVEL2'); % cm^2/s^2
      VVEL2d(:,:,:) = netcdf.getVar(ncid,varid,'double'); 
      VVEL2p = permute(VVEL2d(:,:,:),[2,1,3]);  %reorder to (ny,nx)
      VVEL2p = VVEL2p.*(.01)^2; % convert units to m^2/s^2
      VVEL2 = zeros(ny,nx,nz);
      VVEL2(iocn) = VVEL2p(iocn);
     % maxVVEL2 = max(VVEL2(iocn));
     % minVVEL2 = min(VVEL2(iocn));
      %disp(['max VVEL2 = ',num2str(maxVVEL2,4), ...
      %      ', min VVEL2 = ',num2str(minVVEL2,4)]);
     % if (maxRHO > 1100|| minRHO < 1000)
     %     disp('RHO out of range')
     %     disp (['minRHO = ',num2str(minRHO,4), ...
     %           ', maxRHO = ',num2str(maxRHO,4) ]);
     % end


   %% get the average velocity squared in the w-direction
      varid = netcdf.inqVarID(ncid,'WVEL2'); % cm^2/s^2
      WVEL2d(:,:,:) = netcdf.getVar(ncid,varid,'double'); 
      WVEL2p = permute(WVEL2d(:,:,:),[2,1,3]);  %reorder to (ny,nx)
      WVEL2p = WVEL2p.*(.01)^2; % convert units to m^2/s^2
      WVEL2 = zeros(ny,nx,nz);
      WVEL2(iocn) = WVEL2p(iocn);
      %maxWVEL2 = max(WVEL2(iocn));
      %minWVEL2 = min(WVEL2(iocn));
     % disp(['max WVEL2 = ',num2str(maxWVEL2,4), ...
      %      ', min WVEL2 = ',num2str(minWVEL2,4)]);
     % if (maxRHO > 1100|| minRHO < 1000)
     %     disp('RHO out of range')
     %     disp (['minRHO = ',num2str(minRHO,4), ...
     %           ', maxRHO = ',num2str(maxRHO,4) ]);
     % end

 %% get the average wind-stress in the x-direction
      varid = netcdf.inqVarID(ncid,'TAUX'); % dyne/(cm)^2
      TAUXd(:,:) = netcdf.getVar(ncid,varid,'double'); 
      TAUXp = permute(TAUXd(:,:),[2,1]);  %reorder to (ny,nx)
      TAUXp = TAUXp.*(1.0197e-6)*((100)^2); % convert units to 
                                            % kilogram forcing/m2
      TAUX = zeros(ny,nx);
      TAUX(isurfocn) = TAUXp(isurfocn);
    %  maxTAUX = max(TAUX(isurfocn));
    %  minTAUX = min(TAUX(isurfocn));
     % disp(['max TAUX = ',num2str(maxTAUX,4), ...
     %       ', min TAUX = ',num2str(minTAUX,4)]);
     % if (maxRHO > 1100|| minRHO < 1000)
     %     disp('RHO out of range')
     %     disp (['minRHO = ',num2str(minRHO,4), ...
     %           ', maxRHO = ',num2str(maxRHO,4) ]);
     % end

 %% get the average wind-stress in the y-direction
      varid = netcdf.inqVarID(ncid,'TAUY'); % dyne/(cm)^2
      TAUYd(:,:) = netcdf.getVar(ncid,varid,'double'); 
      TAUYp = permute(TAUYd(:,:),[2,1]);  %reorder to (ny,nx)
      TAUYp = TAUYp.*(1.0197e-6)*((100)^2); % convert units to 
                                            % kilogram forcing/m2
      TAUY = zeros(ny,nx);
      TAUY(isurfocn) = TAUYp(isurfocn);
     % maxTAUY = max(TAUY(isurfocn));
     % minTAUY = min(TAUY(isurfocn));
     % disp(['max TAUY = ',num2str(maxTAUY,4), ...
      %      ', min TAUY = ',num2str(minTAUY,4)]);
     % if (maxRHO > 1100|| minRHO < 1000)
     %     disp('RHO out of range')
     %     disp (['minRHO = ',num2str(minRHO,4), ...
     %           ', maxRHO = ',num2str(maxRHO,4) ]);
     % end

  %% get KAPPA ISOP coefficient
      varid = netcdf.inqVarID(ncid,'KAPPA_ISOP'); % 
      KId(:,:,:) = netcdf.getVar(ncid,varid,'double'); 
      KIp = permute(KId(:,:,:),[2,1,3]);  %reorder to (ny,nx)
      KIp = KIp ./(100^2); % convert cm^2/s to m^2/s
      KAPPA_ISOP = zeros(ny,nx,nz);
      KAPPA_ISOP(iocn) = KIp(iocn);
     % maxKI = max(KAPPA_ISOP(iocn));
      %minKI = min(KAPPA_ISOP(iocn));
     % disp(['max KAPPA ISOP = ',num2str(maxKI,4), ...
      %      ', min KAPPA ISOP = ',num2str(minKI,4)]);
     % if (maxRHO > 1100|| minRHO < 1000)
     %     disp('RHO out of range')
     %     disp (['minRHO = ',num2str(minRHO,4), ...
     %           ', maxRHO = ',num2str(maxRHO,4) ]);
     % end

 %% get KAPPA THIC coefficient
      varid = netcdf.inqVarID(ncid,'KAPPA_THIC'); % 
      KTd(:,:,:) = netcdf.getVar(ncid,varid,'double'); 
      KTp = permute(KTd(:,:,:),[2,1,3]);  %reorder to (ny,nx)
      KTp = KTp ./(100^2); % convert cm^2/s to m^2/s
      KAPPA_THIC = zeros(ny,nx,nz);
      KAPPA_THIC(iocn) = KTp(iocn);
     % maxKT = max(KAPPA_THIC(iocn));
     % minKT = min(KAPPA_THIC(iocn));
     % disp(['max KAPPA THIC = ',num2str(maxKT,4), ...
      %      ', min KAPPA THIC = ',num2str(minKT,4)]);
     % if (maxRHO > 1100|| minRHO < 1000)
     %     disp('RHO out of range')
     %     disp (['minRHO = ',num2str(minRHO,4), ...
     %           ', maxRHO = ',num2str(maxRHO,4) ]);
     % end

 %% get Horizontal diffusion  coefficient
      varid = netcdf.inqVarID(ncid,'HOR_DIFF'); % 
      HDd(:,:,:) = netcdf.getVar(ncid,varid,'double'); 
      HDp = permute(HDd(:,:,:),[2,1,3]);  %reorder to (ny,nx)
      HDp = HDp ./(100^2); % convert cm^2/s to m^2/s
      HOR_DIFF = zeros(ny,nx,nz);
      HOR_DIFF(iocn) = HDp(iocn);
     % maxHD = max(HOR_DIFF(iocn));
      %minHD = min(HOR_DIFF(iocn));
      %disp(['max HOR DIFF = ',num2str(maxHD,4), ...
      %      ', min HOR DIFF = ',num2str(minHD,4)]);
     % if (maxRHO > 1100|| minRHO < 1000)
     %     disp('RHO out of range')
     %     disp (['minRHO = ',num2str(minRHO,4), ...
     %           ', maxRHO = ',num2str(maxRHO,4) ]);
     % end

 %% get UISOP Bolus velocity in  the grid X-direction
      varid = netcdf.inqVarID(ncid,'UISOP'); % 
      UId(:,:,:) = netcdf.getVar(ncid,varid,'double'); 
      UIp = permute(UId(:,:,:),[2,1,3]);  %reorder to (ny,nx)
      UIp = UIp./100; % convert from cm/s to m/s
      UISOP = zeros(ny,nx,nz);
      UISOP(iocn) = UIp(iocn);
      %maxUI = max(UISOP(iocn));
      %minUI = min(UISOP(iocn));
     % disp(['max UISOP = ',num2str(maxUI,4), ...
     %       ', min UISOP = ',num2str(minUI,4)]);
     % if (maxRHO > 1100|| minRHO < 1000)
     %     disp('RHO out of range')
     %     disp (['minRHO = ',num2str(minRHO,4), ...
     %           ', maxRHO = ',num2str(maxRHO,4) ]);
     % end

 %% get VISOP Bolus velocity in  the grid Y-direction
      varid = netcdf.inqVarID(ncid,'VISOP'); % 
      VId(:,:,:) = netcdf.getVar(ncid,varid,'double'); 
      VIp = permute(VId(:,:,:),[2,1,3]);  %reorder to (ny,nx)
      VIp = VIp ./100; % convert from cm/s to m/s
      VISOP = zeros(ny,nx,nz);
      VISOP(iocn) = VIp(iocn);
     % maxVI = max(VISOP(iocn));
     % minVI = min(VISOP(iocn));
     % disp(['max VISOP = ',num2str(maxVI,4), ...
      %      ', min VISOP = ',num2str(minVI,4)]);
     % if (maxRHO > 1100|| minRHO < 1000)
     %     disp('RHO out of range')
     %     disp (['minRHO = ',num2str(minRHO,4), ...
     %           ', maxRHO = ',num2str(maxRHO,4) ]);
     % end

 %% get WISOP Bolus velocity in  the grid Z-direction
      varid = netcdf.inqVarID(ncid,'WISOP'); % 
      WId(:,:,:) = netcdf.getVar(ncid,varid,'double'); 
      WIp = permute(WId(:,:,:),[2,1,3]);  %reorder to (ny,nx)
      WIp = WIp ./100; % convert cm/s to m/s
      WISOP = zeros(ny,nx,nz);
      WISOP(iocn) = WIp(iocn);
     % maxWI = max(WISOP(iocn));
     % minWI = min(WISOP(iocn));
     % disp(['max WISOP = ',num2str(maxWI,4), ...
      %      ', min WISOP = ',num2str(minWI,4)]);
     % if (maxRHO > 1100|| minRHO < 1000)
     %     disp('RHO out of range')
     %     disp (['minRHO = ',num2str(minRHO,4), ...
     %           ', maxRHO = ',num2str(maxRHO,4) ]);
     % end

         %% get the MOC, the water transport in the Northwards direction, [Sv] 

           varid = netcdf.inqVarID(ncid,'MOC'); % 
           % MOC has dimensions of [lat, z, component,t_reg]
           MOCd(:,:,:,:) = netcdf.getVar(ncid,varid,'double');
           MOC = MOCd; 
 
           % separate the basins
           gbl = squeeze(MOC(:,:,:,1) );
           Atl = squeeze(MOC(:,:,:,2) );
           % add the components for a total MOC
           MOC_gbl = permute(squeeze(sum(gbl,3) ),[2,1]);
           MOC_Atl = permute(squeeze(sum(Atl,3) ),[2,1]);
           % Indio-Pacific is the difference
           MOC_IndPac = MOC_gbl - MOC_Atl;
      
           % get the MOC (zonal average) dimentions
           varid = netcdf.inqVarID(ncid,'lat_aux_grid'); % 
           lat_aux_grid = netcdf.getVar(ncid,varid,'double');
           varid = netcdf.inqVarID(ncid,'moc_z'); % 
           moc_z_d = netcdf.getVar(ncid,varid,'double');
           moc_z = moc_z_d ./100;  %cm to m, positive

%% convert the data set u, v, w, to sperical coordinates Umag, 
%    the horizontal angle Ulambda, and the vertical angle  Uphi.
%    The magnitude (Umag) and angle of elevation (Uele) are saved.
%    Both the horizontal angle relative to the grid box East(Ugbazi),
%    and relative to due-East (Uaximuth) are saved.
%    Includes the bolus component UISOP.
     Unet = UVEL + UISOP;
     Vnet = VVEL + VISOP;
     Wnet = WVEL + WISOP;
     % convert to spherical coordinates relative to the grid-box
     [Umag, Ugbazi, Uele] = convert_to_spherical(Unet, Vnet, Wnet); 
     % calculate U for due-East = 0 (ie, adjust for grid-box angle alignment)
     [Uazimuth] = make_true_East(Ugbazi,MET);

  %% save monthly variables to file   
      eval(['save ',F.noc_dir,out_file_base,int2str(seq_num),'.mat ', ...
		   ' SSH TEMP SALT RHO PD IFRAC ', ...         
                   ' UVEL VVEL WVEL UVEL2 VVEL2 WVEL2 ', ...
                   ' TAUX TAUY ', ...
                   ' ECOSYS_IFRAC ECOSYS_XKW ECOSYS_ATM_PRESS ', ...
                   ' O2 O2_PRODUCTION O2_CONSUMPTION ', ...
                   ' SCHMIDT_O2 O2SAT DOP PO4 ', ...                 
                   ' NITRIF SedDenitrif ', ...
                   ' PH ALK SCHMIDT_CO2 ATM_CO2 ', ...
                   ' KAPPA_ISOP KAPPA_THIC HOR_DIFF ', ...
                   ' TLT INT_DEPTH ', ...                  
                   ' UISOP VISOP WISOP ', ...
                   ' Umag Ugbazi Uele Uazimuth ', ...
                   ' IAGE ', ...
                   ' MOC MOC_gbl MOC_Atl MOC_IndPac lat_aux_grid moc_z ', ... 
                   ' -v7.3 ']);
%    ' CFC_IFRAC CFC_XKW CFC_ATM_PRESS CFC11 STF_CFC11 ', ...
%	 ' POC_remin DOC_remin DOCr_remin ', ...
      disp(['saved period ',int2str(seq_num),' state data to file']);

      % close the netcdf file
      netcdf.close(ncid);
 return

end % function build_model_state_data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% convert_to_spherical %%%%%%%%%%%%%%%%%%%%%%%
function [Umag, Ugbazi, Uele] = convert_to_spherical(u,v,w)
% converts relative to the grid-box
% Ugbazi is -pi to +pi
% Uele is -pi/2 to +pi/2
  [Ugbazi,Uele,Umag] = cart2sph(u,v,w);
  return
end % function convert_U_to_spherical
 
%% make true East azimuth %%%%%%%%%%%%%%%%%%%%%%%% 
function  [Uazimuth] = make_true_East(Ugbazi,MET)
  % adjust the Ugbazi to true-East orientation from the orientation of the grid-box
  persistent ANGLE_3D
  if isempty(ANGLE_3D)
            [ny,nx,nz] = size(MET.MASK);
            ANGLE_3D = 0*MET.MASK;
            for k = 1:nz
		      ANGLE_3D(:,:,k) = MET.ANGLE(:,:,1);
            end
  end % if empty
  
  Uazimuth = Ugbazi + ANGLE_3D;
  i_too_big = find(Uazimuth > pi);
  if ~isempty(i_too_big)
      Uazimuth(i_too_big)  = Uazimuth(i_too_big) - 2*pi;
  end % if not empty
  i_too_sm = find( Uazimuth < -pi);
  if ~isempty(i_too_sm)    
      Uazimuth (i_too_sm) = Uazimuth(i_too_sm) + 2*pi;
  end % if not empty
  return
end % function make_true_East
