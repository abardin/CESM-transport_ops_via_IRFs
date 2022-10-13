function Build_mo_climate_state(F)
% Gathers the ocean state variables from the /noc directory, and puts them
% into .mat files climatologically-averaged per month, 
% and an annualized climatological avg.
% Data in each file:
%  SSH, TEMP (potential temprature), PD (potential density)
%  SALT (salinity), RHO (in-situ density), IFRAC (sea-ice fraction)
%  and others, as needed for O2 simultaion 
% The data are in dimensions of ny X nx X nz for 3_D data, and ny X nx for
% 2-D data. 
% The structure F needs to have
% F.ops_dir = directory for the output data set;
% F.noc_dir = directory for the supper-set of the individual months ops and state data reside 
% F.zeroyr  = year before the first year of the data set
% F.nyears = number of sequential years to be made into a climatological average
% F.firstyr = first year to be processed; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  out_dir = F.ops_dir;                  % into the final ops directory
  % make dir if DNE
  if ~exist(out_dir,'dir')
      eval(['!mkdir ',out_dir]);
  end

  state_file_base = 'model_state_data_';

  % set up dimensions from MET
  eval (['load ',F.noc_dir,'MET.mat  MET']);
  eval (['save ',F.ops_dir,'MET.mat MET -v7.3']);
  [ny,nx,nz] = size(MET.MASK);
  iocn = MET.iocn;
  surf_mask = MET.MASK(:,:,1);
  isurfocn = find(surf_mask ~= 0);

  % load 1 state variable file to get the odd dimensions (MOC)
  fn = ([F.noc_dir,state_file_base,'1','.mat']);
  eval(['load ', fn, ' lat_aux_grid moc_z ']);
  n_moc_z   = length(moc_z);
  n_lat_aux = length(lat_aux_grid);
  % save in MET
  MET.moc_z = moc_z;
  MET.lat_aux_grid = lat_aux_grid;
  eval (['save ',F.ops_dir,'MET.mat MET -v7.3']);
  % initialize annual averages
  ALK_sum              = zeros(ny,nx,nz);
  ATM_CO2_sum          = zeros(ny,nx);
%  CFC11_sum            = zeros(ny,nx,nz); 
%  CFC_IFRAC_sum        = zeros(ny,nx);
%  CFC_XKW_sum          = zeros(ny,nx);
%  CFC_ATM_PRESS_sum    = zeros(ny,nx);
%  STF_CFC11_sum        = zeros(ny,nx);
  ECOSYS_ATM_PRESS_sum = zeros(ny,nx);
  ECOSYS_IFRAC_sum     = zeros(ny,nx);
  ECOSYS_XKW_sum       = zeros(ny,nx);
  HOR_DIFF_sum         = zeros(ny,nx,nz);
  IFRAC_sum            = zeros(ny,nx);
  INT_DEPTH_sum        = zeros(ny,nx); 
  KAPPA_ISOP_sum       = zeros(ny,nx,nz);
  KAPPA_THIC_sum       = zeros(ny,nx,nz);  
  O2_sum               = zeros(ny,nx,nz);
  O2SAT_sum            = zeros(ny,nx);
  O2_CONSUMPTION_sum   = zeros(ny,nx,nz);
  O2_PRODUCTION_sum    = zeros(ny,nx,nz);
  PD_sum               = zeros(ny,nx,nz);
  PH_sum               = zeros(ny,nx);
  RHO_sum              = zeros(ny,nx,nz);
  SALT_sum             = zeros(ny,nx,nz);
  SCHMIDT_CO2_sum      = zeros(ny,nx);
  SCHMIDT_O2_sum       = zeros(ny,nx);
  SSH_sum              = zeros(ny,nx);
  TAUX_sum             = zeros(ny,nx); 
  TAUY_sum             = zeros(ny,nx);
  TEMP_sum             = zeros(ny,nx,nz);
  TLT_sum              = zeros(ny,nx);
  UISOP_sum            = zeros(ny,nx,nz);
  UVEL_sum             = zeros(ny,nx,nz);
  UVEL2_sum            = zeros(ny,nx,nz);
  Umag_sum             = zeros(ny,nx,nz);
  Ugbazi_sum           = zeros(ny,nx,nz);
  Uele_sum             = zeros(ny,nx,nz);
  Uazimuth_sum         = zeros(ny,nx,nz); 
  VISOP_sum            = zeros(ny,nx,nz);
  VVEL_sum             = zeros(ny,nx,nz);
  VVEL2_sum            = zeros(ny,nx,nz);
  WISOP_sum            = zeros(ny,nx,nz);
  WVEL_sum             = zeros(ny,nx,nz);
  WVEL2_sum            = zeros(ny,nx,nz);
  MOC_sum              = zeros(n_lat_aux,n_moc_z,3,2);
  MOC_gbl_sum          = zeros(n_moc_z,n_lat_aux);
  MOC_Atl_sum          = zeros(n_moc_z,n_lat_aux);
  MOC_IndPac_sum       = zeros(n_moc_z,n_lat_aux);

  total_days = 0;
  
% for each seasonal month (Jan,Feb, etc)
  for smo =1:12
          F.clim_mo_state_fn = ([F.ops_dir,state_file_base,int2str(smo),'.mat']);
          make_clim_mo_state(smo,F,MET); % outputs file for climatological month
          % read the data from the clim_mo_file;
          eval(['load ',F.clim_mo_state_fn]);

          % get number of days per this month
          eval(['load ',F.noc_dir,'MTM',int2str(smo),'.mat', ' num_days '])

          % add to sum for super-annual
          total_days = total_days + num_days;
          ALK_sum              = ALK_sum              + num_days.*ALK;
          ATM_CO2_sum          = ATM_CO2_sum          + num_days.*ATM_CO2;
%          CFC11_sum            = CFC11_sum            + num_days.*CFC11; 
%          CFC_IFRAC_sum        = CFC_IFRAC_sum        + num_days.*CFC_IFRAC;
%          CFC_XKW_sum          = CFC_XKW_sum          + num_days.*CFC_XKW;
%          CFC_ATM_PRESS_sum    = CFC_ATM_PRESS_sum    + num_days.*CFC_ATM_PRESS;
%          STF_CFC11_sum        = STF_CFC11_sum        + num_days.*STF_CFC11;
          ECOSYS_ATM_PRESS_sum = ECOSYS_ATM_PRESS_sum + num_days.*ECOSYS_ATM_PRESS;
          ECOSYS_IFRAC_sum     = ECOSYS_IFRAC_sum     + num_days.*ECOSYS_IFRAC;
          ECOSYS_XKW_sum       = ECOSYS_XKW_sum       + num_days.*ECOSYS_XKW;
          HOR_DIFF_sum         = HOR_DIFF_sum         + num_days.*HOR_DIFF;
          IFRAC_sum            = IFRAC_sum            + num_days.*IFRAC;
          INT_DEPTH_sum        = INT_DEPTH_sum        + num_days.*INT_DEPTH; 
          KAPPA_ISOP_sum       =  KAPPA_ISOP_sum      + num_days.*KAPPA_ISOP;
          KAPPA_THIC_sum       = KAPPA_THIC_sum       + num_days.*KAPPA_THIC;  
          O2_sum               = O2_sum               + num_days.*O2;
          O2SAT_sum            = O2SAT_sum            + num_days.*O2SAT;
          O2_CONSUMPTION_sum   = O2_CONSUMPTION_sum   + num_days.*O2_CONSUMPTION;
          O2_PRODUCTION_sum    = O2_PRODUCTION_sum    + num_days.*O2_PRODUCTION;
          PD_sum               = PD_sum               + num_days.*PD;
          PH_sum               = PH_sum               + num_days.*PH;
          RHO_sum              = RHO_sum              + num_days.*RHO;
          SALT_sum             = SALT_sum             + num_days.*SALT;
          SCHMIDT_CO2_sum      = SCHMIDT_CO2_sum      + num_days.*SCHMIDT_CO2;
          SCHMIDT_O2_sum       = SCHMIDT_O2_sum       + num_days.*SCHMIDT_O2;
          SSH_sum              = SSH_sum              + num_days.*SSH;
          TAUX_sum             = TAUX_sum             + num_days.*TAUX; 
          TAUY_sum             = TAUY_sum             + num_days.*TAUY;
          TEMP_sum             = TEMP_sum             + num_days.*TEMP;
          TLT_sum              = TLT_sum              + num_days.*TLT;
          UISOP_sum            = UISOP_sum            + num_days.*UISOP;
          UVEL_sum             = UVEL_sum             + num_days.*UVEL;
          UVEL2_sum            = UVEL2_sum            + num_days.*UVEL2;
          Umag_sum             = Umag_sum             + num_days.*Umag;
          Ugbazi_sum           = Ugbazi_sum           + num_days.*Ugbazi;
          Uele_sum             = Uele_sum             + num_days.*Uele;
          Uazimuth_sum         = Uazimuth_sum         + num_days.*Uazimuth; 
          VISOP_sum            = VISOP_sum            + num_days.*VISOP;
          VVEL_sum             = VVEL_sum             + num_days.*VVEL;
          VVEL2_sum            = VVEL2_sum            + num_days.*VVEL2;
          WISOP_sum            = WISOP_sum            + num_days.*WISOP;
          WVEL_sum             = WVEL_sum             + num_days.*WVEL;
          WVEL2_sum            = WVEL2_sum            + num_days.*WVEL2;
          MOC_sum              = MOC_sum              + num_days.*MOC;
          MOC_gbl_sum          = MOC_gbl_sum          + num_days.*MOC_gbl;
          MOC_Atl_sum          = MOC_Atl_sum          + num_days.*MOC_Atl;
          MOC_IndPac_sum       = MOC_IndPac_sum       + num_days.*MOC_IndPac;

  end % for each seasonal month

  % make super_annual_avgs
  if total_days ~= 365
    disp(['wrong number of days per year'])
    keyboard
  end

  ALK              = ALK_sum./total_days;
  ATM_CO2          = ATM_CO2_sum./total_days;
%  CFC11            = CFC11_sum./total_days;
%  CFC_IFRAC        = CFC_IFRAC_sum./total_days;
%  CFC_XKW          = CFC_XKW_sum./total_days;
%  CFC_ATM_PRESS    = CFC_ATM_PRESS_sum./total_days;
%  STF_CFC11        = STF_CFC11_sum./total_days;
  ECOSYS_ATM_PRESS = ECOSYS_ATM_PRESS_sum./total_days;
  ECOSYS_IFRAC     = ECOSYS_IFRAC_sum./total_days;    
  ECOSYS_XKW       = ECOSYS_XKW_sum./total_days;     
  HOR_DIFF         = HOR_DIFF_sum./total_days;       
  IFRAC            = IFRAC_sum./total_days;          
  INT_DEPTH        = INT_DEPTH_sum./total_days;    
  KAPPA_ISOP       = KAPPA_ISOP_sum./total_days;    
  KAPPA_THIC       = KAPPA_THIC_sum./total_days;      
  O2               = O2_sum./total_days;       
  O2SAT            = O2SAT_sum./total_days;          
  O2_CONSUMPTION   = O2_CONSUMPTION_sum./total_days; 
  O2_PRODUCTION    = O2_PRODUCTION_sum./total_days;   
  PD               = PD_sum./total_days;             
  PH               = PH_sum./total_days;             
  RHO              = RHO_sum./total_days;           
  SALT             = SALT_sum./total_days;          
  SCHMIDT_CO2      = SCHMIDT_CO2_sum./total_days;   
  SCHMIDT_O2       = SCHMIDT_O2_sum./total_days;    
  SSH              = SSH_sum./total_days;            
  TAUX             = TAUX_sum./total_days;           
  TAUY             = TAUY_sum./total_days;           
  TEMP             = TEMP_sum./total_days;           
  TLT              = TLT_sum./total_days;           
  UISOP            = UISOP_sum./total_days;          
  UVEL             = UVEL_sum./total_days;            
  UVEL2            = UVEL2_sum./total_days;
  Umag             = Umag_sum./total_days;
  Ugbazi           = Ugbazi_sum./total_days;
  Uele             = Uele_sum./total_days;
  Uazimuth         = Uazimuth_sum./total_days;           
  VISOP            = VISOP_sum./total_days;          
  VVEL             = VVEL_sum./total_days;        
  VVEL2            = VVEL2_sum./total_days;         
  WISOP            = WISOP_sum./total_days;         
  WVEL             = WVEL_sum./total_days;          
  WVEL2            = WVEL2_sum./total_days;         
  MOC              = MOC_sum./total_days;   
  MOC_gbl          = MOC_gbl_sum./total_days;       
  MOC_Atl          = MOC_Atl_sum./total_days;        
  MOC_IndPac       = MOC_IndPac_sum./total_days;
      
  % save as climatological avg
  eval(['save ',F.ops_dir,'model_state_data_0.mat ', ...
	       '   ALK ', ... 
	       '   ATM_CO2 ', ...     
	       '   ECOSYS_ATM_PRESS ', ... 
	       '   ECOSYS_IFRAC ', ...        
	       '   ECOSYS_XKW ', ...            
	       '   HOR_DIFF ', ...              
	       '   IFRAC ', ...                   
	       '   INT_DEPTH ', ...          
	       '   KAPPA_ISOP ', ...        
	       '   KAPPA_THIC ', ...       
	       '   O2  ', ...                  
	       '   O2SAT ', ...                  
	       '   O2_CONSUMPTION ', ...  
	       '   O2_PRODUCTION ', ...     
	       '   PD ', ...                       
	       '   PH ', ...                         
	       '   RHO ', ...                     
	       '   SALT ', ...                  
	       '   SCHMIDT_CO2 ', ...       
	       '   SCHMIDT_O2 ', ...        
	       '   SSH ', ...                      
	       '   TAUX ', ...                    
	       '   TAUY ', ...                     
	       '   TEMP ', ...                     
	       '   TLT ', ...                     
	       '   UISOP ', ...                 
	       '   UVEL ', ...                     
	       '   UVEL2 ', ...
               '   Umag Ugbazi Uele Uazimuth ', ...                   
	       '   VISOP ', ...                   
	       '   VVEL ', ...                  
	       '   VVEL2 ', ...                  
	       '   WISOP ', ...                  
	       '   WVEL ', ...                    
	       '   WVEL2 ', ...
               '   MOC ', ...      
               '   MOC_gbl ', ...       
               '   MOC_Atl ', ...       
               '   MOC_IndPac ', ...
               '   lat_aux_grid ', ...
               '   moc_z ', ...
               '   -v7.3 ' ]);

%               '   CFC11 CFC_IFRAC CFC_XKW CFC_ATM_PRESS STF_CFC11 ', ...   

  return
end % Build_mo_climate_state

%% make_clim_mo_state %%%%%%%%%%%%%%%%%%%%
function make_clim_mo_state(smo,F,MET)
  % get thedimensions
  [ny,nx,nz] = size(MET.MASK);
  n_moc_z = length(MET.moc_z);
  n_lat_aux = length(MET.lat_aux_grid);

 % initialize annual averages
  ALK_sum              = zeros(ny,nx,nz);
  ATM_CO2_sum          = zeros(ny,nx);
%  CFC11_sum            = zeros(ny,nx,nz); 
%  CFC_IFRAC_sum        = zeros(ny,nx);
%  CFC_XKW_sum          = zeros(ny,nx);
%  CFC_ATM_PRESS_sum    = zeros(ny,nx);
%  STF_CFC11_sum        = zeros(ny,nx);
  ECOSYS_ATM_PRESS_sum = zeros(ny,nx);
  ECOSYS_IFRAC_sum     = zeros(ny,nx);
  ECOSYS_XKW_sum       = zeros(ny,nx);
  HOR_DIFF_sum         = zeros(ny,nx,nz);
  IFRAC_sum            = zeros(ny,nx);
  INT_DEPTH_sum        = zeros(ny,nx); 
  KAPPA_ISOP_sum       = zeros(ny,nx,nz);
  KAPPA_THIC_sum       = zeros(ny,nx,nz);  
  O2_sum               = zeros(ny,nx,nz);
  O2SAT_sum            = zeros(ny,nx);
  O2_CONSUMPTION_sum   = zeros(ny,nx,nz);
  O2_PRODUCTION_sum    = zeros(ny,nx,nz);
  PD_sum               = zeros(ny,nx,nz);
  PH_sum               = zeros(ny,nx);
  RHO_sum              = zeros(ny,nx,nz);
  SALT_sum             = zeros(ny,nx,nz);
  SCHMIDT_CO2_sum      = zeros(ny,nx);
  SCHMIDT_O2_sum       = zeros(ny,nx);
  SSH_sum              = zeros(ny,nx);
  TAUX_sum             = zeros(ny,nx); 
  TAUY_sum             = zeros(ny,nx);
  TEMP_sum             = zeros(ny,nx,nz);
  TLT_sum              = zeros(ny,nx);
  UISOP_sum            = zeros(ny,nx,nz);
  UVEL_sum             = zeros(ny,nx,nz);
  UVEL2_sum            = zeros(ny,nx,nz);
  Umag_sum             = zeros(ny,nx,nz);
  Ugbazi_sum           = zeros(ny,nx,nz);
  Uele_sum             = zeros(ny,nx,nz);
  Uazimuth_sum         = zeros(ny,nx,nz); 
  VISOP_sum            = zeros(ny,nx,nz);
  VVEL_sum             = zeros(ny,nx,nz);
  VVEL2_sum            = zeros(ny,nx,nz);
  WISOP_sum            = zeros(ny,nx,nz);
  WVEL_sum             = zeros(ny,nx,nz);
  WVEL2_sum            = zeros(ny,nx,nz);
  MOC_sum              = zeros(n_lat_aux,n_moc_z,3,2);
  MOC_gbl_sum          = zeros(n_moc_z,n_lat_aux);
  MOC_Atl_sum          = zeros(n_moc_z,n_lat_aux);
  MOC_IndPac_sum       = zeros(n_moc_z,n_lat_aux);

  % determine month_number for first year

  % sum data for all years seasonal month
  for year = F.firstyr:1:(F.firstyr + F.nyears - 1)
      mo_num = (year - F.zeroyr - 1)*12 + smo;
      eval(['load ',F.noc_dir,'model_state_data_',int2str(mo_num),'.mat']);

      % sum the data


          ALK_sum              = ALK_sum              + ALK;
          ATM_CO2_sum          = ATM_CO2_sum          + ATM_CO2;
%          CFC11_sum            = CFC11_sum            + CFC11; 
%          CFC_IFRAC_sum        = CFC_IFRAC_sum        + CFC_IFRAC;
%          CFC_XKW_sum          = CFC_XKW_sum          + CFC_XKW;
%          CFC_ATM_PRESS_sum    = CFC_ATM_PRESS_sum    + CFC_ATM_PRESS;
%          STF_CFC11_sum        = STF_CFC11_sum        + STF_CFC11;
          ECOSYS_ATM_PRESS_sum = ECOSYS_ATM_PRESS_sum + ECOSYS_ATM_PRESS;
          ECOSYS_IFRAC_sum     = ECOSYS_IFRAC_sum     + ECOSYS_IFRAC;
          ECOSYS_XKW_sum       = ECOSYS_XKW_sum       + ECOSYS_XKW;
          HOR_DIFF_sum         = HOR_DIFF_sum         + HOR_DIFF;
          IFRAC_sum            = IFRAC_sum            + IFRAC;
          INT_DEPTH_sum        = INT_DEPTH_sum        + INT_DEPTH; 
          KAPPA_ISOP_sum       = KAPPA_ISOP_sum       + KAPPA_ISOP;
          KAPPA_THIC_sum       = KAPPA_THIC_sum       + KAPPA_THIC;  
          O2_sum               = O2_sum               + O2;
          O2SAT_sum            = O2SAT_sum            + O2SAT;
          O2_CONSUMPTION_sum   = O2_CONSUMPTION_sum   + O2_CONSUMPTION;
          O2_PRODUCTION_sum    = O2_PRODUCTION_sum    + O2_PRODUCTION;
          PD_sum               = PD_sum               + PD;
          PH_sum               = PH_sum               + PH;
          RHO_sum              = RHO_sum              + RHO;
          SALT_sum             = SALT_sum             + SALT;
          SCHMIDT_CO2_sum      = SCHMIDT_CO2_sum      + SCHMIDT_CO2;
          SCHMIDT_O2_sum       = SCHMIDT_O2_sum       + SCHMIDT_O2;
          SSH_sum              = SSH_sum              + SSH;
          TAUX_sum             = TAUX_sum             + TAUX; 
          TAUY_sum             = TAUY_sum             + TAUY;
          TEMP_sum             = TEMP_sum             + TEMP;
          TLT_sum              = TLT_sum              + TLT;
          UISOP_sum            = UISOP_sum            + UISOP;
          UVEL_sum             = UVEL_sum             + UVEL;
          UVEL2_sum            = UVEL2_sum            + UVEL2;
          Umag_sum             = Umag_sum             + Umag;
          Ugbazi_sum           = Ugbazi_sum           + Ugbazi;
          Uele_sum             = Uele_sum             + Uele;
          Uazimuth_sum         = Uazimuth_sum         + Uazimuth; 
          VISOP_sum            = VISOP_sum            + VISOP;
          VVEL_sum             = VVEL_sum             + VVEL;
          VVEL2_sum            = VVEL2_sum            + VVEL2;
          WISOP_sum            = WISOP_sum            + WISOP;
          WVEL_sum             = WVEL_sum             + WVEL;
          WVEL2_sum            = WVEL2_sum            + WVEL2;
          MOC_sum              = MOC_sum              + MOC;
	  MOC_gbl_sum          = MOC_gbl_sum          + MOC_gbl;
	  MOC_Atl_sum          = MOC_Atl_sum          + MOC_Atl;
          MOC_IndPac_sum       = MOC_IndPac_sum       + MOC_IndPac;
  end % for each year

  % create data average
  ALK              = ALK_sum./F.nyears;
  ATM_CO2          = ATM_CO2_sum./F.nyears;
  ECOSYS_ATM_PRESS = ECOSYS_ATM_PRESS_sum./F.nyears;
  ECOSYS_IFRAC     = ECOSYS_IFRAC_sum./F.nyears;    
  ECOSYS_XKW       = ECOSYS_XKW_sum./F.nyears;     
  HOR_DIFF         = HOR_DIFF_sum./F.nyears;       
  IFRAC            = IFRAC_sum./F.nyears;          
  INT_DEPTH        = INT_DEPTH_sum./F.nyears;    
  KAPPA_ISOP       = KAPPA_ISOP_sum./F.nyears;    
  KAPPA_THIC       = KAPPA_THIC_sum./F.nyears;      
  O2               = O2_sum./F.nyears;       
  O2SAT            = O2SAT_sum./F.nyears;          
  O2_CONSUMPTION   = O2_CONSUMPTION_sum./F.nyears; 
  O2_PRODUCTION    = O2_PRODUCTION_sum./F.nyears;   
  PD               = PD_sum./F.nyears;             
  PH               = PH_sum./F.nyears;             
  RHO              = RHO_sum./F.nyears;           
  SALT             = SALT_sum./F.nyears;          
  SCHMIDT_CO2      = SCHMIDT_CO2_sum./F.nyears;   
  SCHMIDT_O2       = SCHMIDT_O2_sum./F.nyears;    
  SSH              = SSH_sum./F.nyears;            
  TAUX             = TAUX_sum./F.nyears;           
  TAUY             = TAUY_sum./F.nyears;           
  TEMP             = TEMP_sum./F.nyears;           
  TLT              = TLT_sum./F.nyears;           
  UISOP            = UISOP_sum./F.nyears;          
  UVEL             = UVEL_sum./F.nyears;            
  UVEL2            = UVEL2_sum./F.nyears;          
  VISOP            = VISOP_sum./F.nyears;          
  VVEL             = VVEL_sum./F.nyears;        
  VVEL2            = VVEL2_sum./F.nyears;         
  WISOP            = WISOP_sum./F.nyears;         
  WVEL             = WVEL_sum./F.nyears;          
  WVEL2            = WVEL2_sum./F.nyears;         
  MOC              = MOC_sum./F.nyears;           
  MOC_gbl          = MOC_gbl_sum./F.nyears;         
  MOC_Atl          = MOC_Atl_sum./F.nyears;          
  MOC_IndPac       = MOC_IndPac_sum./F.nyears;      

  % save to file
  % save as climatological avg
  eval(['save ', F.clim_mo_state_fn, ...
	       '   ALK ', ... 
	       '   ATM_CO2 ', ...
	       '   ECOSYS_ATM_PRESS ', ... 
	       '   ECOSYS_IFRAC ', ...        
	       '   ECOSYS_XKW ', ...            
	       '   HOR_DIFF ', ...              
	       '   IFRAC ', ...                   
	       '   INT_DEPTH ', ...          
	       '   KAPPA_ISOP ', ...        
	       '   KAPPA_THIC ', ...       
	       '   O2  ', ...                  
	       '   O2SAT ', ...                  
	       '   O2_CONSUMPTION ', ...  
	       '   O2_PRODUCTION ', ...     
	       '   PD ', ...                       
	       '   PH ', ...                         
	       '   RHO ', ...                     
	       '   SALT ', ...                  
	       '   SCHMIDT_CO2 ', ...       
	       '   SCHMIDT_O2 ', ...        
	       '   SSH ', ...                      
	       '   TAUX ', ...                    
	       '   TAUY ', ...                     
	       '   TEMP ', ...                     
	       '   TLT ', ...                     
	       '   UISOP ', ...                 
	       '   UVEL ', ...                     
	       '   UVEL2 ', ...
               '   Umag Ugbazi Uele Uazimuth ', ...                   
	       '   VISOP ', ...                   
	       '   VVEL ', ...                  
	       '   VVEL2 ', ...                  
	       '   WISOP ', ...                  
	       '   WVEL ', ...                    
	       '   WVEL2 ', ...
               '   MOC ', ...      
               '   MOC_gbl ', ...       
               '   MOC_Atl ', ...       
               '   MOC_IndPac ', ...
               '   lat_aux_grid ', ...
               '   moc_z ', ...
               '   -v7.3 ' ]);

%               '   CFC11 CFC_IFRAC CFC_XKW CFC_ATM_PRESS STF_CFC11 ', ...

  return
end % function make_clim_mo_state

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
