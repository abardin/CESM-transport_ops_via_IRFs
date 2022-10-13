function []=build_stencils_X3(F)
% build the IRF masks (stencils) that will be needed to do the conversion
% of the CESM archived IRF output (.nc) to build the transport operators. 

% THE FOLLOWING WILL NEED TO BE CHANGED FOR A DIFFERENT CONFIGURATION
% OR SET OF RUNS:
   % Directory and file for obtaining needed geometry
   reference_dir = F.ref_dir;
   data_file = F.ref_file;

   % Expected configuration
   ny_OGCM = 116; nx_OGCM = 100; nz_OGCM = 60; 
   % NOTE: nzz = nz + 1; an added layer in the offline model, to assure 
   % land on the boundary

   % Directory in which to put the stencils  
   stencil_dir = F.stencil_dir; 

   % Directory and text file for defining the product regions
   % The format for this file is described at the end of this file.
   % The product_txt_file is created by modifying the file
   % gx1v6_overflow using a text editor, and is therefor specific to the
   % configuration.  The file needs to be created before this program 
   % is run.
   product_txt_dir  = F.product_txt_dir; 
   product_txt_file = F.product_txt_file;

   ovf_cnt          = 4; % number of sill overflow areas defined
   % Hard-coded definitions for the source regions (1 column per region):
    % 1 is Denmark Straight
    % 2 is Faroe Bank Channel
    % 3 is Ross Sea
    % 4 is Weddell Sea
    %                       1    2    3     4
    PARAM.src_imin     = [   8   13   59   98]; % from overflow definition file
    PARAM.src_imax     = [  10   15   62  100]; % ditto
    PARAM.src_jmin     = [ 109   93    2    2]; % ditto
    PARAM.src_jmax     = [ 112   97    3    3]; % ditto
    PARAM.src_kmin     = [  33   38   34   36]; % ditto
    PARAM.src_kmax     = [  33   38   34   36]; % ditto
    PARAM.X_band_s     = [   5    5    6    5]; % width: imax-imin +2+1
    PARAM.Y_band_s     = [   6    7    4    4]; % width: jmax-jmin +2+1
    PARAM.Z_band_neg_s = [   2    2    2    2]; % depth band above kmax-kmin + 2 (within ocn)
    PARAM.Z_band_pos_s = [   2    2    2    2]; % depth band below kmax-kmin + 2 (within ocn)
    PARAM.Y_limit_s    = [-101  100    0    0]; % constraint in Y direction to assure no overlap
    
    PARAM.ent_imin =     [   3    7   59   96];
    PARAM.ent_imax =     [   5   10   61   98];
    PARAM.ent_jmin =     [ 104   93    6    6];
    PARAM.ent_jmax =     [ 107   96    8    8];
    PARAM.ent_kmin =     [  39   40   40   40];
    PARAM.ent_kmax =     [  39   40   40   40];
    PARAM.X_band_e     = [   5    6    6    5];
    PARAM.Y_band_e     = [   6    6    4    5];
    PARAM.Z_band_neg_e = [   2    2    2    2];
    PARAM.Z_band_pos_e = [   2    2    2    2];
    PARAM.Y_limit_e    = [-101  100    0    0];
    
% There are several figures that are "spy"s of the matrices as they are
% being constructed, which may be useful for debugging.  These are 
% commented out. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The purpose of Build_stencils is to build the necessary stencils to 
% support the building of the advection/diffusion operators
% for the offline model.  The stencils are the masks that are used to pick
% up the response data from the IRF impulse.
%
% The stencil building has been split off into a seperate program
% because it takes a significant amount of time to run, and only needs to
% be run once, or when there are changes to the circulation algorithm
% that would affect the stencils. 
%
% The program outputs ("saves") the three stencils:
%  g_stencil - for the global k-level IRF runs
%  s_stencil - for the overflow "source" IRF runs
%  e_stencil - for the overflow "entrainment" IRF runs
% These files can then be loaded into the build-operators program to 
% extract the data. 
%
% NOTE: for CESM DATA, NZZ = 61 = NZ + 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The definitions for the overflow source and entrainment regions have 
% been adapted from IRF_mod.f90 developed by Keith Lindsay of UCAR.
% The definition for the product region is from an edited version of 
% gx1v6_overflow, which is a definition file in the 
% /models/ocn/pop2/input_templates/.
% The names of the directories and files are hardcoded above.
%
% There are two kinds of files referenced: 
%    general geographic information from a history file of the OGCM
%    saved stencils
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%  Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
% Get data only needed to be read one time: area and dz

filename = [reference_dir, data_file]
ncid = netcdf.open(filename,'NC_NOWRITE');

varid = netcdf.inqVarID(ncid,'TAREA');
T_area(:,:) = netcdf.getVar(ncid,varid,'double' );
T_area = T_area ./ (100^2); % convert from sq centimeters to sq meters
[nx,ny] = size(T_area);
T_area = permute (T_area,[2,1]); % reorder to (ny,nx)

varid = netcdf.inqVarID(ncid,'dz'); % depth of grid box at each layer
dz(:,:) = netcdf.getVar(ncid,varid);
dz = dz./ 100; % convert from centimeters to meters
[nz,col_dim] = size(dz);
nzz = nz + 1;
dzz = zeros(nzz,1);
dzz(1:nz) = dz;
dz = dzz;

netcdf.close(ncid)

if ((nx == nx_OGCM) && (ny == ny_OGCM) && (nz == nz_OGCM) )
    % continue
else % unexpected dimensions
    disp ( 'check unexpected grid dimensions');
    keyboard
end % if dimensions as expected

PARAM.matrix_dim = ny*nx*nzz; % dimension of the sparse stencil matrices      
index=[1:PARAM.matrix_dim]';  %' used in get_x_data to calculate the vector index
                              % of a value in the matrix
PARAM.ovf_cnt = ovf_cnt;
% Make stencils for data collection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PARAM.ny = ny; PARAM.nx = nx; PARAM.nzz = nzz;
FILEDEF.product_txt_dir  = product_txt_dir;  
FILEDEF.product_txt_file = product_txt_file;
[g_irf_mask,s_irf_mask,e_irf_mask] = make_stencils_61(PARAM, FILEDEF); 

% Save stencils for use in data extraction program% %%%%%%%%%%%%%%%%%%%%%%
g_stencil_file = 'g_irf_mask.mat';
filename = [stencil_dir,g_stencil_file];
cmd = (['save ',filename, ' g_irf_mask -v7.3 ']);
eval(cmd);
disp(['saved ',filename]);

s_stencil_file = 's_irf_mask.mat';
filename = [stencil_dir,s_stencil_file];
cmd = (['save ', filename,' s_irf_mask -v7.3;']);
eval(cmd);
disp(['saved ',filename]);
 
e_stencil_file = 'e_irf_mask.mat';
filename = [stencil_dir,e_stencil_file];
cmd = (['save ',filename,' e_irf_mask -v7.3;']);
eval(cmd);
disp(['saved ',filename]);

return
end % function build_stencils

%% make_stencils %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [g_stencil,s_stencil,e_stencil] = make_stencils_61(PARAM, FILEDEF) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make_stencils
%
% Consructs the stencils for data collection.  
%
%  There are 3 stencils
%    g_stencil for the general k_level data files
%    s_stencil for the source for sill overflow data files
%    e_stencil for the entrainment region of sill overflow data files
%  The dimensions of the stencils are set in this script.
%
%  Upon entry, PARAM contains:
%    nx, ny, nzz
%    matrix_dim
% 
% Upon exit:
%    the stencils have been built.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Everything is built in YXZ order (Z varies slowest.)
%
% Build a matrix g_stencil with ones and zeros as a mask for picking
%  up data from the regular k_level IRF function.
% Given the parameters of the half-width, on either side of a potential
% impulse. The parameters for the stencils are hardcoded in this routine.
%
% Y_limit is an internal limit to avoid conflict between two regions,
% specifically, Denmak Straight (DS) and Faroe Banks (FB). 0 is no internal
% limit; >0 is limit index value not to exceed; <0 is limit index value must
% be greater than.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read in the file defining the product areas for the overflows
nx = PARAM.nx; ny = PARAM.ny; nzz = PARAM.nzz;  
matrix_dim = PARAM.matrix_dim;

product_array = zeros(151,5); 
product_filename = ([FILEDEF.product_txt_dir,FILEDEF.product_txt_file]);
product_array = importdata(product_filename,' ',0); % reads in the whole file as an array

disp('building g_STENCIL');
  tic;
  
% Set up parameters for the stencil
PARAM.imin   = 1;
PARAM.imax = nx;
PARAM.jmin = 1;
PARAM.jmax = ny;
PARAM.kmin = 1;
PARAM.kmax = nzz;
%
PARAM.X_band = 2;
PARAM.Y_band = 2;
PARAM.Z_band_neg = 2;
PARAM.Z_band_pos = 2;
PARAM.Y_limit = 0;

block_size = (PARAM.X_band*2+1) * (PARAM.Y_band*2+1) * ...
             (PARAM.Z_band_pos+PARAM.Z_band_neg+1);
max_stencil_entries = nx*ny*nzz*block_size;
g_stencil = spalloc(matrix_dim,matrix_dim,max_stencil_entries); %allocate space

[SY,SX,SZ] = setup_stencil_bands_61(PARAM); % set up the banded matrixes for each dimension

% Make a sparse matrix with the ones arranged according to the total
% POP mesh stencil
 g_stencil = kron(SZ,(kron(SX,SY)));  %Y X Z order (y varies fastest)
  toc;
%   figure (2);
%   spy (g_stencil)
%   title ('g stencil')
%   disp(' spy of g_stencil');
%   keyboard;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%% Build the matrix s_STENCIL for picking up data from the source %%%%
%    regions of sill overflows.

disp('building s_STENCIL');
  tic;

% Table of Source Info from IRF_mod.F90, extended with band definitions.
    ovf_cnt = PARAM.ovf_cnt;
    %                 1    2     3    4  
    X_band_s     = PARAM.X_band_s;
    Y_band_s     = PARAM.Y_band_s;
    Z_band_neg_s = PARAM.Z_band_neg_s;
    Z_band_pos_s = PARAM.Z_band_pos_s;    
%%%%  Initialize the s_STENCIL                                      %%%%

% Determine the maximum number of entries 
max_stencil_entries = 3600; % allowance for production area entries
for ovf_i = 1:ovf_cnt
    block_size = (X_band_s(ovf_i)*2+1) * (Y_band_s(ovf_i)*2+1) ...
              * ((Z_band_neg_s(ovf_i)+Z_band_pos_s(ovf_i)+1));
    max_stencil_entries = max_stencil_entries + block_size;
end % for all overflow regions

% initialize an empty sparce matrix for the stencil
s_stencil = spalloc(matrix_dim,matrix_dim, max_stencil_entries);

%%%% For each overflow region, form mask, and add to stencil       %%%%%   
for ovf_i = 1:ovf_cnt
    
    PARAM.imin       = PARAM.src_imin(ovf_i);
    PARAM.imax       = PARAM.src_imax(ovf_i);
    PARAM.jmin       = PARAM.src_jmin(ovf_i);
    PARAM.jmax       = PARAM.src_jmax(ovf_i);
    PARAM.kmin       = PARAM.src_kmin(ovf_i);
    PARAM.kmax       = PARAM.src_kmax(ovf_i);
    PARAM.X_band     = PARAM.X_band_s(ovf_i);
    PARAM.Y_band     = PARAM.Y_band_s(ovf_i);
    PARAM.Z_band_neg = PARAM.Z_band_neg_s(ovf_i);
    PARAM.Z_band_pos = PARAM.Z_band_pos_s(ovf_i);
    PARAM.Y_limit    = PARAM.Y_limit_s(ovf_i);
    
    [SY,SX,SZ] = setup_stencil_bands_61(PARAM); % set up SX, SY, SZ for one region
  
    ovfl_stencil = kron(SZ, (kron(SX,SY)));
    
    [prod_stencil] = make_prod_stencil_61(PARAM, product_array,ovf_i, ...
                                          ovfl_stencil); % set up product stencil for one region
                                           % (sparse)
    s_stencil = s_stencil | ovfl_stencil | prod_stencil;                   
    
end % for each overflow region

  toc;
%   figure (3);
%   spy (s_stencil);
%   title ('s stencil');
%   disp(' spy of s_stencil');
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Build the e_STENCIL for picking up data from the entrainment   %%%% 
%    regions of sill overflows.
disp('building e_STENCIL');
tic;

X_band_e     = PARAM.X_band_e;
Y_band_e     = PARAM.Y_band_e;
Z_band_neg_e = PARAM.Z_band_neg_e;
Z_band_pos_e = PARAM.Z_band_pos_e;

%%%% Initialize the e_STENCIL                                       %%%%

% Determine the maximum number of entries 
max_stencil_entries = 3600; % allowance for production area entries
for ovf_i = 1:ovf_cnt
    block_size = (X_band_e(ovf_i)*2+1) * (Y_band_e(ovf_i)*2+1) ...
              * ((Z_band_neg_e(ovf_i)+Z_band_pos_e(ovf_i)+1));
    max_stencil_entries = max_stencil_entries + block_size;
end % for all overflow regions

% initialize an empty sparce matrix for the stencil
e_stencil = spalloc(matrix_dim,matrix_dim, max_stencil_entries);

%%%% For each overflow region, form mask, and add to stencil      %%%%%   
for ovf_i = 1:ovf_cnt
    
    PARAM.imin       = PARAM.ent_imin(ovf_i);
    PARAM.imax       = PARAM.ent_imax(ovf_i);
    PARAM.jmin       = PARAM.ent_jmin(ovf_i);
    PARAM.jmax       = PARAM.ent_jmax(ovf_i);
    PARAM.kmin       = PARAM.ent_kmin(ovf_i);
    PARAM.kmax       = PARAM.ent_kmax(ovf_i);
    PARAM.X_band     = PARAM.X_band_e(ovf_i);
    PARAM.Y_band     = PARAM.Y_band_e(ovf_i);
    PARAM.Z_band_neg = PARAM.Z_band_neg_e(ovf_i);
    PARAM.Z_band_pos = PARAM.Z_band_pos_e(ovf_i);
    PARAM.Y_limit    = PARAM.Y_limit_e(ovf_i);
    
    [SY,SX,SZ] = setup_stencil_bands_61(PARAM); % set up SX, SY, SZ for one region
    
    ovfl_stencil = kron(SZ,(kron(SX,SY)));
    
    [prod_stencil] = make_prod_stencil_61(PARAM,product_array,ovf_i, ...
                                       ovfl_stencil); % set up product stencil for one region
                                                      % (sparse)
    e_stencil = e_stencil | ovfl_stencil | prod_stencil;         
    
end % for each overflow region

  toc;
%   figure (4);
%   spy (e_stencil);
%   title('e stencil');
%   disp(' spy of e_stencil');
%   % keyboard;
  
return
end % function make_stencils

%% setup_stencil_bands %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [SY,SX,SZ] = setup_stencil_bands_61(PARAM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  setup_stencil_bands  makes the bands for the POP stencils
%
%  X is cyclic; Y and Z are not.
%
% Upon entry:
%  X_band, Y_band and Z_band_neg and Z_band_pos have been set in PARAM  
%    These parameters are the half-width, on either side of a potential 
%    impulse. X_band and Y_band are symetric, but the Z_bands are not.
%  imin and imax give the minumum and maximum x-indexs for the masks. 
%    For the general case, imin = 1, imax = nx
%  jmin, jmax, kmin, kmax are also set in PARAM
%  nx, ny, and nzz have been set in PARAM
%
% Upon exit:
%  SX, SY, and SZ contain the sparce diagonal matrices for the bands.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ny = PARAM.ny; nx = PARAM.nx; nzz = PARAM.nzz;
imin       =  PARAM.imin;
imax       =  PARAM.imax;
jmin       =  PARAM.jmin;
jmax       =  PARAM.jmax;
kmin       =  PARAM.kmin;
kmax       =  PARAM.kmax;
X_band     =  PARAM.X_band;
Y_band     =  PARAM.Y_band;
Z_band_neg =  PARAM.Z_band_neg;
Z_band_pos =  PARAM.Z_band_pos; 
Y_limit    =  PARAM.Y_limit; 


%%% Make nx by nx matrix, + and - X_band wide, on the diagonal %%%%
EX = zeros(nx,(2*X_band+1));
EX(imin:imax,:) = ones((imax-imin+1),(2*X_band+1));
SX = spdiags (EX, [-X_band:1:X_band], nx,nx);

% wrap around for circular geometry if block is at either end  or both of x
% dimension.  Wrap around the diagonal - the impulses are on the diagonal
% and collecting data on the columns. 
if ((imin - X_band) < 1) % if wrap needed from the low end
    for xcyc = 1:(X_band-imin + 1)
        ilim = min((X_band-xcyc+1),imax);
        SX(nx-xcyc+1,(imin:ilim)) = ones(1,(ilim-imin+1));
    end % xcyc
end % if wrap needed

if ((imax + X_band) > nx) % if wrap needed from the high end
    for xcyc = 1:(imax+X_band-nx)
        SX(xcyc,((nx-X_band+xcyc):nx)) = ones(1,(X_band-xcyc+1));
    end % xcyc
end %if wrap needed

%figure (50);
% spy (SX);
% title ('SX');
% disp('keyboard to examine SX');
% keyboard;
    
%%%% Make ny by ny matrix, + and - Y_band wide, on the diagonal %%%%
EY = zeros(ny,(2*Y_band+1));
EY(jmin:jmax,:) = ones;
SY = spdiags (EY, [-Y_band:1:Y_band], ny,ny);

% If there is an internal Y limit, zero out any space beyond the limit. 
% (This takes care of the potential overlap between DS and FB using this
%  simple stencil building scheme.)
if (Y_limit > 0) % pos, nz, cannot exceed
    Y_lim = Y_limit+1;
    SY(Y_lim:ny,:) = zeros((ny-Y_lim+1),ny);
    SY(:,Y_lim:ny) = zeros(ny,(ny-Y_lim+1));
elseif (Y_limit < 0) % neg, nz, cannot be less than
    Y_lim = abs(Y_limit)-1;
    SY(1:Y_lim,:) = zeros(Y_lim,ny);
    SY(:,1:Y_lim) = zeros(ny,Y_lim);
end % if there is an internal Y limit

%%%% Make nzz by nzz matrix, + and - Z_band wide, on the diagonal %%%%
EZ = zeros(nzz,(Z_band_neg+Z_band_pos+1));
EZ(kmin:kmax,:) = ones;
SZ = spdiags (EZ, [-Z_band_neg:1:Z_band_pos], nzz,nzz);

return
end %function setup_stencil_bands

%% make_prod_stencil %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [prod_stencil] = make_prod_stencil_61(PARAM,product,ovf_i, ...
                          ovfl_stencil)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  make_prod_stencil  makes the extension mask for the product area for the
%     region.
%
%  X is cyclic; Y and Z are not.
%
% Upon entry:
%  product has the input data for all the product areas.
%  ovfl_stencil has the mask for the region without the product area
%  extension.
%  imin and imax give the minumum and maximum x-indexs for the masks. 
%    For the general case, imin = 1, imax = nx
%  jmin, jmax, kmin, kmax are also set.
%  nx, ny, and nz have been set.
%
% Upon exit:
%  prod_stencil is the sparse matrix with the mask for the particular
%  region.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ny = PARAM.ny; nx = PARAM.nx; nzz = PARAM.nzz;
matrix_dim =  PARAM.matrix_dim;
ovf_cnt    =  PARAM.ovf_cnt;
% Initiate the 3_D product areas mask and list of columns and rows for the
% sparse matrix
prod_3D = zeros(ny,nx,nzz);
icol = [ ]; irow = [ ];

% Find the data for the region in the product data
  % Second line has the number of overflow regions - better match
  ck_num_ovfl = product(2,1);          % ignore the first line
  if (ck_num_ovfl ~= ovf_cnt)
      disp('overflow count in product defs does not match');
      keyboard
  end % if bad number of overflow regions
  prod_i = 3; % initiate the index into the rows of the product array.
  
  % Next line has the overflow region index
  overfl_i = product(prod_i,1);
  
  prod_i = prod_i + 1;       %next line
  while (overfl_i ~= ovf_i)  % find the overflow region of interest
      num_prod_sets = product(prod_i,1); % skip thru the region not of interest
      prod_i = prod_i + 1;
      for set_i = 1:num_prod_sets 
          num_set_entries = product(prod_i,1); % number of entries in set
          prod_i = prod_i + num_set_entries + 1; % skip the entries
      end % for all sets of entries
      overfl_i = product(prod_i,1); % next overflow region
      prod_i = prod_i + 1;       %next line
  end % while not the region of interest
  
% Found the region of interest - set up the mask for the product areas
  num_prod_sets = product(prod_i,1); % get number of product sets
  prod_i = prod_i + 1;
  for set_i = 1:num_prod_sets % for each product set
      num_set_entries = product(prod_i,1); % number of entries in set
      prod_i = prod_i + 1;
      for entry_i = 1:num_set_entries  % for each entry in set
          i_index = product(prod_i,2);
          j_index = product(prod_i,3);
          k_index = product(prod_i,4);
          orient  = product(prod_i,5); % orientation of wall face
          % we want the cell next to the one designated, according to
          % the orientation given
          if (orient == 1)
              i_index = i_index + 1;
          elseif (orient == 2)
              j_index = j_index + 1;
          elseif (orient == 3)
              i_index = i_index - 1;
          elseif (orient == 4)
              j_index = j_index - 1;
          else
              disp ('unknown orientation for product cell');
              keyboard
          end % orientation

          prod_3D(j_index, i_index, k_index) = 1; % set the point in 
                                                  % the 3D mask   
          prod_i = prod_i + 1;
      end % for each entry
  end % for the product set

% make a sparce matrix with the mask for each of the impulse locations in
% the overflow matrix
  
% make a list of indexes of impulse locations for overflow stencil
  m_eye = speye(matrix_dim);
  nz_imp_m = ovfl_stencil .* m_eye; % impulse locs are on the diagonal
  nz_imp_c = nz_imp_m * (ones(matrix_dim,1)); % make impulse locs into a column vector
  nz_imp = find(nz_imp_c); % get the index of the nz impulse locs
  l_nz_imp = length(nz_imp);
  
% make a list of indexes of product mask locations
  nz_prod = find(prod_3D(:));
  l_nz_prod   = length(nz_prod);
  
  for col_i = 1:l_nz_imp  % for each impulse location 
      clear imp_col;
      imp_col(1:l_nz_prod) = [nz_imp(col_i)]; % set current column index for all prod indexes
      
      icol = [icol; imp_col(:)];
      irow = [irow; nz_prod(:)];
  end % for each impulse location
  ival = ones(length(icol),1);
  
  prod_stencil = sparse(irow,icol,ival,matrix_dim,matrix_dim);
return
end % function make_orod_stencil

%% FORMAT for the PRODUCT region definition file %%%%%%%%
% The definitions are extracted from the overflow definition file that
% comes with the system release : gx1v6_overflow equivalent.   The fields
% in the system-released document are well documented with comments.  Below
% is shown the location of the fields that have been extracted to make a
% file that can be directly read into an array to build the product region
% part of the mask.  
%
% Column 1 is in column 1 of the edited text file.  
% Subsequent columns are blank delimited.
% Line 1 is the first line of the text file.  
% There are to be no comments, or this will not work with this version of
% the code.
% No blank lines. (The blank lines below are to accomodate the
% definitions.)
% All the data except for the first line is constructed by text-editing the
% 'gx1v6_overflow' file.
% 
% 
% 1 2 3 4 5        % Line 1 is added so that the file will be read as
%                     having 5 columne. 
% 4                % Number of overflow regions
% 1                % Overflow region ID number
% m                % Number of subregions
% n                % Number of entries defining a subregion ( 1 entry per
%                       line )
% 1 2 3 4 5        % Entry: Column 1 is the entry number (added)
%                           Column 2 is the i index
%                           Column 3 is the j index
%                           Column 4 is the k index
%                           Column 5 is the orientation
% 1 2 3 4 5        % Entry ... for n entry lines
%  ...
% n                % number of entries in the next subregion
% 1 2 3 4 5        % Entry ... for n entry line
% ...
% n                % repeated sets for m subregions
% ...
% 2                % Overlow region ID number
% ...              % repeated format for each overflow region
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  30  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




