function [] = make_A_periodic(F,P,MET)
%
% Calculates the 2-D modification to the surface grid cells for the advection 
% operator, to make it cyclic on an annual basis, without divergence.
% Checks operators for mass balance and divergence.
% Outputs operator .mat files to "final" directory, by month plus an annual
% set, as period (0).
%
% Refer to "An offline emplicit solver for simulating prebomb radiocarbon"
% by Ann Bardin, Francois Primeau, Keith Lindsay, in Ocean Modelling, 
% in press, for a full mathematical explaination of the adjustment terms
% calculated here. 
% 
% INPUTS
%    MET                     % geometric grid dimensioning parameters
%    F                       % file locations
%     ops_dir                % final operator output directory
%    P                       % control parameters
%     first_month_in_annual % first month to be included in the annual
%                            % average - provides for the annual average
%                            % to include months from neighboring years.
%      months_per_annual     % total number of months in the multi-year avg
%                            % different resolution in the months of
%                            % a year. 
% OUTPUTS
%    No variables returned.
%    Operator .mat files are saved in the "final" directory, by calendar 
%    month, plus an annual set, with file name for month (0).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp('making A periodic on multi-annual cycle');
    [ny,nx,nz] = size(MET.MASK);
    iocn = find(MET.MASK > 0);           % ocn cells
    V = d0(MET.VOL(iocn));               % volume, iocn-sized
   
  %% Do the annual average adjustment 
  
  % load in the annual average transport operators, iocn size 
   
    ann_ops_fn = ([F.noc_dir,'MTM0.mat']);
    eval(['load ', ann_ops_fn ]); % [A,H,D,T,num_days]
  
  % compute the surface  mass convergence assumed due to the
  % changing sea surface height, precipitation, evaporation, river inputs
    disp('eliminating the divergence in the surface layer of A')
    rho = A*(0*iocn+1);                % divergence
    A2d = mkAtilde2d(rho,MET);         % adjustment to eliminate the divergence    
                                       % in the surface layer
  % Apply the correction to A and T
    A = A + A2d;
    T = T + A2d; 
  % Display the resulting operator divergence
    % [P] = display_op_div(A,H,D,T,P,MET);
  
  % Save completed operator files; this is the start of the completed files
   % disp(['saving ',ann_ops_fn]);
    dxidt = 0*iocn;                        % zeros for annual average ops
    out_file_name = ([F.ops_dir,'MTM0.mat']);
    eval(['save ',out_file_name,' A H D T dxidt num_days -v7.3;']);  
    disp([ out_file_name,' saved']);

  % run mass-balance checks
    disp(['running mass balance and divergence checks']);
    C_period = (['0']);
    [P] = ck_ops_massb_div(F,P,MET,C_period); 
  
    disp('completed making A annually periodic');

  %%  monthly adjustments
  %  A2d, annual adjustment to make A annually periodic has already been calculated
  %  Need to calculate an additional SSH-change adjustment based on the 
  %  monthly Advection operator.

  for period = 1:F.num_periods      % for each period to be created
     fperiod_num = P.fperiod + period -1;
     lperiod_num = fperiod_num + (F.nyears-1).*F.num_periods;
     % init totals for averaging
     iocn = MET.iocn;
     iocn_matrix_dim = length(iocn);

     % Initialize the totals matrixes
     col_zeros   = (zeros(iocn_matrix_dim,1));
     Adv_total   = spdiags(col_zeros, 0, iocn_matrix_dim, iocn_matrix_dim);
     Hdiff_total = spdiags(col_zeros, 0, iocn_matrix_dim, iocn_matrix_dim);
     Dv_total    = spdiags(col_zeros, 0, iocn_matrix_dim, iocn_matrix_dim);
     T_op_total  = spdiags(col_zeros, 0, iocn_matrix_dim, iocn_matrix_dim);

     for period_num = fperiod_num:F.num_periods:lperiod_num
         op_file = (['MTM',int2str(period_num),'.mat']);
         filename = [F.noc_dir, op_file];
         load (filename); %  A H D T num_days
         % days_total = days_total + num_days;
         T_op_total = T_op_total + T; 
         Adv_total = Adv_total + A;
         Hdiff_total = Hdiff_total + H;
         Dv_total = Dv_total + D;
     end % for month in multi-annual average
     T     = T_op_total  ./ F.nyears;
     D     = Dv_total    ./ F.nyears;
     H     = Hdiff_total ./ F.nyears;
     A     = Adv_total   ./ F.nyears; 
  
     % apply the annual adjustment to the seasonal operators
     A = A+A2d;
     T = T+A2d;
     Ms = MET.MASK;                 % isolate layer 1 for SSH-change adjustment
     Ms(:,:,2:end) =0;
     is = find(Ms(iocn));
     dV = MET.VOL(iocn);
     dA = MET.TAREA(iocn);
     rho = (A)*(0*iocn+1);          % rho has the adjustment for the surface 
     dxidt = 0*iocn;
     % version j:
     % mo_corr = 0*iocn;
     % mo_corr(is) = -( rho(is) -(sum(rho(is).*dV(is) )/sum(dV(is)) ) );
     % A = A + d0(mo_corr);
     % T = T + d0(mo_corr);
     % version k
     % mo_corr = 0*iocn;
     % mo_corr(is) = -( rho(is) );
     % A = A + d0(mo_corr);
     % T = T + d0(mo_corr);
     % version c:
     % dxidt(is) = -( rho(is) -(sum(rho(is).*dV(is) )/sum(dV(is)) ) )./dV(is);
     % version d:
     dxidt(is) = -( rho(is) -(sum(rho(is).*dV(is) )/sum(dV(is)) ) )./dA(is);
     % version e:
     % dxidt(is) = -( rho(is) -(sum(rho(is).*dV(is) )/sum(dV(is)) ) );
     % version g:
     % dxidt(is) = -rho(is)./dA(is);
     % version h:
     % dxidt(is) = -rho(is);
       
     % Save completed operator files;   
     out_file_name = ([F.ops_dir,'MTM',int2str(period),'.mat']);
     eval(['save ',out_file_name,' A H D T dxidt num_days -v7.3;']);
     disp(['saved ',out_file_name])  
   
     % run mass-balance checks 
     disp(['running mass balance and divergence checks']);
     C_period = int2str(period);
     [P] = ck_ops_massb_div(F,P,MET, C_period); 

     disp(['completed making adjusted operators for period ',int2str(period)]);
  end % for each period
  
end % function make_A_periodic

%% mkAtilde2d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Atilde] = mkAtilde2d(rho,MET)
% Makes the periodic-adjustment operator, used to make Advection 
% periodic on an annual basis.
%
% INPUTS:
%    MET                     % geometric grid dimensioning parameters  
%    rho                     % row divergence vector, for A (advection)
%                            % iocn sized
%
% OUTPUTS
%    Atilde                  % periodic adjustment operator for Advection
%                            % iocn sized
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [GRAD,GRD] = grad(MET);            % 3-D gradient operator GRADD, with the 
                                       % individual direction components in
                                       % the structure GRD
    [DIV,DV]   = div(MET);             % 3-D divergence operator DIV, with the
                                       % individual direction components in
                                       % the structure DV
                                       
    [ny,nx,nz] = size(MET.MASK);
    iocn       = find(MET.MASK(:));

    Ms = MET.MASK;                     % surface only grid cells 
    Ms(:,:,2:end) = 0;
    is = find(Ms(iocn)==1);

  % compute the velocity correction potential function
  % ***SYMMETRIC*** volume weighted Laplacian operator
  %
  
    DX = DV.DX;
    GX = GRD.GX;
    WX = DV.WX*GRD.WX;
    D2X = DX*(WX*GX);
    D2X = [D2X(:,iocn)]'; D2X = [D2X(:,iocn)]';
  %
    DY = DV.DY;
    GY = GRD.GY;
    WY = DV.WY*GRD.WY;
    D2Y = DY*(WY*GY);
    D2Y = [D2Y(:,iocn)]'; D2Y = [D2Y(:,iocn)]';
  %
    VDEL2 = D2X+D2Y;
  
  % remove the null space by eliminating one unknown and one equation
  % (i.e. set the gauge to zero)
    VDEL2 = VDEL2(is,:); VDEL2 = VDEL2(:,is);
    xVDEL2 = VDEL2(1:end-1,:); xVDEL2 = xVDEL2(:,1:end-1);
  
    fprintf('factoring the horizontal laplacian...');
    xrho = MET.VOL(iocn).*rho;
    xrho = xrho(is);
  % make sure that rho sums to zero
    xrho = xrho-mean(xrho);
    xrho = xrho(1:end-1);
    tic
      FxVDEL2 = mfactor(xVDEL2);
    toc
    phi = 0*MET.MASK;
  
    phi(iocn(is(1:end-1))) = mfactor(FxVDEL2,xrho);
    phi(iocn(is(end))) = 0;
  %
    NDX = zeros(ny,nx,nz);             % make index for all cells
    n = ny*nx*nz;
    NDX(:) = 1:n;
  
  % permute the rows of the identity matrix to get the neighbor operators
    I = speye(ny*nx*nz);
    in = NDX([2:ny,ny],:,:);   IN = I(:,in(:))';
    ie = NDX(:,[2:nx,1],:);    IE = I(:,ie(:))';
    it = NDX(:,:,[1,1:nz-1]);  IT = I(:,it(:))';
  
  % create operators to average onto the eastern, northern, and top
  % gridbox faces
    Ae = (IE+I)/2;
    An = (IN+I)/2;
    At = (IT+I)/2;
  
  % compute the velocity corrections
  % making sure there is no flow through the basin boundaries
    utilde = zeros(ny,nx,nz);   
    vtilde = zeros(ny,nx,nz);   
    wtilde = zeros(ny,nx,nz);
    utilde(:) = GRD.WX*(GRD.GX*phi(:));  
    utilde(:,:,2:end) = 0;
    vtilde(:) = GRD.WY*(GRD.GY*phi(:));  
    vtilde(:,:,2:end) = 0;
  %
    Atilde = -DIV*[d0(utilde(:))*Ae;...
                   d0(vtilde(:))*An;...
                   d0(wtilde(:))*At];  % wtilde is zeros
  % resize to iocn           
    Atilde = [Atilde(:,iocn)]'; Atilde = [Atilde(:,iocn)]';

    return
end % function make Atilde 2D

%% grad %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [GRAD,varargout] = grad(MET,alpha)
% Makes the 3-D gradient operator, GRAD
%
% INPUTS: 
%    MET                     % geometric definitions
%    alpha                   % weighting coefficient for horizontal vs.
%                            % vertical directions; default = 1. 
%
% OUTPUTS:
%     GRAD                   % 3-D gradient operator, with full-sized
%                            % dimensions taken from MET
%                            % stacked in X,Y,Z order
%     varargout (optional)   % if a second output variable is supplied,
%                            % a structure is returned with components of 
%                            % the total gradient, GX, GY,and GZ,
%                            % and WX, WY ,WZ.  GX, GY, and GZ give the
%                            % difference between the appropriate
%                            % neighboring cell and the subject cell for each
%                            % direction (ex: IE-I for the x-direction);
%                            % and WX, WY, and WZ give weighed 1/distance  
%                            % between their centers (ex:alpha*1/DZU for the
%                            % x-direction).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  if (nargin<2)
      alpha = 1;
  end
  
  % extract the metric coefficients

  M = MET.MASK;
  DXU = MET.DXU;             % horizontal length of the North face of the U_cell
  DYU = MET.DYU;             % horizontal length of the East face of  the U-cell
  DZU = MET.DZU;             % vertical thickness of the U_cell 
  [ny,nx,nz] = size(M);
  NDX = zeros(ny,nx,nz);     % 1-D index for all the cells
  NDX(:) = 1:ny*nx*nz;

  % permute the rows of the identity matrix to get the neighbor operators

  I = speye(ny*nx*nz);
  in = NDX([2:ny,ny],:,:);   IN = I(:,in(:))';
  ie = NDX(:,[2:nx,1],:);    IE = I(:,ie(:))';
  it = NDX(:,:,[1,1:nz-1]);  IT = I(:,it(:))';

  % grad operator with Neuman b.c. at the basin boundaries

  WX = d0(M(:).*(IE*M(:))./DXU(:)); GX = (IE-I);
  WY = d0(M(:).*(IN*M(:))./DYU(:)); GY = (IN-I);
  WZ = d0(M(:).*(IT*M(:))./DZU(:)); GZ = (IT-I);
 
  GRAD = [alpha*WX*GX;alpha*WY*GY;WZ*GZ];      % assemble grad operator

  if (nargout>1)
      % supply the components of the grad operator  
      GRD.GX = GX; GRD.GY = GY; GRD.GZ = GZ;
      GRD.WX = alpha*WX; GRD.WY = alpha*WY; GRD.WZ = WZ;
      varargout{1} = GRD;
  end
  return
end % function grad 

%% div %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DIV,varargout] = div(MET)
% make the divergence operator
%
% INPUTS: 
%    MET                     % geometric definitions
%
% OUTPUTS:
%     DIV                    % 3-D divergence operator, with full-sized
%                            % dimensions taken from MET
%                            % laid out in X,Y,Z order
%     varargout (optional)   % if a second output variable is supplied,
%                            % a structure is returned with components of 
%                            % the total divergence.
%                            % IV is the inverse of the volume.
%                            % DX, DY,and DZ,
%                            %  give the
%                            % difference matrix between the subject cell and 
%                            % the appropriate cell for each
%                            % direction (ex: I-IW for the x-direction).
%                            % WX, WY, and WZ give 1/(area of the 
%                            % perpendicular face of the grid cell),
%                            % (ex: 1/EFACE for the
%                            % x-direction).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % extract the metric coefficients
  M = MET.MASK;
  EFACE = MET.EFACE;
  NFACE = MET.NFACE;
  TFACE = MET.TFACE;
  Vol = MET.VOL;
  %
  [ny,nx,nz] = size(M);
  NDX = zeros(ny,nx,nz);
  NDX(:) = 1:ny*nx*nz;
  
  % permute the rows of the identity matrix to get the neighbor operators
  I = speye(ny*nx*nz);
  is = NDX([1,1:ny-1],:,:); IS = I(:,is)';
  iw = NDX(:,[nx,1:nx-1],:); IW = I(:,iw)';
  ib = NDX(:,:,[2:nz,nz]); IB = I(:,ib)';

  % DIV operator 
  DIV = d0(1./Vol(:))*[(I-IW)*d0(EFACE(:)),...
                       (I-IS)*d0(NFACE(:)),...
                       (I-IB)*d0(TFACE(:))];
  DV.IV = d0(1./Vol(:));
  DV.DX = I-IW; 
  DV.DY = I-IS;
  DV.DZ = I-IB;
  DV.WX = d0(EFACE(:));
  DV.WY = d0(NFACE(:));
  DV.WZ = d0(TFACE(:));
  if (nargout>1)
    varargout{1} = DV;
  end
  return
end %function div

%% load_iocn_ops %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,H,D,T, num_days] = load_iocn_ops(F, month)
% loads A, H, D, T iocn-sized operators for the specified month
% from the tmp_ops_dir.  Stops on a keyboard if the file DNE.
% Month 0 = annual.
%
% INPUTS:
%    F                       % structure defining file locations
%      tmp_ops_dir           % directory to from which to load
%    month                   % month number for the set to load
%
% OUTPUTS;
%    A,H,D,T operators       % iocn sized
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  C_mo = int2str(month);
  fn = ([F.ops_dir,'MTM',C_mo,'.mat']);
  if exist (fn, 'file')
     load(fn);  % A H D T num_days
  else 
     disp(['file DNE:',fn]);
     disp('in load_iocn_ops / make_A_periodic');
     keyboard
  end % if file exists
  return
end % function load_iocn_ops

%% display_op_div %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P] = display_op_div(A,H,D,T,P,MET)
% calculates and makes displays for operator divergence 
%
% INPUTS
%    A,H,D,T                 % operators, iocn sized
%    P.fignum                % control parameters, first figure number to use
%    MET                     % geometric grid dimensioning parameters
%
% OUTPUTS
%    P.fignum                % next figure number to use
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    iocn = find(MET.MASK == 1);
    fignum = P.fignum + 1;
    rhoA = (A)*(0*iocn+1);
    figure(fignum)
    plot(rhoA)
    title('rho (divergence) of A')

    fignum = fignum+1;
    rhoT = (T)*(0*iocn+1);
    figure(fignum)
    plot(rhoT)
    title('rho (divergence) of T')

    fignum = fignum+1;
    rhoH = (H) *(0*iocn+1);
    figure(fignum)
    plot(rhoH)
    title('rho of H')

    fignum = fignum+1;
    rhoD = (D) *(0*iocn+1);
    figure(fignum)
    plot(rhoD)
    title('rho of D')

    P.fignum = fignum+1;
    return
end % function display_op_div

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
