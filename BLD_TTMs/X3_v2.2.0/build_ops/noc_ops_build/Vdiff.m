function [Dv] = Vdiff(MET, GEO_PARAM,KTF, KBF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Builds the sparce diagonal array of operators for vertical diffusion 
%  from the 3-D arrays extracted from the IRF-run history file.
%
% Inputs:
%   KTF           % Kappa coefficient (3D), value at the top of the grid 
%                 % cell, in m^-2 s^-1, extracted from the history file
%                 % KTF  = KPP1 + GM
%
%   KBF           % kappa coefficient (3D), vlaue at the BOTTOM of the grid
%                 % 
%
%   NPP_NLK (obs) % coefficient(3D)for surface flux, in m^-1
%      
%   GEO_PARAMS: 
%      ny,nx,nzz  % dimensions
%      matrix_dim % the sparse array dimension
%      
%   MET:
%      MASK
%      DZT        % layer thickness (dz in POP documentation)
%      DZU        % vertical distance from the middle of one layer 
%                 % to the middle of the next one (dzw in POP
%                 % documentation)
% Outputs:
%      Dv         % vertical diffusion operator, a sparce matrix 
%                 % based on KTF. The units are s^-1. The matrix is 
%                 % full sized, including land points.
%
%      DvS        % vertical diffusion surface flux operator, a sparce 
%                 % matrix of KPP_NONLOC_KERN coefficients, in 
%                 % units of m^-1.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% From the POP Reference Manual (revision of March 2010),eq. 3.19, the vertical 
% diffusion operator, Dv, in spatially discretized form, is given by:
%
% Dv(C) = 1/dz(k) * ( (Kappa(k-1/2) * (C(k-1) - C(k))/dz(k-1/2)) ...
%                     - kappa(k+1/2) * (C(k) - C(k+1))/dz(k+1/2) ).
%  where C is the tracer concentration.
% 
% Kappa is maintained in the model as the value at the top face of the 
% T-grid cell,  = kappa(k-1/2), which we will denote by KTF.  
% Similarly kappa(k+1/2) can be denoted by KBF.
% The dz(k-1/2) is equivalent to dzw(k), the distance between the 
% center point of the kth grid cell and the grid cell below it. 
% If the recipricals of the depth dimensions are denoted as R_dz and R_dzw,
% then the equation becomes:
%
% Dv(C) = R_dz(k) * ( R_dzw(k-1) * ( (KTF(k) * C(k-1) ...
%         - C(k)) - R_dzw(k) * KBF(k) * (C(k) - C(k+1)) ) ). 
% 
% Segregate the terms that operate on C(k), C(k-1), and C(k+1), such that
%
%  Dv = Dvm1(C(k-1)) + Dv0(C(k)) + Dvp1(C(k+1))
% 
%  and the coefficients are spatially related to the location of the 
%  concentration, gives the terms
%
%    Dvm1 =R_dz_kp1(k-1) * R_dzw_kp1(k-1) * KBF(k-1)
%
%    Dv0 = - R_dz(k) *( R_dzw(k)* KTF(k) + R_dzw_kp1(k) * KBF(k))
%
%    Dvp1 = R_dz_km1(k+1) * R_dzw(k+1) * KTF(k+1)
%
%  where R_dz_kp1(k-1) is R_dz(k).  
% 
%  revised for variable thickness
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ny = GEO_PARAM.ny;  nx = GEO_PARAM.nx; nzz = GEO_PARAM.nzz; nz = nzz-1;
  ocn_mask   = MET.MASK;
  dzw        = MET.DZU;
  dz         = MET.DZT;
  matrix_dim = GEO_PARAM.matrix_dim;
  % Set up shifts up and down a column (nz dimension)
  shftkm1 = [1,1:(nzz-1)]; % shift (contents) from k-1 to k
  shftkp1 = [2:nzz,nzz];  % shift (contents) from k+1 to k 
    
  %  Calculate reciprocals R_dzw = 1/dzw and R_dz and associates.
  %
  %  Note that contrary to documentation, dzw(k=1) = 1/2 * dz(k=1),
  %  and dzw(k=2) is the distance from the middle of grid box 1 to
  %  the middle of grid box 2.
   
  R_dzw = 1 ./ dzw;
  R_dzw_km1 = R_dzw(:,:,shftkm1); % at k - 1
  R_dzw_km1(:,:,1) = 0;
  R_dzw_kp1 = R_dzw(:,:,shftkp1); % at k + 1
  R_dzw_kp1(:,:,nz) = 2./dz(:,:,nz);  % dzw = 1/2 dz

  R_dz = 1 ./ dz;
  R_dz_km1 = R_dz(:,:,shftkm1);
  R_dz_km1(:,:,1) = 0;
  R_dz_kp1 = R_dz(:,:,shftkp1);
  R_dz_kp1(:,:,nz) = 0;  

  % Calculate each term for the diffusion operator %%

  Dvm1 = zeros(ny,nx,nzz); % dimension the term
  Dv0  = zeros(ny,nx,nzz);
  Dvp1 = zeros(ny,nx,nzz);
 
  for k = 1:nz % leave nzz level = 0

    Dvm1(:,:,k) = KBF(:,:,k) .* R_dzw_kp1(:,:,k) .* R_dz_kp1(:,:,k);
    
    Dv0(:,:,k) = - R_dz(:,:,k) .* (R_dzw(:,:,k) .* KTF(:,:,k) + ... 
       R_dzw_kp1(:,:,k) .* KBF(:,:,k));

    Dvp1(:,:,k) = R_dz_km1(:,:,k) .* R_dzw(:,:,k) .* KTF(:,:,k);

  end % for k = 1:nz

  % Make into a sparce diagonal array in (j,i,k) order
  xy_dim = nx*ny; % spacing for z terms

  Dv = spdiags( [Dvm1(:) Dv0(:) Dvp1(:)], ... 
       [-xy_dim 0 +xy_dim], matrix_dim, matrix_dim);
   
  return
end % function Vdiff
 %%%%%%%%%%%%%%%%%%%%%%%%% -30-%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

