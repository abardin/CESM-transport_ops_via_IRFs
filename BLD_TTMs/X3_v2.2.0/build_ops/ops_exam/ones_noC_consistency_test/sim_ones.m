function varargout = sim_ones(varargin)
% Routines to integrate the periodic tracer
% This version has been adapted to run as part of a ones consistency test
%
%   1. For UCI adaptation, set up the global gbec to transfer parameters,
%      file names, and the initial value vectors that are provided in the
%      driver script.
%
%   2. Run sim_ones('test_ones') to integrate the
%      system for one year and return the result at the end of the year
%      in gbec.
%
  
% Original. 15 Sep 2014. AMB ambradl@sandia.gov
% Modified for C14 Dec 2015 AMB abardin@uci.edu
% Modified to use mfactor.m instead of C++ routines Aug. 2021 AMB

% Inputs are via global gbec:
%  gbec.P    parameters for the simulation

% Output is set of per period values, in gbec.P.out_dir.
% Log has the important data summary values.

  [varargout{1:nargout}] = feval(varargin{:});
end % function sim_ones

function ee = get_env ()
% Set up the data input/output and code environments. This routine is meant
% to be edited to match the system.  Adapted to get the parameters 
% passed in global gbec.
global gbec
persistent e
  % Default number of OpenMP threads to use. Set it to the number of cores
  % available on your system.
  e.nthreads =16;                     % not used in the mfactor version
  if ~isfield(e,'dir_data')           % only get this data one time
     P = gbec.P;  
     % Location of operator .mat files:
     e.dir_data = P.ops_dir; 
                                   
     % Name of the .mat file containing the year-level matrices.
     %e.mat_annual = [P.ops_dir,'MTM0'];

    % Names of the periodic operator .mat files.
    e.mat_period = {};
    mat_dir = gbec.P.ops_dir;
    for period = gbec.P.first_period:gbec.P.period_span:gbec.P.last_period    
       e.mat_period{period} = sprintf('%sMTM%d',mat_dir, period);
    end
   
    % Location to save data.
    e.dir_save = gbec.P.out_dir;
  end % if need to read parameters
  ee = e;
  return
end % function ee

function p = set_parms ()
  global gbec
  p = gbec.P;
  e = get_env();
  % Set save_every to 0 to turn off file I/O.
  p.save_every = 0;
  p.verbosity = 0;
end % set_parms

function test_ones()
% run multi-year simulation for concentration of ones
  global gbec
  p = set_parms();
  x0  = gbec.x0;               % input and output through global gbec
  x_std = x0;

  eval(['load ',p.ops_dir,'MET.mat MET ']);
  iocn = MET.iocn;
  vol  = MET.VOL(iocn);
  total_vol = sum(vol);

  for yr = 1:gbec.P.nyears
      t0 = 0;  
      for period = p.first_period:p.period_span:p.last_period
          [x0, p] = e_period_integrate(p, x0, t0, period);
  
          % save period output
          dx = x0 - x_std;
          dxV = (dx.*vol)./total_vol;
          sum_dxV = sum(dx.*vol)/total_vol;
          fn = ([p.out_dir,'dx_mo',int2str(period),'_yr',int2str(yr),'.mat']);
          eval(['save ',fn,' dx dxV sum_dxV']);
         
          disp(['sum dxV = ',num2str(sum_dxV,4),', period ',int2str(period), ...
                ' yr',int2str(yr)]);
      end % for each period 	    
  end % for each year
  
  gbec.x1 = x0;
end % function test_ones

function [x1, p] = e_period_integrate (p, x0, t0, period)
  [n, nc] = size(x0);
  e   = get_env();       % structure with filenames
  xi0 = 0.*x0; 
  [C]     = deal(zeros(n, 4) );
  C(:,1)  = x0; 
  xi      = deal(zeros(n, 2) );
  xi(:,1) = xi0;
  xip1 = 1 + xi(:,1);
  [w, w1, w2] = deal(zeros(n, 1) );
  
  eval(['load ', e.mat_period{period}, ' A H D dxidt num_days ']);  
  tsperDay = 24/p.timestep_hrs; 
  steps_per_period = tsperDay.*num_days;
  dt = p.timestep_hrs*60*60;
  xi(:,2) = xi(:,1) + dt*dxidt; 
% 
  LHS  = d0(1+xi(:,2)) + dt.*(-D);   % LHS on all tracers
  FLHS = mfactor(LHS);
%
  for it = 1:steps_per_period
    k = 4; if (it <= 2), k = it+1; end    

    % RHS
    % source/sink term   
    S = 0*x0;
     
    % Advection Adams-Bashforth method
    switch (it)
      case 1
        w1(:,1) = C(:,k-1); 
      case 2
        w1(:,1) = (1/2).*(3*C(:,k-1) - C(:,k-2) );
      otherwise
        w1(:,1) = (1/12).*(23*C(:,k-1) - 16*C(:,k-2) + 5*C(:,k-3) );
    end	% switch		 

    RHS = xip1.*C(:,k-1) + dt.*(A*w1(:,1) + H*C(:,k-1)+ S);
 
    % Solve for C: (I - dt.*D)*C(:,k) = RHS
    C(:,k) = mfactor(FLHS,RHS); 
    
    % Update time steps.
    if (it > 2)
      C(:,1) = C(:,2); C(:,2) = C(:,3); C(:,3) = C(:,4);
    end
  end % step in period

  x1  = squeeze(C(:,k) ); % return value
  return
end % period integrate


function [m] = load_mesh (p) % not used in the mfactor version
  e = get_env();
  m.mask = p.MASK;    
  
  % Indices into the cube of the volume ...
  m.is_v = find(m.mask);
  % ... and surface elements.
  m.is_s = find(m.mask(:,:,1));
  m.l_isurf = length(m.is_s);
  
  % Collect indices for each water column. The diffusion equations in sm.D
  % and sm.PFD decouple by water column.
  m.vwc.i = zeros(1, numel(m.is_v));
  m.vwc.p = zeros(1, numel(m.is_s)+1);
  Is = nan(size(m.mask));
  Is(m.is_v) = 1:numel(m.is_v);
  for (k = 1 : numel(m.vwc.p)-1)
    [i j] = ind2sub(size(m.mask(:,:,1)), m.is_s(k));
    is = Is(i,j,:);
    is = squeeze(is(~isnan(is)));
    m.vwc.i(m.vwc.p(k) + (1:numel(is))) = is - 1;
    m.vwc.p(k+1) = m.vwc.p(k) + numel(is);
  end %k
end % load_mesh
% ------------------------------------------------------------------------------
% Utils.

function pr (varargin)
  fprintf(1, varargin{:});
end
function s = rmsp (s)
  s(s == ' ') = '';
end

