function [ good ] = pulse_ck(P,MET, F)
%
% Check that each ocean point received one and only one pulse.
% Displays a message, and stops on a keyboard if there are missing or
% too many pulses.
%
% INPUTS:
%    MET                     % geometric grid dimensional parameters
%    P                       % control parameters
%    F                       % file location specification
%      tmp_ops_dir           % directory with the pulse matrices, saved
%                            % from the IRF run history file
% OUTPUT:
%    good                    % true false setting for successful check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  [ny,nx,nz] = size(MET.MASK);
  KMT = MET.KMT;
  good = true; 
  ops_dir = F.tmp_ops_dir;
  load ([ops_dir,'Pulse_k1.mat']);
  Pulse_k1 = Pulse_kx;
  load ([ops_dir,'Pulse_k2.mat']);
  Pulse_k2 = Pulse_kx;
  load ([ops_dir,'Pulse_k3.mat']);
  Pulse_k3 = Pulse_kx;
  load ([ops_dir,'Pulse_k4.mat']);
  Pulse_k4 = Pulse_kx;
  load ([ops_dir,'Pulse_k5.mat']);
  Pulse_k5 = Pulse_kx;
  load ([ops_dir,'Pulse_s.mat']);
  load ([ops_dir,'Pulse_e.mat']);

  Pulse_oc = Pulse_k1 + Pulse_k2 + Pulse_k3 + Pulse_k4 + Pulse_k5 ...
    + Pulse_s + Pulse_e;

  Land_mask = zeros(ny,nx,nz);           % This will include the marginal seas
  for k = 1:nz
    Land_mask (:,:,k) = (k > KMT(:,:));
  end % for all k levels
  
  Total_earth = Pulse_oc + Land_mask;

  Min_total = min (Total_earth(:)); % check that each oc point had a pulse
  if Min_total ~= 1
    good = false;
    disp('pulse check found missing points');
    disp('value  j   i   k    index');

    for nopulse = 1:10
     [min_val,index_I]= min(Total_earth(:));
         if min_val == 1       % if no more holes
             break
         end 
         % determine the 3-D location
         kk = floor((index_I-1) / (ny*nx)) + 1;
         ii = floor(((index_I-1) - ((kk-1)*ny*nx))/ny) + 1;
         jj = (index_I-1) - ((kk-1)*ny*nx) - ((ii-1)*ny) + 1;
         disp(['  ', ... 
             int2str(Total_earth(index_I)),'    ', ... 
             int2str(jj),'   ', ... 
             int2str(ii),'   ', ... 
             int2str(kk),'   ', ...
             int2str(index_I)]);
     Total_earth(index_I) = 1; % eliminate the worst one
    end % for up tp 10 grid cells
    disp('missing pulse locations')
    keyboard
  end % if missing point
 
  Max_total = max (Total_earth (:)); % check no oc point had multiple pulses
  if Max_total ~= 1
    good = false;
    disp('pulse check found duplicate points')
    if max( ~Land_mask .* Total_earth) <= 1
       disp ('OK - in masked-off marginal seas, harmless');
    else
    disp('value   j   i   k     index');

    for morepulse = 1:10
     [max_val,index_I]= max(Total_earth(:));
     if max_val <= 1  % check if already found all duplicate points
         break
     end
         % index_I = worst;
         % determine the 3-D location
         kk = floor((index_I-1) / (ny*nx)) + 1;
         ii = floor(((index_I-1) - ((kk-1)*ny*nx))/ny) + 1;
         jj = (index_I-1) - ((kk-1)*ny*nx) - ((ii-1)*ny) + 1;
         disp(['  ', ... 
             num2str(Total_earth(index_I)),'     ', ... 
             int2str(jj),'   ', ... 
             int2str(ii),'   ', ... 
             int2str(kk),'    ', ...
             int2str(index_I)]);
        Total_earth(index_I) = 1; % eliminate the worst one
     end % for up to 10 points
     disp('error locations')
     disp (' should be 1');
     keyboard
    end % if extra pulse not masked by land
  end % if more than 1 pulse
  return
end % function pulse_ck
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

