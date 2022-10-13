function [MET] = get_MET(ops_dir)
% Returns the MET geometry structure 
  eval(['load ',ops_dir,'MET.mat MET']);
  return
end % get_METbygrain



