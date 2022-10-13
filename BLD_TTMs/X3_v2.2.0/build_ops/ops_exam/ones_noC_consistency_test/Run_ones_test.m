function Run_ones_test(F)
% Run_ones_test_mo sets up the external parameters that are used to run
% the losing/gaining tracer test
clear global
global gbec
gbec.P.ops_dir = F.ops_dir;
gbec.P.out_dir = F.out_dir;
gbec.P.nyears  = F.nyears;
gbec.P.timestep_hrs = F.timestep_hrs;

gbec.P.first_period = 1;
gbec.P.period_span  = 1;
gbec.P.last_period  = F.num_periods;


eval (['load ',F.ops_dir,'MET.mat MET ']);
gbec.x0 = 0*MET.iocn + 1;
sim_ones('test_ones')

return
end % function Run_ones_test_mo
