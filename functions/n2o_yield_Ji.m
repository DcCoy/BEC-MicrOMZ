function [yaoa_n2o_nh4,yaoa_n2o_hyb,yaoa_no2] = n2o_yield_Ji(o2,params);
% Parameterization from Ji et al.
% Need to update this to make yaoa_n2o* in fractions corresponding to mol N / mol N
%

disp('RECODE THIS MODULE');
kill kill kill

yaoa_n2o_nh4 = 1./(2 + 2/((params.aoa.Ji_a/tr.o2(k) + params.aoa.Ji_b)/100));    
yaoa_n2o_hyb = 0;    
yaoa_no2 = 1./(1 * ((params.aoa.Ji_a/tr.o2(k) + params.aoa.Ji_b)/100) + 1);


return
%yaoa_n2o = 1./(2 + 2/((Ji_a/O2 + Ji_b)/100));    
%yaoa_no2 = 1./(1 * ((Ji_a/O2 + Ji_b)/100) + 1);
