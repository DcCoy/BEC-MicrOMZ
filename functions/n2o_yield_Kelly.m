function [yaoa_n2o_nh4,yaoa_n2o_hyb,yaoa_no2] = n2o_yield_Kelly(o2,params);
% Function to calculate the partioning of NH4-oxidation into NO2 and N2O
% Via Figure 9, https://doi.org/10.5194/bg-21-3215-2024 (Kelly et al., 2024)
%
% We have NH4-oxidation in the model, which excretes N (mol N/s)
% NH4 + O2 --> Nexc
%
% We calculate Nexc, but need to partition that into NO2 and N2O (both in mol N/s)
% via 3 pathways:
% p1 = NO2 excretion               (NH4 + O2 --> NO2)
% p2 = N2O excretion from NH4 only (NH4 + O2 --> 0.5 N2O)
% p3 = N2O excretion from hybrid   (NH4 + NO2 + O2 --> N2O)
%
% So we need to get the fraction of these pathways in terms of NH4 consumed (mol NH4/s) 
%

% Get fraction of (mol N/s) routed to p2
f2 = (2.2*exp(-1.5*o2))./100;   

% Get fraction of (mol N/s) routed to p3
f3 = (20.4*exp(-0.58*o2))./100; 

% Allow N2O production pathways to have small background yield
if (f2 + f3) < (params.aoa.Ji_b./100)
	target_fraction = (params.aoa.Ji_b./100);
	f2_new = target_fraction.*(f2/(f2+f3));
	f3_new = target_fraction.*(f3/(f2+f3));
	f2 = f2_new;
	f3 = f3_new;
end

% Get fraction of (mol N/s) routed to p1
f1 = 1 - (f2 + f3);

% Collect yields (%) in terms of mol N/s 
yaoa_no2     = f1; % mol N/s goes to NO2
yaoa_n2o_nh4 = f2; % mol N/s goes to N2O
yaoa_n2o_hyb = f3; % mol N/s goes to N2O
