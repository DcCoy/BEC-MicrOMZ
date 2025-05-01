function y_het = yield_sinsabaugh(y_max,B_CN,OM_CN,K_CN,EEA_CN);
% Estimate the yield of bacterial heterotrophy
%
% Inputs:
%	y_max  : The maximum yield
%	B_CN   : The C:N stoichiometric ratio of the biomass produced by heterotrophs
%	OM_CN  : The C:N stoichiometric ratio of the organic matter fuelling heterotrophy
%	K_CN   : Half-saturation coefficient for the assimialtion of C:N precursor molecules
%	EEA_CN : Eco-Enzymatic Activity (EEA) rate of carbon and nitrogen processing
%
% Outputs:
%	y_het  : The yield of bacterial heterotrophy (0-->Y_Max) in mol Biomass C per mol Organic C

S_CN  = B_CN / (OM_CN * EEA_CN);
y_het = y_max * S_CN / (S_CN + K_CN);
