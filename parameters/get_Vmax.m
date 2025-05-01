function params = get_Vmax(params,opt)
% Calculates Vmax (mol substrate / mol B(C) s)
% using max growth rate estimates
%
% umax = (Vmax)/Y
% Vmax = umax*Y
%
% Essentially, at infinite substrate concentrations ... (R / R + Kr) = 1
% Vmax will allow for max growth rate to occur
%
% For anaerobic heterotrophs, growth using alternative oxidants will be all be reduced
% For anaerobic heterotrophs, growth using reductants will be reduced (except 'nos' type)
%
% To allow anaerboic heterotrophic metabolisms to grow at max rate, use MAX_ANA_GROWTH switch

% -----------------------
% Chemoautotrophs
% -----------------------
for i = 1:length(opt.chemo_ind)
	params.(opt.chemo_ind{i}).Vmax_oxy = params.(opt.chemo_ind{i}).u_max.*params.(opt.chemo_ind{i}).y_oxy;
	params.(opt.chemo_ind{i}).Vmax_red = params.(opt.chemo_ind{i}).u_max.*params.(opt.chemo_ind{i}).y_red;
end

% -----------------------
% Heterotrophs
% -----------------------
params.aer.Vmax_oxy_aer = params.aer.u_max.*params.aer.y_oxy_aer;
params.aer.Vmax_red_aer = params.aer.u_max.*params.aer.y_red_aer;
params.aer.Vmax_oxy_ana = params.aer.Vmax_oxy_aer; % same as aerobic settings
params.aer.Vmax_red_ana = params.aer.Vmax_red_aer; % same as aerobic settings
for i = 1:length(opt.hetero_ind)
	if opt.MAX_ANA_GROWTH
		params.(opt.hetero_ind{i}).Vmax_oxy_aer = params.(opt.hetero_ind{i}).u_max.*params.(opt.hetero_ind{i}).y_oxy_aer; 
		params.(opt.hetero_ind{i}).Vmax_red_aer = params.(opt.hetero_ind{i}).u_max.*params.(opt.hetero_ind{i}).y_red_aer;
		params.(opt.hetero_ind{i}).Vmax_oxy_ana = params.(opt.hetero_ind{i}).u_max.*params.(opt.hetero_ind{i}).y_oxy_ana;
		params.(opt.hetero_ind{i}).Vmax_red_ana = params.(opt.hetero_ind{i}).u_max.*params.(opt.hetero_ind{i}).y_red_ana;
	else
		params.(opt.hetero_ind{i}).Vmax_oxy_aer = params.aer.Vmax_oxy_aer; % Same as aerobic metabolism 
		params.(opt.hetero_ind{i}).Vmax_red_aer = params.aer.Vmax_red_aer; % Same as aerobic metabolism
		params.(opt.hetero_ind{i}).Vmax_oxy_ana = params.aer.Vmax_oxy_ana; % Same as aerobic metabolism
		params.(opt.hetero_ind{i}).Vmax_red_ana = params.aer.Vmax_red_ana; % Same as aerobic metabolism
	end
end
