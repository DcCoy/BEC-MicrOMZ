function params = MicrOMZ_overrides(params,opt);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MicrOMZ model parameter overrides
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dps = params.bec.dps;

% ------------ %
% Heterotrophs %
% ------------ %
% General overrides
for i = 1:length(opt.hetero_ind)
	% HIGH UMAX SWITCH
	% Raise maximum growth rates to match Buchanan et al., 2024
	if opt.HIGH_UMAX
		params.(opt.hetero_ind{i}).u_max = 1.0*dps; % from 0.50
	end
	% RAPID EXCLUSION SWITCH
	% Raise linear, lower quadratic mortality rates
	if opt.RAPID_EXCLUSION
		params.(opt.hetero_ind{i}).m_l = 0.0145*dps; % from 0.01 (gets same loss assuming B = 0.05 uM C)
		params.(opt.hetero_ind{i}).m_q = 0.01*dps;   % from 0.1  (gets same loss assuming B = 0.05 uM C)
	end
end
% Specific overrides

% --------------- %
% Chemoautotrophs %
% --------------- %
% General overrides
for i = 1:length(opt.chemo_ind)
	% RAPID EXCLUSION SWITCH
	% Raise linear, lower quadratic mortality rates
	if opt.RAPID_EXCLUSION
		params.(opt.chemo_ind{i}).m_l = 0.0145*dps; % from 0.01 (gets same loss assuming B = 0.05 uM C)
		params.(opt.chemo_ind{i}).m_q = 0.01*dps;   % from 0.1  (gets same loss assuming B = 0.05 uM C)
	end
end
% Specific overrides
% HIGH UMAX SWITCH
% Raise maximum growth rates to match Buchanan et al., 2024
if opt.HIGH_UMAX
	params.aoa.u_max = 1.0*dps; % from 0.50
	params.nob.u_max = 2.0*dps; % from 1.00
	params.aox.u_max = 0.5*dps; % from 0.25
end
