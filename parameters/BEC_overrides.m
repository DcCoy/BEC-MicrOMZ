function params = BEC_overrides(params,opt);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEC model parameter overrides
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Needed to route DOM to refractory pools
if opt.EXPLICIT_MICROBES
    params.bec.docrefract = 0.3; % fraction of DOC to DOCr
    params.bec.donrefract = 0.3; % fraction of DON to DONr
    params.bec.doprefract = 0.3; % fraction of DOP to DOPr
end

% If false, remove routing of mortality and grazing to DIC 
if ~opt.INSTANT_REMIN
	% Portion of grazed OrgC to DIC now sent to DOC
	%for i = 1:length(opt.auto_ind)	
	%	params.(opt.auto_ind{i}).graze_doc = ...
	%		(1 - params.(opt.auto_ind{i}).graze_zoo - params.(opt.auto_ind{i}).graze_poc);
	%end
	% Turn off instant remineralization (DIC production) of losses not routed to POC
	%params.bec.labile_ratio = 0;
end
