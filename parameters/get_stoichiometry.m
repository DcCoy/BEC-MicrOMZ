function params = get_stoichiometry(params,OM,BM);
% ----------------------------------------------
% Set stoichiometry of OM and bacteria
% ----------------------------------------------
%      C    H    O    N 
params.stoich.OM = [OM]; % set in options.m
params.stoich.BM = [BM]; % set in options.m

% Normalize all to C (since the model currency is in mol C/m3)
params.stoich.OM = params.stoich.OM ./ (params.stoich.OM(1));
params.stoich.BM = params.stoich.BM ./ (params.stoich.BM(1));

% Get composition of OM based on inputs
params.stoich.cOM = params.stoich.OM(1);
params.stoich.hOM = params.stoich.OM(2);
params.stoich.oOM = params.stoich.OM(3);
params.stoich.nOM = params.stoich.OM(4);
params.stoich.pOM = params.stoich.OM(5);

% Get composition of bacteria based on inputs
params.stoich.cBM = params.stoich.BM(1);
params.stoich.hBM = params.stoich.BM(2);
params.stoich.oBM = params.stoich.BM(3);
params.stoich.nBM = params.stoich.BM(4);
params.stoich.pBM = params.stoich.BM(5);

% Get number of electrons required to oxidize organic matter and bacterial biomass
% with one mole of phosphorous
params.stoich.dOM = ...
    4*params.stoich.cOM + ...
    1*params.stoich.hOM - ...
    2*params.stoich.oOM - ...
    3*params.stoich.nOM;
params.stoich.dBM = ...
    4*params.stoich.cBM + ...
    1*params.stoich.hBM - ...
    2*params.stoich.oBM - ...
    3*params.stoich.nBM; 
