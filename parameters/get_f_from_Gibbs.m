function params = get_f_from_Gibbs(params,opt)
% Calculates the fraction of electrons used for biomass synthesis
% following Rittman & McCarty (Environmental Biotechnology) 

% Extract stoichiometry
OM = params.stoich.OM;
BM = params.stoich.BM;

% Get ratios
% Organice matter
dOM = params.stoich.dOM;
cOM = params.stoich.cOM;
hOM = params.stoich.hOM;
oOM = params.stoich.oOM;
nOM = params.stoich.nOM;
pOM = params.stoich.pOM;

% Bacteria
dBM = params.stoich.dBM;
cBM = params.stoich.cBM;
hBM = params.stoich.hBM;
oBM = params.stoich.oBM;
nBM = params.stoich.nBM;
pBM = params.stoich.pBM;

% ----------------------------------------------
% Get Gibbs free energies of various reactions
% ----------------------------------------------
% Define environmental parameters for calculations
T  = 25 + 273.15; % K
pH = 7.5;         % pH scale
R  = 8.3145;      % J/(mol K), ideal gas constant
H  = 10.^(-pH);   % Conversion of pH to mol/L

% Get Gibbs free energy of formation for relevant compounds and phases
% All values from Amend & Shock (2001) assuming 25C
deltaGf.NO3 = -110.91 * 1e3; % J/mol (aq) [NO3-]
deltaGf.NO2 = -32.22  * 1e3; % J/mol (aq) [NO2-]
deltaGf.H2O = -237.18 * 1e3; % J/mol (l)
deltaGf.N2O =  113.38 * 1e3; % J/mol (aq)
deltaGf.H   =    0.00 * 1e3; % J/mol (aq) [H+]
deltaGf.NO  =  102.06 * 1e3; % J/mol (aq)
deltaGf.N2  =   18.18 * 1e3; % J/mol (aq)
deltaGf.O2  =   16.54 * 1e3; % J/mol (aq) 
deltaGf.NH4 =  -79.45 * 1e3; % J/mol (aq)        

% Get Gibbs standard free energy for half reactions
% Oxidation with O2: (O2 --> H2O)
deltaGf.oxy_aer = ...
    (1/2)*(deltaGf.H2O) - ... % product
    (1/4)*(deltaGf.O2)  - ... % reactant
    (1/1)*(deltaGf.H);        % reactant
% Denitrification 1: (NO3 --> NO2) 
deltaGf.oxy_den1 = ...
    (1/2)*(deltaGf.NO2) + ... % product 
    (1/2)*(deltaGf.H2O) - ... % product 
    (1/2)*(deltaGf.NO3) - ... % reactant
    (1/1)*(deltaGf.H);        % reactant
% Denitrification 2: (NO2 --> N2O)
deltaGf.oxy_den2 = ...
    (1/4)*(deltaGf.N2O) + ... % product
    (3/4)*(deltaGf.H2O) - ... % product
    (1/2)*(deltaGf.NO2) - ... % reactant
    (3/2)*(deltaGf.H);        % reactant
% Denitrification 3: (N2O --> N2)
deltaGf.oxy_den3 = ...
    (1/2)*(deltaGf.N2)  + ... % product
    (1/2)*(deltaGf.H2O) - ... % product
    (1/2)*(deltaGf.N2O) - ... % reactant
    (1/1)*(deltaGf.H);        % reactant
% Denitrification 4: (NO3 --> N2O)
deltaGf.oxy_den4 = ...
    (1/8)*(deltaGf.N2O) + ... % product
    (5/8)*(deltaGf.H2O) - ... % product
    (1/4)*(deltaGf.NO3) - ... % reactant
    (5/4)*(deltaGf.H);        % reactant
% Denitrification 5: (NO2 --> N2)
deltaGf.oxy_den5 = ...
    (1/6)*(deltaGf.N2)  + ... % product
    (2/3)*(deltaGf.H2O) - ... % product
    (1/3)*(deltaGf.NO2) - ... % reactant
    (4/3)*(deltaGf.H);        % reactant
% Denitrification 6: (NO3 --> N2)
deltaGf.oxy_den6 = ...
    (1/10)*(deltaGf.N2)  + ... % product
     (3/5)*(deltaGf.H2O) - ... % product 
     (1/5)*(deltaGf.NO3) - ... % reactant
     (6/5)*(deltaGf.H);        % reactant
if opt.DNRA1
	% DNRA 1: (NO3 --> NO2 --> NH4)
	deltaGf.oxy_dnra1 = ...
		(1/8)*(deltaGf.NH4) + ... % product
		(3/8)*(deltaGf.H2O) - ... % product
		(1/8)*(deltaGf.NO3) - ... % reactant
		(5/4)*(deltaGf.H);        % reactant
end
if opt.DNRA2
	% DNRA 2: (NO2 --> NH4)
	deltaGf.oxy_dnra2 = ...
		(1/6)*(deltaGf.NH4) + ... % product
		(1/3)*(deltaGf.H2O) - ... % product
		(1/6)*(deltaGf.NO2) - ... % reactant
		(4/3)*(deltaGf.H);        % reactant
end
	
% Set average concentrations for calculations
ini.NO3 = 30e-6; % mol/l
ini.NO2 =  1e-6; % mol/l
ini.N2O =  1e-8; % mol/l
ini.N2  =  1e-4; % mol/l
ini.O2  =  1e-6; % mol/l
ini.NH4 =  1e-6; % mol/l

% Get nonstandard Gibbs free energy
deltaGf.oxy_aer   = deltaGf.oxy_aer   + R*T*log(1/ini.O2^(1/4)/H);
deltaGf.oxy_den1  = deltaGf.oxy_den1  + R*T*log(ini.NO2^(1/2)/ini.NO3^(1/2)/H);
deltaGf.oxy_den2  = deltaGf.oxy_den2  + R*T*log(ini.N2O^(1/4)/ini.NO2^(1/2)/H^(3/2));
deltaGf.oxy_den3  = deltaGf.oxy_den3  + R*T*log(ini.N2^(1/2)/ini.N2O^(1/2)/H);
deltaGf.oxy_den4  = deltaGf.oxy_den4  + R*T*log(ini.N2O^(1/8)/ini.NO3^(1/4)/H^(5/4));
deltaGf.oxy_den5  = deltaGf.oxy_den5  + R*T*log(ini.N2^(1/6)/ini.NO2^(1/3)/H^(4/3));
deltaGf.oxy_den6  = deltaGf.oxy_den6  + R*T*log(ini.N2^(1/10)/ini.NO3^(1/5)/H^(6/5));
if opt.DNRA1
	deltaGf.oxy_dnra1 = deltaGf.oxy_dnra1 + R*T*log(ini.NH4^(1/8)/ini.NO3^(1/8)/H^(5/4));  
end
if opt.DNRA2
	deltaGf.oxy_dnra2 = deltaGf.oxy_dnra2 + R*T*log(ini.NH4^(1/6)/ini.NO2^(1/6)/H^(4/3));  
end

% Get Gibbs reaction energy for oxidation of organic matter
om_energy = 3.33e3; % J/g cells (Rittman & McCarty)
deltaGf.red_OM = om_energy*(cOM*12 + hOM*1 + oOM*16 + nOM*14)/dOM;

% Get Gibbs reaction energy to form bacterial biomass from pyruvate (Rittman & McCarty)
b_energy = 3.33e3; % J/g cells (Rittman & McCarty)
deltaGf.pyr_BM = b_energy*(cBM*12 + hBM*1 + oBM*16 + nBM*14)/dBM;

% Set efficiency of electron transfer
ep = 0.6;

% Get Gibbs reaction energy for conversion of OM to pyruvate
% Energy required to conver the carbon source to activated pyruvate (168 M&R)
% 35.09 via Table 5.4, O-21
deltaGf.OM_pyr = 35.09e3 - deltaGf.red_OM;

% Get total energy needed for cell synthesis
n = 1; % assuming deltaGf.OM_pry > 1
deltaGf.syn_BM = deltaGf.OM_pyr/ep^n + deltaGf.pyr_BM/ep; 

% Override aerobic functional type 'f' to backout energy needed for cell synthesis
if opt.SET_AER_F
    % Set 'f' of aerobic heterotrophs to literature values
    tmp.f = 0.2; % Robinson et al., 2008 (f = 0.1 -- 0.2)
    % Back-out deltaGf.syn_BM based on 'f'
    deltaGf.syn_BM = ...
        -(ep*(deltaGf.oxy_aer-deltaGf.red_OM) - ...
          ep*tmp.f*(deltaGf.oxy_aer-deltaGf.red_OM));
    deltaGf.syn_BM = deltaGf.syn_BM/tmp.f;
end

% Override energy needed for cell synthesis for anaerobic bacteria
if opt.FACULTATIVE_MICROBES & opt.FACULTATIVE_PENALTY
	step1facpen = 0.9; % aerobic penalty 
	step2facpen = 0.8; % aerobic penalty
	step3facpen = 0.7; % aerobic penalty
	anafacpen   = 0.6; % anaerobic penalty
else
    step1facpen = 1; % aerobic penalty
    step2facpen = 1; % aerobic penalty
    step3facpen = 1; % aerobic penalty
    anafacpen   = 1; % anaerobic penalty
end

% Get penalty for multi-step denitrification
P = 0.16;

% ---------------------------------------
% Heterotrophic bacteria
% ---------------------------------------

% AER (OM + O2 --> B)
f = calc_f(deltaGf.oxy_aer,deltaGf.red_OM,deltaGf.syn_BM,ep);
params.aer.f = f; % NOTE: This is set to 0.2 above via deltaGf.syn_BM calculation

% NAR (OM + NO3 --> B + NO2)
f = calc_f(deltaGf.oxy_den1,deltaGf.red_OM,deltaGf.syn_BM,ep);
f = f*(1-P*(1-1)); % 1 step, no penalty
params.nar.f     = f .* anafacpen;
params.nar.f_fac = params.aer.f .* step1facpen;

% NIR (OM + NO2 --> B + N2O)
f = calc_f(deltaGf.oxy_den2,deltaGf.red_OM,deltaGf.syn_BM,ep);
f = f*(1-P*(1-1)); % 1 step, no penalty
params.nir.f = f .* anafacpen;
params.nir.f_fac = params.aer.f .* step1facpen;

% NOS (OM + N2O --> B + N2)
f = calc_f(deltaGf.oxy_den3,deltaGf.red_OM,deltaGf.syn_BM,ep);
f = f*(1-P*(1-1)); % 1 step, no penalty
params.nos.f = f .* anafacpen;
params.nos.f_fac = params.aer.f .* step1facpen;

% NAI (OM + NO3 --> B + N2O)
f = calc_f(deltaGf.oxy_den4,deltaGf.red_OM,deltaGf.syn_BM,ep);
f = f*(1-P*(2-1)); % 2 steps
params.nai.f = f .* anafacpen;
params.nai.f_fac = params.aer.f .* step2facpen;

% NIO (OM + NO2 --> B + N2)
f = calc_f(deltaGf.oxy_den5,deltaGf.red_OM,deltaGf.syn_BM,ep);
f = f*(1-P*(2-1)); % 2 steps
params.nio.f = f .* anafacpen;
params.nio.f_fac = params.aer.f .* step2facpen;

% NAO (OM + NO3 --> B + N2)
f = calc_f(deltaGf.oxy_den6,deltaGf.red_OM,deltaGf.syn_BM,ep);
f = f*(1-P*(3-1)); % 3 steps
params.nao.f = f .* anafacpen;
params.nao.f_fac = params.aer.f .* step3facpen;

% DNRA1 (OM + NO3 --> B + NH4)
if opt.DNRA1
	f = calc_f(deltaGf.oxy_dnra1,deltaGf.red_OM,deltaGf.syn_BM,ep);
	f = f*(1-P*(2-1)); % 2 steps
	params.dnra1.f = f .* anafacpen;
	params.dnra1.f_fac = params.aer.f .* step2facpen;
end

% DNRA2 (OM + NO2 --> B + NH4)
if opt.DNRA2
	f = calc_f(deltaGf.oxy_dnra2,deltaGf.red_OM,deltaGf.syn_BM,ep);
	f = f*(1-P*(1-1)); % 1 steps
	params.dnra2.f = f .* anafacpen;
	params.dnra2.f_fac = params.aer.f .* step2facpen;
end

% ---------------------------------------
% Chemoautotrophic bacteria
% ---------------------------------------

% AOA
yNH4_aoa = 0.0245;               % (mol Bn  / mol NH4) Bayer et al., 2022; Zakem et al., 2022
yNH4_aoa = yNH4_aoa*(cBM/nBM);   % (mol Bc  / mol NH4)
yNH4_aoa = 1/yNH4_aoa;           % (mol NH4 / mol Bc)
f = dBM/(6*(yNH4_aoa - nBM));    % Solve for f, yNH4_aoa = (dBM/f)*(1/6 + f*nBM/dBM)
params.aoa.f = f;

% NOB
yNO2_nob = 0.0126;             % (mol Bn  / mol NO2) Bayer et al., 2022
yNO2_nob = yNO2_nob*(cBM/nBM); % (mol Bc  / mol NO2)
yNO2_nob = 1/yNO2_nob;         % (mol NO2 / mol Bc)
f = (dBM)/(2*yNO2_nob);        % Solve for f, yNO2_nob = (dBM/f)*(1/2)
params.nob.f = f;

% AOX
f = 0.07; % Lotti et al., 2014
x = 0.64; % Solve for X (via anammox stoich of products NO3/NO2) Daims 2014, Chesson 2000
params.aox.f = f;
params.aox.x = x;
