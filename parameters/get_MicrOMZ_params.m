function params = get_MicrOMZ_params(params,opt)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization of microbial ecosystem parameters
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get shortcut to DPS
days2secs = params.bec.dps;

% Small zoo
params.szoo.u_max   = 4.5*days2secs;  % Max growth rate (s-1)
params.szoo.m_l     = 0*days2secs;    % Linear mortality (uM C-1 s-1) 
params.szoo.m_q     = 0.5*days2secs;  % Quadratic mortality (uM C-1 s-1)
params.szoo.Y       = 0.5;            % Fraction of digestion
params.szoo.bmin    = 1e-6;           % Minimum biomass where losses stop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Heterotrophic %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Aerobic heterotroph
% NOTE: ana = aer for simplification
params.aer.u_max     = 0.5*days2secs;  % Max growth rate (1/s) following Buchanan (Peligibacter estimate) 
params.aer.CN        = 5;              % CN ratio
params.aer.CP        = 55;             % CP ratio
params.aer.CFe       = 1/20e-6;        % CFe ratio
params.aer.K_oxy_aer = 0.20;           % Half-saturation uptake of O2
params.aer.K_oxy_ana = 0.20;           % Half-saturation uptake of O2
params.aer.K_doc     = 10.0;           % Half-saturation uptake of DOC
params.aer.K_docr    = 1000;           % Half-saturation uptake of DOCr
params.aer.m_l       = 0.01*days2secs; % Linear mortality
params.aer.m_q       = 0.1*days2secs;  % Quadratic mortality
params.aer.bmin      = 1e-10;          % Minimum biomass (mmol C / m3)
params.aer.kszoo     = 2.4;            % Half saturation for grazing
params.aer.mortdoc   = 0.5;            % Fraction of mortality to DOC (%)
params.aer.mortpoc   = 0.3;            % Fraction of mortality to POC (%)
params.aer.mortdic   = 0.2;            % Fraction of mortality to DIC (%)
params.aer.oxyup_aer = 'michaelis';    % Switch for diffusive oxidant uptake (aerobic)
params.aer.oxyup_ana = 'michaelis';    % Switch for diffusive oxidant uptake (anaerboic)
params.aer.oxy_aer   = 'o2';           % Oxidant (aerobic)
params.aer.oxy_ana   = 'o2';           % Oxidant (anaerobic)

% NO3 to NO2 denitrifier
params.nar.u_max     = 0.5*days2secs;  % Max growth rate (1/s) 
params.nar.CN        = 5;              % CN ratio
params.nar.CP        = 55;             % CP ratio
params.nar.CFe       = 1/20e-6;        % CFe ratio
params.nar.K_oxy_aer = 0.20;           % Half-saturation uptake of O2
params.nar.K_oxy_ana = 10.0;           % Half-saturation uptake of NO3 (Penn 2016, Parsonage 1985, Betlach 1981)
params.nar.K_doc     = 10.0;           % Half-saturation uptake of DOC
params.nar.K_docr    = 1000;           % Half-saturation uptake of DOCr
params.nar.m_l       = 0.01*days2secs; % Linear mortality
params.nar.m_q       = 0.1*days2secs;  % Quadratic mortality
params.nar.bmin      = 1e-10;          % Minimum biomass (mmol C / m3)
params.nar.kszoo     = 2.4;            % Half saturation for grazing
params.nar.mortdoc   = 0.5;            % Fraction of mortality to DOC (%)
params.nar.mortpoc   = 0.3;            % Fraction of mortality to POC (%)
params.nar.mortdic   = 0.2;            % Fraction of mortality to DIC (%)
params.nar.oxyup_aer = 'michaelis';    % Switch for diffusive oxidant uptake (aerobic)
params.nar.oxyup_ana = 'michaelis';    % Switch for diffusive oxidant uptake (anaerboic)
params.nar.oxy_aer   = 'o2';           % Oxidant (aerobic)
params.nar.oxy_ana   = 'no3';          % Oxidant (anaerobic)

% NO3 to N2O denitrifier
params.nai.u_max     = 0.5*days2secs;  % Max growth rate (1/s) 
params.nai.CN        = 5;              % CN ratio
params.nai.CP        = 55;             % CP ratio
params.nai.CFe       = 1/20e-6;        % CFe ratio
params.nai.K_oxy_aer = 0.20;           % Half-saturation uptake of O2
params.nai.K_oxy_ana = 10.0;           % Half-saturation uptake of NO3 (Penn 2016, Parsonage 1985, Betlach 1981)
params.nai.K_doc     = 10.0;           % Half-saturation uptake of DOC
params.nai.K_docr    = 1000;           % Half-saturation uptake of DOCr
params.nai.m_l       = 0.01*days2secs; % Linear mortality
params.nai.m_q       = 0.1*days2secs;  % Quadratic mortality
params.nai.bmin      = 1e-10;          % Minimum biomass (mmol C / m3)
params.nai.kszoo     = 2.4;            % Half saturation for grazing
params.nai.mortdoc   = 0.5;            % Fraction of mortality to DOC (%)
params.nai.mortpoc   = 0.3;            % Fraction of mortality to POC (%)
params.nai.mortdic   = 0.2;            % Fraction of mortality to DIC (%)
params.nai.oxyup_aer = 'michaelis';    % Switch for diffusive oxidant uptake (aerobic)
params.nai.oxyup_ana = 'michaelis';    % Switch for diffusive oxidant uptake (anaerboic)
params.nai.oxy_aer   = 'o2';           % Oxidant (aerobic)
params.nai.oxy_ana   = 'no3';          % Oxidant (anaerobic)

% NO3 to N2 denitrifier
params.nao.u_max     = 0.5*days2secs;  % Max growth rate (1/s) 
params.nao.CN        = 5;              % CN ratio
params.nao.CP        = 55;             % CP ratio
params.nao.CFe       = 1/20e-6;        % CFe ratio
params.nao.K_oxy_aer = 0.20;           % Half-saturation uptake of O2
params.nao.K_oxy_ana = 10.0;           % Half-saturation uptake of NO3 (Penn 2016, Parsonage 1985, Betlach 1981)
params.nao.K_doc     = 10.0;           % Half-saturation uptake of DOC
params.nao.K_docr    = 1000;           % Half-saturation uptake of DOCr
params.nao.m_l       = 0.01*days2secs; % Linear mortality
params.nao.m_q       = 0.1*days2secs;  % Quadratic mortality
params.nao.bmin      = 1e-10;          % Minimum biomass (mmol C / m3)
params.nao.kszoo     = 2.4;            % Half saturation for grazing
params.nao.mortdoc   = 0.5;            % Fraction of mortality to DOC (%)
params.nao.mortpoc   = 0.3;            % Fraction of mortality to POC (%)
params.nao.mortdic   = 0.2;            % Fraction of mortality to DIC (%)
params.nao.oxyup_aer = 'michaelis';    % Switch for diffusive oxidant uptake (aerobic)
params.nao.oxyup_ana = 'michaelis';    % Switch for diffusive oxidant uptake (anaerboic)
params.nao.oxy_aer   = 'o2';           % Oxidant (aerobic)
params.nao.oxy_ana   = 'no3';          % Oxidant (anaerobic)

% NO2 to N2O denitrifier
params.nir.u_max     = 0.5*days2secs;  % Max growth rate (1/s) 
params.nir.CN        = 5;              % CN ratio
params.nir.CP        = 55;             % CP ratio
params.nir.CFe       = 1/20e-6;        % CFe ratio
params.nir.K_oxy_aer = 0.20;           % Half-saturation uptake of O2
params.nir.K_oxy_ana = 10.0;           % Half-saturation uptake of NO2 (Penn 2016, Parsonage 1985, Betlach 1981)
params.nir.K_doc     = 10.0;           % Half-saturation uptake of DOC
params.nir.K_docr    = 1000;           % Half-saturation uptake of DOCr
params.nir.m_l       = 0.01*days2secs; % Linear mortality
params.nir.m_q       = 0.1*days2secs;  % Quadratic mortality
params.nir.bmin      = 1e-10;          % Minimum biomass (mmol C / m3)
params.nir.kszoo     = 2.4;            % Half saturation for grazing
params.nir.mortdoc   = 0.5;            % Fraction of mortality to DOC (%)
params.nir.mortpoc   = 0.3;            % Fraction of mortality to POC (%)
params.nir.mortdic   = 0.2;            % Fraction of mortality to DIC (%)
params.nir.oxyup_aer = 'michaelis';    % Switch for diffusive oxidant uptake (aerobic)
params.nir.oxyup_ana = 'michaelis';    % Switch for diffusive oxidant uptake (anaerboic)
params.nir.oxy_aer   = 'o2';           % Oxidant (aerobic)
params.nir.oxy_ana   = 'no2';          % Oxidant (anaerobic)

% NO2 to N2 denitrifier
params.nio.u_max     = 0.5*days2secs;  % Max growth rate (1/s) 
params.nio.CN        = 5;              % CN ratio
params.nio.CP        = 55;             % CP ratio
params.nio.CFe       = 1/20e-6;        % CFe ratio
params.nio.K_oxy_aer = 0.20;           % Half-saturation uptake of O2
params.nio.K_oxy_ana = 10.0;           % Half-saturation uptake of NO2 (Penn 2016, Parsonage 1985, Betlach 1981)
params.nio.K_doc     = 10.0;           % Half-saturation uptake of DOC
params.nio.K_docr    = 1000;           % Half-saturation uptake of DOCr
params.nio.m_l       = 0.01*days2secs; % Linear mortality
params.nio.m_q       = 0.1*days2secs;  % Quadratic mortality
params.nio.bmin      = 1e-10;          % Minimum biomass (mmol C / m3)
params.nio.kszoo     = 2.4;            % Half saturation for grazing
params.nio.mortdoc   = 0.5;            % Fraction of mortality to DOC (%)
params.nio.mortpoc   = 0.3;            % Fraction of mortality to POC (%)
params.nio.mortdic   = 0.2;            % Fraction of mortality to DIC (%)
params.nio.oxyup_aer = 'michaelis';    % Switch for diffusive oxidant uptake (aerobic)
params.nio.oxyup_ana = 'michaelis';    % Switch for diffusive oxidant uptake (anaerboic)
params.nio.oxy_aer   = 'o2';           % Oxidant (aerobic)
params.nio.oxy_ana   = 'no2';          % Oxidant (anaerobic)

% N2O to N2 denitrifier
params.nos.u_max     = 0.5*days2secs;  % Max growth rate (1/s) 
params.nos.CN        = 5;              % CN ratio
params.nos.CP        = 55;             % CP ratio
params.nos.CFe       = 1/20e-6;        % CFe ratio
params.nos.K_oxy_aer = 0.20;           % Half-saturation uptake of O2
params.nos.K_oxy_ana = 0.4;            % Half-saturation uptake of N2O (Penn 2016, Parsonage 1985, Betlach 1981)
params.nos.K_doc     = 10.0;           % Half-saturation uptake of DOC
params.nos.K_docr    = 1000;           % Half-saturation uptake of DOCr
params.nos.m_l       = 0.01*days2secs; % Linear mortality
params.nos.m_q       = 0.1*days2secs;  % Quadratic mortality
params.nos.bmin      = 1e-10;          % Minimum biomass (mmol C / m3)
params.nos.kszoo     = 2.4;            % Half saturation for grazing
params.nos.mortdoc   = 0.5;            % Fraction of mortality to DOC (%)
params.nos.mortpoc   = 0.3;            % Fraction of mortality to POC (%)
params.nos.mortdic   = 0.2;            % Fraction of mortality to DIC (%)
params.nos.oxyup_aer = 'michaelis';    % Switch for diffusive oxidant uptake (aerobic)
params.nos.oxyup_ana = 'michaelis';    % Switch for diffusive oxidant uptake (anaerboic)
params.nos.oxy_aer   = 'o2';           % Oxidant (aerobic)
params.nos.oxy_ana   = 'n2o';          % Oxidant (anaerobic)

% NO3 to NH4 DNRA1
if opt.DNRA1
	params.dnra1.u_max     = 0.5*days2secs;  % Max growth rate (1/s) 
	params.dnra1.CN        = 5;              % CN ratio
	params.dnra1.CP        = 55;             % CP ratio
	params.dnra1.CFe       = 1/20e-6;        % CFe ratio
	params.dnra1.K_oxy_aer = 0.20;           % Half-saturation uptake of O2
	params.dnra1.K_oxy_ana = 10.0;           % Half-saturation uptake of NO3
	params.dnra1.K_doc     = 10.0;           % Half-saturation uptake of DOC
	params.dnra1.K_docr    = 1000;           % Half-saturation uptake of DOCr
	params.dnra1.m_l       = 0.01*days2secs; % Linear mortality
	params.dnra1.m_q       = 0.1*days2secs;  % Quadratic mortality
	params.dnra1.bmin      = 1e-10;          % Minimum biomass (mmol C / m3)
	params.dnra1.kszoo     = 2.4;            % Half saturation for grazing
	params.dnra1.mortdoc   = 0.5;            % Fraction of mortality to DOC (%)
	params.dnra1.mortpoc   = 0.3;            % Fraction of mortality to POC (%)
	params.dnra1.mortdic   = 0.2;            % Fraction of mortality to DIC (%)
	params.dnra1.oxyup_aer = 'michaelis';    % Switch for diffusive oxidant uptake (aerobic)
	params.dnra1.oxyup_ana = 'michaelis';    % Switch for diffusive oxidant uptake (anaerboic)
	params.dnra1.oxy_aer   = 'o2';           % Oxidant (aerobic)
	params.dnra1.oxy_ana   = 'no3';          % Oxidant (anaerobic)
end

% NO2 to NH4 DNRA2
if opt.DNRA2
	params.dnra2.u_max     = 0.5*days2secs;  % Max growth rate (1/s) 
	params.dnra2.CN        = 5;              % CN ratio
	params.dnra2.CP        = 55;             % CP ratio
	params.dnra2.CFe       = 1/20e-6;        % CFe ratio
	params.dnra2.K_oxy_aer = 0.20;           % Half-saturation uptake of O2
	params.dnra2.K_oxy_ana = 10.0;           % Half-saturation uptake of NO3
	params.dnra2.K_doc     = 10.0;           % Half-saturation uptake of DOC
	params.dnra2.K_docr    = 1000;           % Half-saturation uptake of DOCr
	params.dnra2.m_l       = 0.01*days2secs; % Linear mortality
	params.dnra2.m_q       = 0.1*days2secs;  % Quadratic mortality
	params.dnra2.bmin      = 1e-10;          % Minimum biomass (mmol C / m3)
	params.dnra2.kszoo     = 2.4;            % Half saturation for grazing
	params.dnra2.mortdoc   = 0.5;            % Fraction of mortality to DOC (%)
	params.dnra2.mortpoc   = 0.3;            % Fraction of mortality to POC (%)
	params.dnra2.mortdic   = 0.2;            % Fraction of mortality to DIC (%)
	params.dnra2.oxyup_aer = 'michaelis';    % Switch for diffusive oxidant uptake (aerobic)
	params.dnra2.oxyup_ana = 'michaelis';    % Switch for diffusive oxidant uptake (anaerboic)
	params.dnra2.oxy_aer   = 'o2';           % Oxidant (aerobic)
	params.dnra2.oxy_ana   = 'no2';          % Oxidant (anaerobic)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Chemoautotrophic %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% AOA
params.aoa.u_max     = 1.0*days2secs;  % Max growth rate (1/s) 
params.aoa.CN        = 5;              % CN ratio
params.aoa.CP        = 55;             % CP ratio
params.aoa.CFe       = 1/20e-6;        % CFe ratio
params.aoa.K_oxy     = 0.333;          % Half-saturation uptake of O2 (Bristow et al.)
params.aoa.K_red     = 0.134;          % Half-saturation uptake of NH4 (Bristow et al.)
params.aoa.m_l       = 0.01*days2secs; % Linear mortality
params.aoa.m_q       = 0.1*days2secs;  % Quadratic mortality
params.aoa.Ji_a      = 0.2;            % non-dimensional (Ji et al., 2018) 
params.aoa.Ji_b      = 0.08;           % 1/(mmolO2/m3)   (Ji et al., 2018)
params.aoa.bmin      = 1e-10;          % Minimum biomass (mmol C / m3)
params.aoa.kzoo      = 2.4;            % Half saturation for grazing (mmol C / m3)
params.aoa.mortdoc   = 0.35;           % Fraction of mortality to DOC (%)
params.aoa.mortdic   = 0.35;           % Fraction of mortality to DOC (%)
params.aoa.mortpoc   = 0.3;            % Fraction of mortality to POC (%)
params.aoa.oxyup     = 'michaelis';    % Switch for diffusive oxidant uptake
params.aoa.oxy       = 'o2';           % Oxidant
params.aoa.red       = 'nh4';          % Reductant

% NOB
params.nob.u_max     = 1.5*days2secs;  % Max growth rate (1/s) 
params.nob.CN        = 5;              % CN ratio
params.nob.CP        = 55;             % CP ratio
params.nob.CFe       = 1/20e-6;        % CFe ratio
params.nob.K_oxy     = 0.778;          % Half-saturation uptake of O2 (Bristow et al.)
params.nob.K_red     = 0.254;          % Half-saturation uptake of NO2 (Sun et al.)
params.nob.m_l       = 0.01*days2secs; % Linear mortality
params.nob.m_q       = 0.1*days2secs;  % Quadratic mortality
params.nob.bmin      = 1e-10;          % Minimum biomass (mmol C / m3)
params.nob.kzoo      = 2.4;            % Half saturation for grazing (mmol C / m3)
params.nob.mortdoc   = 0.35;           % Fraction of mortality to DOC (%)
params.nob.mortdic   = 0.35;           % Fraction of mortality to DOC (%)
params.nob.mortpoc   = 0.3;            % Fraction of mortality to POC (%)
params.nob.oxyup     = 'michaelis';    % Switch for diffusive oxidant uptake
params.nob.oxy       = 'o2';           % Oxidant
params.nob.red       = 'no2';          % Reductant

% AOX
params.aox.u_max     = 0.25*days2secs; % Max growth rate (1/s) 
params.aox.CN        = 5;              % CN ratio
params.aox.CP        = 55;             % CP ratio
params.aox.CFe       = 1/20e-6;        % CFe ratio
params.aox.K_oxy     = 0.45;           % Half-saturation uptake of NH4
params.aox.K_red     = 0.45;           % Half-saturation uptake of NO2
params.aox.m_l       = 0.01*days2secs; % Linear mortality
params.aox.m_q       = 0.1*days2secs;  % Quadratic mortality
params.aox.bmin      = 1e-10;          % Minimum biomass (mmol C / m3)
params.aox.kzoo      = 2.4;            % Half saturation for grazing (mmol C / m3)
params.aox.mortdoc   = 0.35;           % Fraction of mortality to DOC (%)
params.aox.mortdic   = 0.35;           % Fraction of mortality to DOC (%)
params.aox.mortpoc   = 0.3;            % Fraction of mortality to POC (%)
params.aox.oxyup     = 'michaelis';    % Switch for diffusive oxidant uptake
params.aox.oxy       = 'no2';          % Oxidant
params.aox.red       = 'nh4';          % Reductant
