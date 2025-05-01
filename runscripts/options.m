% Set model options below

% User specific
opt.root = my_root; % Get paths           

% Model switches
opt.REGION               = 'ETSP';   % ETSP or ETNP
opt.ADVECTION            = 'UPWIND'; % 1st order Upwind advection (consider adding other methods) 
opt.DEPTHVAR_WUP         = false;    % If true, use depth-dependent upwelling velocity 
opt.DEPTHVAR_KV          = true;     % If true, use depth-dependent diffusion profile
opt.CONSTANT_TSTEP       = true;     % If true, use constant time-stepping
opt.RESTORING            = true;     % If true, turns on restoring of variables
opt.TAUZVAR              = true;     % If true, use a depth-dependent restoring time-scale
opt.HIST_VERBOSE         = true;     % prompts a message at each saving timestep
opt.REST_VERBOSE         = true;     % prompts a message when loading restart
opt.AVERAGING            = true;     % If true, save averages over 'hist' timesteps 
opt.TDEP_REMIN           = false;    % Use temperature dependent length scale for POC remin 
opt.O2_REMIN             = false;    % Use O2 dependent dissassociation length scale for POC remin
opt.N2_FULL              = false;    % If false, only track excess N2 
opt.SEDIMENTS            = false;    % If true, run sediment module via 'compute_particulate_terms'
opt.GASEXCHANGE          = true;     % If true, run airsea flux module (turns off surface bry conditions)

% Testing switches
opt.PHYSICS_OFF = false; % If true, turn off advection module
opt.BIOLOGY_OFF = false; % If true, turn off biology module

% MicrOMZ switches
opt.EXPLICIT_MICROBES = true; % Use explicit chemoautotrophs and heterotrophs
if opt.EXPLICIT_MICROBES
	opt.FACULTATIVE_MICROBES = false; % If true, heterotrophs can switch between aerobic and anaerobic metabolisms
	opt.HETERO_GRAZING       = true;  % If true, include small zooplankton type that consumes heterotrophs
	opt.HETERO_GRAZING_O2    = true;  % If true, limiting grazing at low O2 (same as zooplankton)
	opt.CHEMO_GRAZING        = true;  % If true, allow chemoautotrophs to be grazed by zooplankton
	opt.CHEMO_GRAZING_O2     = true;  % If true, limit grazing at low O2
	opt.O2_PULSE             = false; % If true, add eddy-fluxes of O2 in OMZ at distinct periods
	opt.KELLY_N2O_YIELD      = true;  % If true, use combined N2O yields from hybrid & nh4 (Kelly et al., 2024)
	opt.DNRA1                = true;  % If true, include DNRA type (starting at NO3)
	opt.DNRA2                = true;  % If true, include DNRA type (starting at NO2)
	opt.SET_AER_F            = true;  % If true, set 'f' for obligate aerobic bacteria (aer) to literature values
	opt.MAX_ANA_GROWTH       = false; % If true, allow anaerobic heterotrophs to reach their maximum growth rate
	opt.FACULTATIVE_PENALTY  = false; % If true, assign a penalty to facultative anaerobes during 'f' calculation
	opt.SINSABAUGH           = false; % If true, solve for heterotroph OM yield using approach of Sinsabaugh (2013)
	opt.INSTANT_REMIN        = false; % If true, route some mortality/grazing directly to DIC (instant remineralization)
	opt.RAPID_EXCLUSION      = false; % If true, override linear and quadratic mortality to encourage competitive exclusion
	opt.HIGH_UMAX            = false; % If true, raise default maximum growth rates (Buchanan et al., 2024)
end

% RunName and RestartFile settings
opt.RunName      = ['MicrOMZ_default_dnra2_rst']; % Set name of run
opt.RESTART      = true;                         % If true, initialize from restart
opt.RestartFile  = ['MicrOMZ_default_dnra2.mat']; % Restart file (include .mat)

% RESTORING switches
opt.rest_o2   = false; % Restore to Farfield profile
opt.rest_no3  = true;  % Restore to Farfield profile (ROMS override)
opt.rest_po4  = true;  % Restore to Farfield profile (ROMS override)
opt.rest_no2  = false; % Restore to Farfield profile
opt.rest_n2o  = false; % Restore to Farfield profile
opt.rest_fe   = true;  % Restore to initial conditions (ROMS override)
opt.rest_sio3 = true;  % Restore to initial conditions (ROMS override)

% FIXED CONCENTRATION switches
% The below, if true, will set constant nutrient concentrations
opt.fix_fe   = false; % Fix concentrations to initial conditions
opt.fix_sio3 = false; % Fix concentrations to initial conditions
opt.fix_po4  = false; % Fix concentrations to initial conditions
opt.fix_no3  = false; % Fix concentrations to initial conditions (surface only)
opt.fix_nh4  = false; % Fix concentrations to initial conditions (surface only)

% Data sources for wup, Tau, Restart
opt.wup_profile  = 'vertical_CESM.mat';                     % Vertical velocities
opt.Tau_profiles = 'Tau_restoring.mat';                     % Depth dependent Restoring timescale
opt.Farfield     = ['farfield_',opt.REGION,'_gridded.mat']; % Restoring profiles
opt.InitialFile  = [opt.REGION,'_ROMS_restart.mat'];        % Initial conditions file
opt.ForcingFile  = [opt.REGION,'_ROMS_forcing.mat'];        % Annual (repeat) forcing file

% Vertical grid
opt.ztop    = 0;     % Top depth (m)
opt.zbottom = -2000; % Bottom depth (m)
opt.dz      = 10;    % Width of each cell (m)

% Timestepping:
if (1) % yearly (5-day output)
    opt.nyears = 20;                   % Number of years
    opt.dt     = 3600*6;               % Timestep (seconds)
    opt.hist   = ((86400.*5)/opt.dt);  % Save a snapshot every 'hist' timesteps
	opt.avg    = opt.hist;             % Average over 'hist' timesteps
elseif (1) % daily (daily output)
    opt.nyears = (50/365);             % Number of days (in years) 
    opt.dt     = 3600;                 % Timestep (seconds)
    opt.hist   = (86400/opt.dt);       % Save a snapshot every 'hist' timesteps
	opt.avg    = opt.hist;             % Average over 'hist' timesteps
else % debug (10 day, hourly output)
    opt.nyears = (1/365);             % Number of days (in years) 
    opt.dt     = 3600*4;               % Timestep (seconds)
    opt.hist   = (86400/opt.dt);       % Save a snapshot every 'hist' timesteps
	opt.avg    = opt.hist;             % Average over 'hist' timesteps
end
opt.nt = ((86400*365)/opt.dt) .* opt.nyears;    % Simulation length (by # of timesteps);

% Prognostic tracers 
opt.tracers = {'po4','no3','sio3','nh4','fe',...
               'no2','n2','n2o','o2',...
               'dic','alk','doc','don',...
               'dop','dopr','donr',...
               'zoo','sp','sp_chl','sp_fe','sp_caco3',...
               'diat','diat_chl','diat_fe','diat_si',...
               'diaz','diaz_chl','diaz_fe'};
opt.titles = {'PO$^{3-}_4$','NO$^{-}_3$','SiO$_3$','NH$^{+}_4$','Fe',...
              'NO$^{-}_2$','N$_2$','N$_2$O','O$_2$',...
              'DIC','Alk','DOC','DON',...
              'DOP','DOPr','DONr',...
              'Zoo','Sp','Sp$_{chl}$','Sp$_{fe}$','Sp$_{caco3}$',...
              'Diat','Diat$_{chl}$','Diat$_{fe}$','Diat$_{sio3}$',...
              'Diaz','Diaz$_{chl}$','Diaz$_{fe}$'};

% If other modules are requested, add the extra tracers
if opt.EXPLICIT_MICROBES
	% Add default tracers (MicrOMZ + BEC tracers)
    tracers_to_add = {'docr','szoo','aoa','nob','aox',...
                      'aer','nar','nai','nao','nir','nio','nos'};
    titles_to_add = {'DOCr','sZoo','AOA','NOB','AOX',...
                     'AER','NAR','NAI','NAO','NIR','NIO','NOS'};
	% Include DNRA1 type
	if opt.DNRA1
		tracers_to_add{end+1} = 'dnra1';
		titles_to_add{end+1}  = 'DNRA1';
	end
	% Include DNRA2 type
	if opt.DNRA2
		tracers_to_add{end+1} = 'dnra2';
		titles_to_add{end+1}  = 'DNRA2';
	end
	% Uptake tracers and titles
    for i = 1:length(tracers_to_add)
        opt.tracers{end+1} = tracers_to_add{i};
        opt.titles{end+1} = titles_to_add{i};
    end

	% Add auto, hetero and chemo indexes
	opt.auto_ind   = {'sp','diat','diaz'};
	opt.chemo_ind  = {'aoa','nob','aox'};
	opt.hetero_ind = {'aer','nar','nai','nao','nir','nio','nos'};
	if opt.DNRA1
		opt.hetero_ind{end+1} = 'dnra1';
	end
	if opt.DNRA2
		opt.hetero_ind{end+1} = 'dnra2';
	end
	
	% Get indices for each functional type
	% Autotrophs
	opt.sp_ind   = find(strcmp('sp',opt.auto_ind)==1);
	opt.diat_ind = find(strcmp('diat',opt.auto_ind)==1);
	opt.diaz_ind = find(strcmp('diaz',opt.auto_ind)==1);
	% Chemoautotrophs
	opt.aoa_ind = find(strcmp('aoa',opt.chemo_ind)==1);
	opt.nob_ind = find(strcmp('nob',opt.chemo_ind)==1);
	opt.aox_ind = find(strcmp('aox',opt.chemo_ind)==1);
	% Heterotrophs
	opt.aer_ind = find(strcmp('aer',opt.hetero_ind)==1);
	opt.nar_ind = find(strcmp('nar',opt.hetero_ind)==1);
	opt.nai_ind = find(strcmp('nai',opt.hetero_ind)==1);
	opt.nao_ind = find(strcmp('nao',opt.hetero_ind)==1);
	opt.nir_ind = find(strcmp('nir',opt.hetero_ind)==1);
	opt.nio_ind = find(strcmp('nio',opt.hetero_ind)==1);
	opt.nos_ind = find(strcmp('nos',opt.hetero_ind)==1);
	if opt.DNRA1
		opt.dnra1_ind = find(strcmp('dnra1',opt.hetero_ind)==1); 
	end
	if opt.DNRA2
		opt.dnra2_ind = find(strcmp('dnra2',opt.hetero_ind)==1); 
	end
end

% Stoichiometry of organic matter and bacteria
%          C    H    O    N    P
opt.OM = [117  297   85   16   1];      % Anderson & Sarmiento, 1995 ... as in BEC
if opt.EXPLICIT_MICROBES
%              C     H     O     N    P
    opt.BM = [55.0  77.0  22.0  11.0  1];  % Zimmerman et al., 2014
end
