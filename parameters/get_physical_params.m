function params = get_physical_params(params,opt);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default physical parameters for the model
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Time (days to seconds)
days2secs  = 1/(86400);
years2secs = 1/(86400 * 365);

% Particle sinking (m/s)
params.phy.wsink = -20 .* days2secs; % constant speed (opt.varsink==0)

% Upwelling speed (m/s)
% Depth-dependent velocity requires a forcing file 
params.phy.wup = 10 .* years2secs; % convert m/yr to m/s (here, 10 m/yr)

% Diffusion (m^2/s^2)
params.phy.Kv       = 2.0 * 1.701e-5; % constant vertical diffusion coefficient in m^2/s
params.phy.Kv_top   = 0.70 * 2.0 * 1.701e-5;
params.phy.Kv_bot   = 1.00 * 2.0 * 1.701e-5;
params.phy.Kv_flex  = -250;
params.phy.Kv_width = 300;

% Physical scalings for restoring
params.phy.Rh = 1.0;          % unitless scaling for sensitivity analysis. Default is 1.0
params.phy.Lh = 4000.0 * 1e3; % m - horizontal scale
% if you chose constant restoring timescales
if ~opt.TAUZVAR
    params.phy.Kh = 1000; % m2/s - horizontal diffusion
    params.phy.Uh = 0.05; % m/s - horizontal advection
end

% Temperature
params.phy.TCoefArr       = 0.8;    % Arrhenius temperature coefficient
params.phy.TAeArr         = -4000;  % Arrhenius temperature coefficient
params.phy.TrefArr        = 293.15; % Arrhenius temperature
params.phy.Tkel           = 273.15; % C to Kelvin
params.phy.tcline_shallow = 150;    % top of thermocline 
params.phy.tcline_deep    = 500;    % bottom of thermocline
params.phy.temp_top       = 12;     % Bry cond. for surface temperature (C)
params.phy.temp_bot       = 2;      % Bry cond. for deep temperature (C)

% Light
params.phy.par_in = 700;        % Bry cond. for surface irradiance      (W / m2)
