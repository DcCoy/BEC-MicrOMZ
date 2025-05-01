function [inputs,grid,phy,opt,params] = initialize_model(inputs,opt)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization of model based on options.m 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Timestepping via Constant dt
if opt.CONSTANT_TSTEP
    endTimey  = opt.nt*opt.dt   ./ (365*86400); % end time of simulation (years)
    histTimey = opt.hist*opt.dt ./ (365*86400); % history timestep (years)
    [inputs.dt_vec inputs.time_vec inputs.hist_time_vec inputs.hist_time_ind inputs.hist_time] = ...
        initialize_time_stepping(opt.dt,endTimey,histTimey);
    inputs.nt = length(inputs.dt_vec);
    inputs.nt_hist = length(inputs.hist_time_ind);
end

% Get bgc and ecosystem parameters
params = get_BEC_params;                     % BEC parameters 
params = BEC_overrides(params,opt);          % Overrides and additional params
params = get_physical_params(params,opt);    % Physical environmental parameters
params = physical_overrides(params,opt);     % Overrides physical parameters
if opt.EXPLICIT_MICROBES                    
    params = get_stoichiometry(params,opt.OM,opt.BM); % Elemental ratios (C:H:O:N)
    params = get_MicrOMZ_params(params,opt);          % Bacterial and small zoo parameters
    params = MicrOMZ_overrides(params,opt);           % Overrides MicrOMZ parameters based on option
	params = get_f_from_Gibbs(params,opt);            % Fraction of electrons routed to biomass synth
    params = get_growth_yields(params,opt);           % Calculates growth yields
    params = get_diffusive_uptake(params);            % Calculates diffusive uptake of O2/N2O into cell
	params = get_Vmax(params,opt);                    % Get Vmax parameters, account for stoich
end

% Get bgc tracers
inputs.nvar = length(opt.tracers);

% Vertical grid
grid.H   = abs(opt.ztop - opt.zbottom);
grid.z_w = (opt.ztop:-opt.dz:opt.zbottom)';
grid.z_r = (grid.z_w(1:end-1) + grid.z_w(2:end)) ./ 2; 
grid.Hz  = abs(diff(grid.z_w));
grid.nz  = length(grid.z_r);

% Vert Diffusion 
if opt.DEPTHVAR_KV
    phy.Kv = 0.5*(params.phy.Kv_top + params.phy.Kv_bot) + 0.5*(params.phy.Kv_top - params.phy.Kv_bot) * ...
          tanh((grid.z_w - params.phy.Kv_flex)/(0.5*params.phy.Kv_width));
else
    phy.Kv = params.phy.Kv_param.*ones(1,length(grid.z_w));
end

% Upwelling
if opt.DEPTHVAR_WUP     
    % Apply a spline to the model output w(z) profile
    load([opt.root,'/data/',opt.wup_profile]);
    tdepth = -abs(vert.depth);
    if strcmp(opt.region,'ETNP')
        tconc  = -vert.wvelZ_etnp;     
    elseif strcmp(opt.region,'ETSP')
        tconc  = -vert.wvelZ_etsp;
    end
    ibad = find(isnan(tconc));
    tdepth(ibad) = [];
    tconc(ibad) = [];     
    wupfcn = spline(tdepth,tconc);
    phy.wup = ppval(wupfcn,grid.z_w)/100.0; % convert to m/s;
    phy.wup = phy.wup';
else
    % Constant upwelling
    phy.wup = params.phy.wup .* ones(1,length(grid.z_w));
    phy.wup = phy.wup';
end
