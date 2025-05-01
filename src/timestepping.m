function [out,bgc] = timestepping(inputs,grid,phy,bgc,opt,params,frc)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Timestepping module
%
% Steps:
% (0) - Set boundary conditions
% (1) - Advect and diffuse tracers 
% (2) - Calculate BGC sources-minus-sinks
% OPTIONAL:
% (3) - Calculate air-sea flux of gas tracers (o2, n2, n2o)
% (4) - Lateral restoring of tracers
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath([opt.root,'/src/']);

% Initialize solutions
out.vars = zeros(inputs.nt_hist,inputs.nvar,grid.nz);
out.adv  = out.vars;
out.dfz  = out.vars;
out.sms  = out.vars;
out.flux = zeros(inputs.nt_hist,3,1); % only applied at surface for o2, n2o, n2
if opt.RESTORING
    out.rest = out.vars;
end

% Initialize averaged solutions
% Diag averaging set at timestep == 1 below
if opt.AVERAGING
	out.vars_avg = out.vars;
	out.adv_avg  = out.adv;
	out.dfz_avg  = out.dfz;
	out.sms_avg  = out.sms;
	out.flux_avg = out.flux;
	if opt.RESTORING
		out.rest_avg = out.rest;
	end
end

% Initialize tracer matrices 
for i = 1:length(opt.tracers)
    % Dimensions = 2 (tracer in, tracer out)
    tmp.(opt.tracers{i}) = zeros(grid.nz,2);
    % Set in == restart
    if opt.RESTART
        tmp.(opt.tracers{i})(:,1) = bgc.rst(:,i);
    else
        tmp.(opt.tracers{i})(:,1) = bgc.(opt.tracers{i});  
    end
end

% Set up pulsing test
if opt.O2_PULSE
	inputs.o2_pulse_period = 30;
	inputs.o2_pulse_conc   = 0.5;
    dt       = inputs.dt_vec(1);                % timestep (seconds)
    period   = inputs.o2_pulse_period .* 86400; % days 2 seconds
    sarray   = (1:inputs.nt).*dt;               % array of all times sampled in the model
    pulse_dt = [period:period:max(sarray)];     % array of timesteps corresponding to period
    pulse_dt = pulse_dt ./ dt;                  % index of 'i' corresponding to period
end

% Initialize tracer average arrays
if opt.AVERAGING
	for i = 1:length(opt.tracers)
		tmp_avg.(opt.tracers{i}) = tmp.(opt.tracers{i})(:,1); 
	end
end

% Get forcing fields    
forcing_fields = fields(frc);

% Start time-stepping
for i = 1:inputs.nt

    % Gets current timestep in seconds
    dt = inputs.dt_vec(i);

    % Get current day (for forcing, which is at daily frequency)
    this_day = floor((inputs.time_vec(i)./86400)+1);
    for j = 1:10000 
        if this_day > 365
            this_day = this_day - 365; 
        end
        if this_day <= 365
            break
        end
    end

    % Set input forcing to daily value
    for j = 1:length(forcing_fields)
        this_frc.(forcing_fields{j}) = frc.(forcing_fields{j})(this_day);
    end

    % Evaluate tendencies at current timestep 
    for j = 1:length(opt.tracers)
        tr.(opt.tracers{j}) = tmp.(opt.tracers{j})(:,1);
    end

    % Check for pulsing of O2
    if opt.O2_PULSE
        if ismember(i,pulse_dt)
            % Set O2 to pulse magnitude
			ind = find(tr.o2 <= 0);
			for j = 1:length(ind)
				tr.o2(ind(j)) = rand(1).*inputs.o2_pulse_conc;
			end
        end
    end

    % ------ %
    % Step 1 %
    % ------ %
    % Collect gas exchange at current timestep 
    if opt.GASEXCHANGE
        % Calculate airsea flux (mmol/m2s), setting other tracers to 0
        [flux] = airsea_flux(tr,phy,this_frc,params,grid,opt);
        flux_fields = fields(flux);
        for j = 1:length(opt.tracers)
            if ~ismember(opt.tracers{j},flux_fields)
                flux.(opt.tracers{j}) = 0;
            end
        end
    else
        % Set all fluxes to 0
        for j = 1:length(opt.tracers)
            flux.(opt.tracers{j}) = 0;
        end
    end

    % ------ %
    % STEP 2 %
    % ------ %
    % Evaluate sources-and-sinks at current time-step
    % NOTE: Updates 'tr' to remove negative numbers
	if ~opt.BIOLOGY_OFF
		if ismember(i,inputs.hist_time_ind) | opt.AVERAGING;
			% Save diagnostics at timestep    
			[sms,diag,tr] = sources_sinks(tr,phy,grid,opt,params,this_frc,1);
		else
			% Skip building of diagnostic structure
			[sms,~,tr] = sources_sinks(tr,phy,grid,opt,params,this_frc,0);
		end
	else
		for j = 1:length(opt.tracers)
			sms.(opt.tracers{j}) = zeros(grid.nz,1);
			diag = struct;
		end
	end

    % ------ %
    % STEP 3 %
    % ------ %
    % Advect and diffuse tracers at current timestep 
    if strcmp(opt.ADVECTION,'UPWIND');
        for j = 1:length(opt.tracers)
            % Get advection
            adv.(opt.tracers{j}) = upwind(tr.(opt.tracers{j}),phy.wup,grid.Hz);
            % Get diffusion
            dfz.(opt.tracers{j}) = diffusion(tr.(opt.tracers{j}),phy.Kv,grid.Hz);
        end
    elseif strcmp(opt.ADVECTION,'CRANK-NICHOLSON');
        % Call '[c1 c2 c3] = crank-nicholson(phy,opt)'    
        % Need to reorganize the code though, since it already applies the 'dt'
        % Basically, need to collect tendency*dt at each of these 'STEP'(s)
        % Also, look into how to solve this at boundaries
    else    
        disp('Check options.m for opt.ADVECTION');
        kill
    end

    % ------ %
    % STEP 4 %
    % ------ %
    % Restore tracers at current timestep 
    if opt.RESTORING
        % Add lateral restoring, set in options.m 
        rest = restoring(opt,inputs,grid,tr);
        rest_fields = fields(rest);
        for j = 1:length(opt.tracers)
            if ~ismember(opt.tracers{j},rest_fields)
                rest.(opt.tracers{j}) = zeros(size(grid.z_r));
            end
        end
    else
        % Set restoring to 0
        for j = 1:length(opt.tracers)
            rest.(opt.tracers{j}) = zeros(size(grid.z_r));
        end
    end

    % ------ %
    % STEP 5 %
    % ------ %
    % Update tracers at forward timestep, including tendencies from:
    % (1) gas exchange, distributed over the cell (k == 1 only)
    % (2) sources-minus-sinks
    % (3) advection and diffusion
    % (4) lateral restoring
    for j = 1:length(opt.tracers)
        % SMS, ADV, REST
        tmp.(opt.tracers{j})(:,2) = tr.(opt.tracers{j}) + ((...
            sms.(opt.tracers{j})  + ...
            adv.(opt.tracers{j})  + ...
            dfz.(opt.tracers{j})  + ...
            rest.(opt.tracers{j})).*dt);
        % AIR-SEA (surface only)
        tmp.(opt.tracers{j})(1,2) = tmp.(opt.tracers{j})(1,2) + ...
             (flux.(opt.tracers{j})/grid.Hz(1)).*dt;
    end

    % ------ %
    % STEP 6 %
    % ------ %
    % Add bottom boundary conditions to forward timestep
    % (constant source to be upwelled) 
    for j = 1:length(opt.tracers)
        tmp.(opt.tracers{j})(end,2) = bgc.(opt.tracers{j})(end);
        % If no gas exchange, also set surface boundary condition
        if ~opt.GASEXCHANGE
            tmp.(opt.tracers{j})(1,2) = bgc.(opt.tracers{j})(1);
        end
    end

    % ------ %
    % STEP 7 %
    % ------ %
	% (OPTIONAL) Fix tracer concentrations?
	if opt.fix_fe
		tmp.fe(:,2) = tmp.fe(:,1);
	end
	if opt.fix_po4
		tind = find(strcmp('po4',opt.tracers)==1);
		tmp.po4(:,2) = tmp.po4(:,1);
	end
	if opt.fix_sio3
		tmp.sio3(:,2) = tmp.sio3(:,1);
	end
	if opt.fix_no3
		zidx = find(grid.z_r>-50); % only fix surface concentrations
		tmp.no3(zidx,2) = tmp.no3(zidx,1);
	end
	if opt.fix_nh4
		zidx = find(grid.z_r>-50); % only fix surface concentrations
		tmp.nh4(zidx,2) = tmp.nh4(zidx,1);
	end
		

    % Old tracer equals new tracer
    for j = 1:length(opt.tracers)
        tmp.(opt.tracers{j})(:,1) = tmp.(opt.tracers{j})(:,2);
		% Sum tracers into average structure
		if opt.AVERAGING
			tmp_avg.(opt.tracers{j}) = tmp_avg.(opt.tracers{j}) + tmp.(opt.tracers{j})(:,2); 				
		end
    end

	% Also save averages from other terms
	if opt.AVERAGING
		if i == 1
			% Get fields of diagnostic structure
			ff = fields(diag);
			% Set first entries
			diag_avg = diag;
			adv_avg  = adv;
			dfz_avg  = dfz;
			sms_avg  = sms;
			flux_avg = flux;
			if opt.RESTORING
				rest_avg = rest;
			end
			% Get timesteps counter
			NT = 1;
		else
			% Add diagnostics
			for j = 1:length(ff)
				diag_avg.(ff{j}) = diag_avg.(ff{j}) + diag.(ff{j});
			end
			% Add advection, diffusion, sms, fluxes, restoring
			for j = 1:length(opt.tracers)
				adv_avg.(opt.tracers{j})  = adv_avg.(opt.tracers{j}) + adv.(opt.tracers{j});
				dfz_avg.(opt.tracers{j})  = dfz_avg.(opt.tracers{j}) + dfz.(opt.tracers{j});
				sms_avg.(opt.tracers{j})  = sms_avg.(opt.tracers{j}) + sms.(opt.tracers{j});
				flux_avg.(opt.tracers{j}) = flux_avg.(opt.tracers{j}) + flux.(opt.tracers{j}); 
				if opt.RESTORING
					rest_avg.(opt.tracers{j}) = rest_avg.(opt.tracers{j}) + rest.(opt.tracers{j});
				end
			end
			NT = NT + 1;
		end
	end

    % Save history files and diagnostics? 
    if ismember(i,inputs.hist_time_ind);
        % Record progress
        idx = find(inputs.hist_time_ind == i);
        if opt.HIST_VERBOSE
            disp(['Saving step #' num2str(idx) '/' num2str(inputs.nt_hist)]);
        end
        % Go through tracers and save fields
        % Start flux counter
        fluxcnt = 1;
        for j = 1:length(opt.tracers)

            % Concentration
            out.vars(idx,j,:) = tmp.(opt.tracers{j})(:,1);
			if opt.AVERAGING	
				out.vars_avg(idx,j,:) = tmp_avg.(opt.tracers{j})./NT;
			end

            % Advection (uM/s) 
            if ~opt.GASEXCHANGE
                out.adv(idx,j,2:end-1) = adv.(opt.tracers{j})(2:end-1);
                out.dfz(idx,j,2:end-1) = dfz.(opt.tracers{j})(2:end-1); 
            else
                out.adv(idx,j,1:end-1) = adv.(opt.tracers{j})(1:end-1);
                out.dfz(idx,j,1:end-1) = dfz.(opt.tracers{j})(1:end-1); 
            end
			if opt.AVERAGING
				out.adv_avg(idx,j,:) = adv_avg.(opt.tracers{j})./(NT-1); 
			end

            % Sources minus sinks (uM/s)
            if ~opt.GASEXCHANGE
                out.sms(idx,j,2:end-1) = sms.(opt.tracers{j})(2:end-1);
            else
                out.sms(idx,j,1:end-1) = sms.(opt.tracers{j})(1:end-1);
            end
			if opt.AVERAGING
				out.sms_avg(idx,j,:) = sms_avg.(opt.tracers{j})./(NT-1);
			end

            % Air-sea fluxes (uM/s)
			if opt.GASEXCHANGE
				if ismember(opt.tracers{j},flux_fields);
					out.flux(idx,fluxcnt,1) = (flux.(opt.tracers{j})/grid.Hz(1));
					if opt.AVERAGING
						out.flux_avg(idx,fluxcnt) = flux_avg.(opt.tracers{j})./(grid.Hz(1)*(NT-1));
					end
					fluxcnt = fluxcnt + 1;
				end
			end

            % Save restoring (uM/s)
            if opt.RESTORING
                if ~opt.GASEXCHANGE
                    out.rest(idx,j,2:end-1) = rest.(opt.tracers{j})(2:end-1);
                else
                    out.rest(idx,j,1:end-1) = rest.(opt.tracers{j})(1:end-1);
                end
				out.rest_avg(idx,j,:) = rest_avg.(opt.tracers{j})./(NT-1);
            end
        end

        % Also save diagnostics into structure
        out.diag(idx) = diag;
		if opt.AVERAGING
			for j = 1:length(ff)
				out.diag_avg(idx).(ff{j}) = diag_avg.(ff{j})./(NT-1);
				% Reset diagnostics
				diag_avg.(ff{j}) = zeros(size(diag_avg.(ff{j})));
			end
		end

		% Reset remaining averaging structures 
		if opt.AVERAGING
			for j = 1:length(opt.tracers)
				tmp_avg.(opt.tracers{j})  = tmp.(opt.tracers{j})(:,1);
				adv_avg.(opt.tracers{j})  = zeros(size(adv_avg.(opt.tracers{j})));
				dfz_avg.(opt.tracers{j})  = zeros(size(dfz_avg.(opt.tracers{j})));
				sms_avg.(opt.tracers{j})  = zeros(size(sms_avg.(opt.tracers{j})));
				if opt.GASEXCHANGE
					flux_avg.(opt.tracers{j}) = zeros(size(flux_avg.(opt.tracers{j}))); 
				end
				if opt.RESTORING
					rest_avg.(opt.tracers{j}) = zeros(size(rest_avg.(opt.tracers{j})));
				end
			end
			% Reset averaging counter
			NT = 1;
		end
    end 
end  

