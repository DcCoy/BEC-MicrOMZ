%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Template BEC-MicrOMZ runscript 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;

% Add root path and paths to model directories
paths_init

% Load run options
options
disp(['Running ',opt.RunName]);

% Initialize the model
disp('Initializing model based on options.m');
inputs = struct;
[inputs,grid,phy,opt,params] = initialize_model(inputs,opt); 

% Initialize tracers 
disp('Initializing tracers based on options.m');
[bgc,phy] = initialize_tracers(inputs,grid,phy,opt,params);

% Initialiaze forcing
[frc] = initialize_forcing(opt);

% Initialize Restoring
if opt.RESTORING
    if ~opt.TAUZVAR
        phy.tauh = ((params.phy.Uh/params.phy.Lh + 2*params.phy.Kh/params.phy.Lh^2)^-1)*params.phy.Rh;
    else
        phy.tauh = [];
    end
    inputs = initialize_restoring(inputs,grid,phy,opt,params,bgc);
end

% Run the model 
disp('Start timestepping');
disp(' ');
tic;
[out,bgc] = timestepping(inputs,grid,phy,bgc,opt,params,frc);
bgc.RunTime = toc;
disp('Finished timestepping');
disp(['Runtime : ' num2str(bgc.RunTime)]);
disp(' ');

% Update bgc and diag structures with output
for i = 1:length(opt.tracers);
    % Concentrations
    bgc.(opt.tracers{i}) = single(squeeze(out.vars(:,i,:)));
	if opt.AVERAGING
		bgc.([opt.tracers{i},'_avg']) = single(squeeze(out.vars_avg(:,i,:)));
	end
    % Advective transport
    diag.([opt.tracers{i},'_adv']) = single(squeeze(out.adv(:,i,:)));
	if opt.AVERAGING
		diag.([opt.tracers{i},'_adv_avg']) = single(squeeze(out.adv_avg(:,i,:)));
	end
    % Diffusive transport
    diag.([opt.tracers{i},'_dfz']) = single(squeeze(out.dfz(:,i,:)));
	if opt.AVERAGING
		diag.([opt.tracers{i},'_dfz_avg']) = single(squeeze(out.dfz_avg(:,i,:)));
	end
    % Sources-minus-sinks
    diag.([opt.tracers{i},'_sms']) = single(squeeze(out.sms(:,i,:)));
	if opt.AVERAGING
		diag.([opt.tracers{i},'_sms_avg']) = single(squeeze(out.sms_avg(:,i,:)));
	end
    % Air-sea fluxes
    if opt.GASEXCHANGE
        if strcmp(opt.tracers{i},'n2o');
            diag.n2o_flux = single(out.flux(:,1));
			if opt.AVERAGING
				diag.n2o_flux_avg = single(out.flux_avg(:,1));
			end
        elseif strcmp(opt.tracers{i},'n2');
            diag.n2_flux = single(out.flux(:,2));
			if opt.AVERAGING
				diag.n2_flux_avg = single(out.flux_avg(:,2));
			end
        elseif strcmp(opt.tracers{i},'o2');
            diag.o2_flux = single(out.flux(:,3));
			if opt.AVERAGING
				diag.o2_flux_avg = single(out.flux_avg(:,3));
			end
        end
    end
    % Restoring
    if opt.RESTORING
        diag.([opt.tracers{i},'_rest']) = single(squeeze(out.rest(:,i,:)));
		if opt.AVERAGING
			diag.([opt.tracers{i},'_rest_avg']) = single(squeeze(out.rest_avg(:,i,:)));
		end
    end
end

% Update diagnostic with output
diag_fields = fields(out.diag(1));
if ~isempty(diag_fields)
    for i = 1:length(diag_fields)
        for j = 1:inputs.nt_hist;
            diag.(diag_fields{i})(j,:) = single(out.diag(j).(diag_fields{i}));
			if opt.AVERAGING
				diag.([diag_fields{i},'_avg'])(j,:) = single(out.diag_avg(j).(diag_fields{i}));
			end
        end
    end
end
clear out

% Save output
disp('Saving output');
disp(' ');
t = datetime('today');
fname = [opt.RunName,'.mat'];
disp(['Saving ',fname]);
save(['../output/',fname],'bgc','phy','inputs','params','opt','grid','diag','frc');

% Announce finish
disp('done!');
