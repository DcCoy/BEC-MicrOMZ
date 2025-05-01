function [bgc,phy] = initialize_tracers(inputs,grid,phy,opt,params)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set initial conditions or use restart file based on options
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A very small number
eps = 1e-10;

% Load final timestep from previous run
if opt.RESTART 
    % Load previous BGC output from restart file 
    disp(['Loading restart : ' opt.root '/restart/' opt.RestartFile]);
    tmp = load([opt.root,'/output/',opt.RestartFile]);
    for i = 1:length(opt.tracers)
        bgc.(opt.tracers{i}) = tmp.bgc.(opt.tracers{i})(end,:);
        bgc.rst(:,i)         = tmp.bgc.(opt.tracers{i})(end,:);
    end
    phy.temp = tmp.phy.temp;
    phy.salt = tmp.phy.salt;
% Load initial conditions from ROMS 
elseif strcmp(opt.REGION,'ETSP');
    % Get ROMS initial conditions file
    fname = ['../data/',opt.InitialFile];
    load(fname);
    fname = ['../data/ETSP_ROMS_output.mat'];
    load(fname);
    % Get depth    
    z_r = tmp.z_r;
    % Interpolate temp and salt to grid
    phy.temp = double(interp1(z_r,tmp.temp,grid.z_r));
    phy.salt = double(interp1(z_r,tmp.salt,grid.z_r));
    % Go through each variable and attempt to load IC, else set to eps
    for i = 1:length(opt.tracers)
        % Try to load ROMS initial conditions
        try
            this_var = tmp.(opt.tracers{i});
            bgc.(opt.tracers{i}) = double(interp1(z_r,this_var,grid.z_r));
            bgc.rst(:,i) = double(interp1(z_r,this_var,grid.z_r));
        % If not available (i.e., for microbes), set to a very small number
        catch
            bgc.(opt.tracers{i}) = ones(size(grid.z_r)).*eps;
            bgc.rst(:,i) = ones(size(grid.z_r)).*eps;
            % Further reduce initial anaerboic microbes
            % This way, the aerobic heterotroph functional type (aer) will grow more rapidly
            % and competitively exclude the others in instances where 'opt.FACULTATIVE' is true 
            if ismember(opt.tracers{i},{'nar','nai','nao','nir','nio','nos'});
                bgc.(opt.tracers{i}) = bgc.(opt.tracers{i})./100;
                bgc.rst(:,i) = bgc.rst(:,i)./100;
            end
			% Set szoo to 1/10th of zoo 
			if strcmp(opt.tracers{i},'szoo');
				bgc.(opt.tracers{i}) = bgc.zoo./10;
			end
        end
        % Override N2O and N2 to atmospheric saturation values
        if strcmp(opt.tracers{i},'n2o');
            bgc.n2o = n2osatu(phy.temp,phy.salt,params) .* params.bec.xn2o;
            bgc.rst(:,i) = bgc.n2o;
        elseif strcmp(opt.tracers{i},'n2');
            if opt.N2_FULL
                bgc.n2  = n2satu(phy.temp,phy.salt,params);
                bgc.rst(:,i) = bgc.n2;
            else
                bgc.n2 = zeros(size(grid.z_r));
                bgc.rst(:,i) = zeros(size(grid.z_r));
            end
        end
        % Override surface O2 to atmospheric saturation values
        if strcmp(opt.tracers{i},'o2');
            o2sfc = o2satu(phy.temp(1),phy.salt(1));
            bgc.o2(1:3) = o2sfc;
            bgc.rst(1:3,i) = o2sfc;
        end
        % Override initial conditions for trace nutrients
        if ismember(opt.tracers{i},{'fe','sio3','po4'})
            bgc.(opt.tracers{i}) = nanmean(roms.(opt.tracers{i}),2); % annual average from monthly output
        end
        % Override initial conditions for N2O 
        if ismember(opt.tracers{i},{'n2o'})
            bgc.(opt.tracers{i}) = nanmean(roms.(opt.tracers{i}),2); % annual average from monthly output
        end
    end
	% Override initial DOM ratios
	bgc.don  = bgc.doc.*(params.stoich.nOM./params.stoich.cOM);
	bgc.dop  = bgc.doc.*(params.stoich.pOM./params.stoich.cOM);
	bgc.donr = bgc.docr.*(params.stoich.nOM./params.stoich.cOM);
	bgc.dopr = bgc.docr.*(params.stoich.pOM./params.stoich.cOM);
end

% OVERRIDE FOR NO BIOLOGY
if opt.BIOLOGY_OFF
	for i = 1:length(opt.tracers)
		bgc.(opt.tracers{i}) = ones(size(grid.z_r)).*eps;
		bgc.(opt.tracers{i})(end) = 1000;
	end
end
