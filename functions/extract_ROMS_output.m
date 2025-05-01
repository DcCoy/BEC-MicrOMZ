% Script to extract ROMS output for comparison purposes
clear all
fname = ['/data/project2/model_output/peru_chile_0p1/peru_chile_0p1_dccoy_VKV4_tune7_fixAx/avg/annual/avg_2049.nc'];

% List coordinate to extract profile
xidx = 247; % 275W
yidx = 382; % -8S

% Get list of variables to extract from ROMS
vars  = {'O2','NO2','NO3','N2O','NH4','PO4',...
		 'SiO3','Fe','DIC','Alk','DOC','DON',...
		 'DOP','DOPR','DONR',...
		 'ZOOC','SPC','SPCHL','SPFE','SPCACO3',...
		 'DIATC','DIATCHL','DIATFE','z_r'};
diag_vars = {'POC_REMIN','POC_FLUX_IN',...
		 'DOC_REMIN','DON_REMIN','O2_CONSUMPTION','O2_PRODUCTION',...
		 'AMMOX','NITROX','DENITRIF1','DENITRIF2','DENITRIF3','ANAMMOX','N2OAMMOX'};

% Get updated var names (to match MicrOMZ)
myvars = {'o2','no2','no3','n2o','nh4','po4',...
		  'sio3','fe','dic','alk','doc','don',...
		  'dop','dopr','donr',...
		  'zoo','sp','sp_chl','sp_fe','sp_caco3',...
		  'diat','diat_chl','diat_fe','z_r'};
mydiag  = {'poc_remin','poc_flux_in','doc_remin','don_remin','o2_consumption','o2_production',...
		   'ammox','nitrox','denitrif1','denitrif2','denitrif3','anammox','n2oammox'};

% Grid to depths (assumes 0:10:2000 grid)
out_grid = (-1995:10:-5);
out_grid = out_grid';
out_grid = flipud(out_grid);

% Load and restrict variables
for i = 1:length(vars)
	disp(['Loading ',vars{i}]);
    tmp.(myvars{i}) = flipud(squeeze(ncread(fname,vars{i},[xidx yidx 1 1],[1 1 inf inf])));
	tinfo = ncinfo(fname,vars{i});
	units.(myvars{i}) = tinfo.Attributes(2).Value;
end

% Interpolate
for i = 1:length(vars)-1
	disp(['Interpolating ',vars{i}]);
	for t = 1:12
		roms.(myvars{i})(:,t) = interp1(tmp.z_r(:,t),tmp.(myvars{i})(:,t),out_grid);
	end
	roms.(myvars{i})(roms.(myvars{i})<0) = 0;
end
roms.z_r = out_grid;

% Load and restrict diagnostics
for i = 1:length(diag_vars)
	disp(['Loading ',diag_vars{i}]);
	tmp.(mydiag{i}) = flipud(squeeze(ncread(fname,diag_vars{i},[xidx yidx 1 1],[1 1 inf inf])));
	tinfo = ncinfo(fname,diag_vars{i});
	diag_units.(mydiag{i}) = tinfo.Attributes(2).Value;
end

% Interpolate
for i = 1:length(diag_vars)
	disp(['Interpolating ',diag_vars{i}]);
	for t = 1:12
		roms_diag.(mydiag{i})(:,t) = interp1(tmp.z_r(:,t),tmp.(mydiag{i})(:,t),out_grid);
	end
	roms_diag.(mydiag{i})(roms_diag.(mydiag{i})<0) = 0;
end

% Save output
save('../data/ETSP_ROMS_output.mat','roms','roms_diag','units','diag_units');
