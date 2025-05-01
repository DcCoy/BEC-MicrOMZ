% Script to extract surface forcing from ROMS files
clear all
fname1 = ['/data/project2/demccoy/ROMS_configs/peru_chile_0p1/frc/peru_chile_0p1_365days_1979_normal_frc_corr_dfs_river.nc'];
fname2 = ['/data/project2/demccoy/ROMS_configs/peru_chile_0p1/frc/peru_chile_0p1_era_1979-2016_clim_frc_reduced.nc'];

% List coordinate to extract profile
xidx = 247; % 275W
yidx = 382; % -8S

% Get list of variables to extract from fname1
fvars  = {'sustr','svstr','shflux','swrad','swflux'};
myvars = fvars;

% Load and restrict variables
for i = 1:length(fvars)
    frc.(myvars{i}) = flipud(squeeze(ncread(fname1,fvars{i},[xidx yidx 1],[1 1 inf])));
end

% Get list of variables to extract from fname2
fvars  = {'pco2_air','iron','dust'};
myvars = fvars;

% Load and restrict variables
for i = 1:length(fvars)
    frc.(myvars{i}) = flipud(squeeze(ncread(fname2,fvars{i},[xidx yidx 1],[1 1 inf])));
end

% Get total list of variables
fvars = {'sustr','svstr','shflux','swrad','swflux','pco2_air','iron','dust'};

% Interpolate monthly forcing to daily
days_to_avg_norm = { 1:31  ; 32:59  ; 60:90  ; 91:120 ; 121:151; 152:181;
                    182:212; 213:243; 244:273; 274:304; 305:334; 335:365};
for i = 1:length(fvars)
    if length(frc.(fvars{i}))==12
        for j = 1:365
            for k = 1:12
                idx = ismember(j,days_to_avg_norm{k});
                if idx == 1
                    break
                end
            end
            TMP.(fvars{i})(j) = frc.(fvars{i})(idx);
        end
        % Override original array with daily output
        frc.(fvars{i}) = TMP.(fvars{i})';
    else
        continue
    end
end
    
% Save
save('../data/ETSP_ROMS_forcing.mat','frc');
