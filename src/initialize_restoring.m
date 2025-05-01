function  inputs = initialize_restoring(inputs,grid,phy,opt,params,bgc);
% -------------------------------------------------------------------- %
% Interpolate profile by creating a cubic spline 
% NOTE: modify substituting [0;tconc;0] if flat-slope endings are needed
% WARNING: avoid extrapolation by using profiles from 0-4000 m (even if sparse)
% -------------------------------------------------------------------- %

% Load restoring data (profiles, Tau)
load([opt.root,'/data/',opt.Farfield]);
load([opt.root,'/data/',opt.Tau_profiles]);

% Get farfield file
inputs.restore = [];
eval(['Farfield = farfield_',opt.REGION,'_gridded;']);

% Grab tauZvar from data files
if isempty(phy.tauh)
    if strcmp(opt.REGION,'ETSP');
        rest.tauh = ((rest.currentZ_etsp./(params.phy.Lh) + ...
            rest.kappaZ_etsp.*2./((params.phy.Lh)^2)).^-1);
    elseif strcmp(opt.REGION,'ETNP');
        rest.tauh = ((rest.currentZ_etnp./(params.phy.Lh) + ...
            rest.kappaZ_etnp.*2./((params.phy.Lh)^2)).^-1);
    end
    tdepth = -abs(rest.depth);
    tconc  = rest.tauh;     
    ibad = find(isnan(tconc));
    tdepth(ibad) = [];
    tconc(ibad) = [];     
    inputs.restore.tauh_pre = spline(tdepth,tconc);
end

% ---------------- %
% FARFIELD restoring
% ---------------- %
% PO4
if opt.rest_po4
    tdepth = -abs(Farfield.zgrid);
    tconc  = Farfield.po4;
    ibad = find(isnan(tconc));
    tdepth(ibad) = [];
    tconc(ibad) = [];
    inputs.restore.po4_cs = spline(tdepth,tconc);
    inputs.restore.po4_cout = ppval(inputs.restore.po4_cs,grid.z_r);    
else
    inputs.restore.po4_cout = zeros(size(grid.z_r));    
end

% NO3
if opt.rest_no3 
    tdepth = -abs(Farfield.zgrid);
    tconc  = Farfield.no3;
    ibad = find(isnan(tconc));
    tdepth(ibad) = [];
    tconc(ibad) = [];
    inputs.restore.no3_cs = spline(tdepth,tconc);
    inputs.restore.no3_cout = ppval(inputs.restore.no3_cs,grid.z_r);    
else
    inputs.restore.no3_cout = zeros(size(grid.z_r));    
end

% O2
if opt.rest_o2 
    tdepth = -abs(Farfield.zgrid);
    tconc  = Farfield.o2;
    ibad = find(isnan(tconc));
    tdepth(ibad) = [];
    tconc(ibad) = [];
    inputs.restore.o2_cs = spline(tdepth,tconc);
    inputs.restore.o2_cout = ppval(inputs.restore.o2_cs,grid.z_r);    
else
    inputs.restore.o2_cout = zeros(size(grid.z_r));    
end

% N2O
if opt.rest_n2o
    tdepth = -abs(Farfield.zgrid);
    tconc  = Farfield.n2o;
    ibad = find(isnan(tconc));
    tdepth(ibad) = [];
    tconc(ibad) = [];
    inputs.restore.n2o_cs = spline(tdepth,tconc);
    inputs.restore.n2o_cout = ppval(inputs.restore.n2o_cs,grid.z_r);    
else
    inputs.restore.n2o_cout = zeros(size(grid.z_r));    
end

% NO2
if opt.rest_no2
    tdepth = -abs(Farfield.zgrid);
    tconc  = Farfield.no2;
    ibad = find(isnan(tconc));
    tdepth(ibad) = [];
    tconc(ibad) = [];
    inputs.restore.no2_cs = spline(tdepth,tconc);
    inputs.restore.no2_cout = ppval(inputs.restore.no2_cs,grid.z_r);    
else
    inputs.restore.no2_cout = zeros(size(grid.z_r));    
end

% --------------- %
% Initial restoring
% --------------- %
% Fe
if opt.rest_fe
    tdepth = -abs(grid.z_r);
    tconc  = bgc.fe;
    inputs.restore.fe_cs = spline(tdepth,tconc);
    inputs.restore.fe_cout = ppval(inputs.restore.fe_cs,grid.z_r);    
end    

% SiO3
if opt.rest_sio3
    tdepth = -abs(grid.z_r);
    tconc  = bgc.sio3;
    inputs.restore.sio3_cs = spline(tdepth,tconc);
    inputs.restore.sio3_cout = ppval(inputs.restore.sio3_cs,grid.z_r);    
end    

% Override restoring to ROMS solutions
fname = ['../data/ETSP_ROMS_output.mat'];
load(fname);
if opt.rest_sio3
    inputs.restore.sio3_cout = nanmean(roms.sio3,2);
end
if opt.rest_fe
    inputs.restore.fe_cout   = nanmean(roms.fe,2);
end
if opt.rest_po4
    inputs.restore.po4_cout  = nanmean(roms.po4,2);
end
if opt.rest_no3
    inputs.restore.no3_cout = nanmean(roms.no3,2);
end
