% Script to extract initial conditions from ROMS files

% List initial condition file
fname = ['/data/project2/demccoy/ROMS_configs/peru_chile_0p1/ini/peru_chile_0p1_ini_hybrid_ini.nc'];
gname = ['/data/project2/demccoy/ROMS_configs/peru_chile_0p1/grid/peru_chile_0p1_grd.nc'];

% List coordinate to extract profile
xidx = 247; % 275W
yidx = 382; % -8S

% Load vertical coordinates
zeta    = ncread(fname,'zeta');
hc      = ncread(fname,'hc');
h       = ncread(gname,'h');
hc      = ncread(fname,'hc');
theta_s = ncread(fname,'theta_s'); 
theta_b = ncread(fname,'theta_b');
s_rho   = 42;

% Get depth
z_r = zlevs4(h,zeta,theta_s,theta_b,hc,s_rho,'r','new2012');
z_r = permute(z_r,[2 3 1]);
tmp.z_r = flipud(squeeze(z_r(xidx,yidx,:)));

% List variables to extract and their renamed counterparts
fvars = {'temp','salt','O2','PO4','NO3','SiO3','Fe','Alk','DIC','SPC','SPCHL','SPFE','SPCACO3',...
         'DIATC','DIATCHL','DIATSI','DIATFE','DIAZC','DIAZCHL','DIAZFE','ZOOC','DON','DONR','DOP','DOPR',...
         'DOFE','PIC'};
myvars = {'temp','salt','o2','po4','no3','sio3','fe','alk','dic','sp','sp_chl','sp_fe','sp_caco3',...
          'diat','diat_chl','diat_si','diat_fe','diaz','diaz_chl','diaz_fe','zoo','don','donr','dop','dopr',...
          'dofe','pic'};

% Load and restrict variables
for i = 1:length(fvars)
    tmp.(myvars{i}) = flipud(squeeze(ncread(fname,fvars{i},[xidx yidx 1 1],[1 1 inf inf])));
end

% Save
save('../data/ETSP_ROMS_restart.mat','tmp');
