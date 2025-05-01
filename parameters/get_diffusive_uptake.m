function params = get_diffusive_uptake(params)
% Calculates diffusive uptake of O2 across cell membrane
% used in place of Michaelis-Menten uptake kinetics
% since O2 is a small, non-polar molecule that does not
% require active enzyme-mediated transport into cells

% Diffusion coefficient of O2 @ 12C, 35PSU, 50bar 
D_o2  = 1.5776e-5; % (cm2/s to m2/s) 
D_n2o = D_o2*1.0049; 

% Convert to m2/s
D_o2  = D_o2/1e4;
D_n2o = D_n2o/1e4;

% Volume of cells (um^3 / cell)
vol_aer = 0.05;    % Based on SAR11
vol_aer = pi*((0.16*0.5)^2)*0.6;
vol_nar = vol_aer; % Based on SAR11
vol_nai = vol_aer; % Based on SAR11
vol_nao = vol_aer; % Based on SAR11
vol_nir = vol_aer; % Based on SAR11
vol_nio = vol_aer; % Based on SAR11
vol_nos = vol_aer; % Based on SAR11
vol_aoa = pi*((0.2*0.5)^2)*0.8;
vol_nob = pi*((0.3*0.5)^2)*3;
vol_aox = (4/3)*pi*((0.8*0.5)^3);

% Diameter of cells (um)
diam_aer = 2*(((3*vol_aer)/(4*pi))^(1/3));
diam_nar = 2*(((3*vol_nar)/(4*pi))^(1/3));
diam_nai = 2*(((3*vol_nai)/(4*pi))^(1/3));
diam_nao = 2*(((3*vol_nao)/(4*pi))^(1/3));
diam_nir = 2*(((3*vol_nir)/(4*pi))^(1/3));
diam_nio = 2*(((3*vol_nio)/(4*pi))^(1/3));
diam_nos = 2*(((3*vol_nos)/(4*pi))^(1/3));
diam_aoa = 2*(((3*vol_aoa)/(4*pi))^(1/3));
diam_nob = 2*(((3*vol_nob)/(4*pi))^(1/3));
diam_aox = 2*(((3*vol_aox)/(4*pi))^(1/3));

% Grams carbon per cell (assumes 0.1g DW/WW for all microbial types)
Ccell_aer = 0.1 * (12*params.aer.CN / (12*params.aer.CN + 7   + 16*2    + 14)) / (1e12 / vol_aer);
Ccell_nar = 0.1 * (12*params.nar.CN / (12*params.nar.CN + 7   + 16*2    + 14)) / (1e12 / vol_nar);
Ccell_nai = 0.1 * (12*params.nai.CN / (12*params.nai.CN + 7   + 16*2    + 14)) / (1e12 / vol_nai);
Ccell_nao = 0.1 * (12*params.nao.CN / (12*params.nao.CN + 7   + 16*2    + 14)) / (1e12 / vol_nao);
Ccell_nir = 0.1 * (12*params.nir.CN / (12*params.nir.CN + 7   + 16*2    + 14)) / (1e12 / vol_nir);
Ccell_nio = 0.1 * (12*params.nio.CN / (12*params.nio.CN + 7   + 16*2    + 14)) / (1e12 / vol_nio);
Ccell_nos = 0.1 * (12*params.nos.CN / (12*params.nos.CN + 7   + 16*2    + 14)) / (1e12 / vol_nos);
Ccell_aoa = 0.1 * (12*params.aoa.CN / (12*params.aoa.CN + 7   + 16*2    + 14)) / (1e12 / vol_aoa);
Ccell_nob = 0.1 * (12*params.nob.CN / (12*params.nob.CN + 7   + 16*2    + 14)) / (1e12 / vol_nob);
Ccell_aox = 0.1 * (12*params.aox.CN / (12*params.aox.CN + 8.7 + 16*1.55 + 14)) / (1e12 / vol_aox);

% Cell quotas (mol C / um^3) scaled by White et al., 2019 estimate of SAR11
Qc_aer = Ccell_aer / vol_aer / 12.0 * (6.5e-15 / Ccell_aer); % normalise to 6.5 fg C measured by White et al 2019
Qc_nar = Ccell_nar / vol_nar / 12.0 * (6.5e-15 / Ccell_aer); % normalise to 6.5 fg C measured by White et al 2019
Qc_nai = Ccell_nai / vol_nai / 12.0 * (6.5e-15 / Ccell_aer); % normalise to 6.5 fg C measured by White et al 2019
Qc_nao = Ccell_nao / vol_nao / 12.0 * (6.5e-15 / Ccell_aer); % normalise to 6.5 fg C measured by White et al 2019
Qc_nir = Ccell_nir / vol_nir / 12.0 * (6.5e-15 / Ccell_aer); % normalise to 6.5 fg C measured by White et al 2019
Qc_nio = Ccell_nio / vol_nio / 12.0 * (6.5e-15 / Ccell_aer); % normalise to 6.5 fg C measured by White et al 2019
Qc_nos = Ccell_nos / vol_nos / 12.0 * (6.5e-15 / Ccell_aer); % normalise to 6.5 fg C measured by White et al 2019
Qc_aoa = Ccell_aoa / vol_aoa / 12.0 * (6.5e-15 / Ccell_aer); % normalise to 6.5 fg C measured by White et al 2019
Qc_nob = Ccell_nob / vol_nob / 12.0 * (6.5e-15 / Ccell_aer); % normalise to 6.5 fg C measured by White et al 2019
Qc_aox = Ccell_aox / vol_aox / 12.0 * (6.5e-15 / Ccell_aer); % normalise to 6.5 fg C measured by White et al 2019

% Convert to mol C / cell 
% (mol C / um3) * (um3 / cell) * (1000 mmol / mol)
Qc_aer = Qc_aer * vol_aer * 1e3;
Qc_nar = Qc_nar * vol_nar * 1e3;
Qc_nai = Qc_nai * vol_nai * 1e3;
Qc_nao = Qc_nao * vol_nao * 1e3;
Qc_nir = Qc_nir * vol_nir * 1e3;
Qc_nio = Qc_nio * vol_nio * 1e3;
Qc_nos = Qc_nos * vol_nos * 1e3;
Qc_aoa = Qc_aoa * vol_aoa * 1e3;
Qc_nob = Qc_nob * vol_nob * 1e3;
Qc_aox = Qc_aox * vol_aox * 1e3;

% Get diffusion coeffient for O2
params.aer.pcoef = (4*pi*D_o2*(diam_aer*1e-6/2)) / Qc_aer; % m3/s/mmol BC
params.nar.pcoef = (4*pi*D_o2*(diam_nar*1e-6/2)) / Qc_nar; % m3/s/mmol BC
params.nai.pcoef = (4*pi*D_o2*(diam_nai*1e-6/2)) / Qc_nai; % m3/s/mmol BC
params.nao.pcoef = (4*pi*D_o2*(diam_nao*1e-6/2)) / Qc_nao; % m3/s/mmol BC
params.nir.pcoef = (4*pi*D_o2*(diam_nir*1e-6/2)) / Qc_nir; % m3/s/mmol BC
params.nio.pcoef = (4*pi*D_o2*(diam_nio*1e-6/2)) / Qc_nio; % m3/s/mmol BC
params.nos.pcoef = (4*pi*D_o2*(diam_nos*1e-6/2)) / Qc_nos; % m3/s/mmol BC
params.aoa.pcoef = (4*pi*D_o2*(diam_aoa*1e-6/2)) / Qc_aoa; % m3/s/mmol BC
params.nob.pcoef = (4*pi*D_o2*(diam_nob*1e-6/2)) / Qc_nob; % m3/s/mmol BC
params.aox.pcoef = (4*pi*D_o2*(diam_aox*1e-6/2)) / Qc_aox; % m3/s/mmol BC

% Get diffusion coefficient for N2O
params.nos.pana = (4*pi*D_n2o*(diam_nos*1e-6/2)) / Qc_nos;

% Set 'aer' functional type to have equal pana (trick to always get aerobic metabolism)
params.aer.pana = params.aer.pcoef;

