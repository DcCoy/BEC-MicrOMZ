function params = get_BEC_params; 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default BGC parameters for the BEC model
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: If you want to override BEC parameters, do so in BEC_overrides.m 

% General constants (time, simple conversions) 
params.bec.spd       = 86400;                    % seconds per day
params.bec.dps       = (1/params.bec.spd);       % days per second
params.bec.yps       = (1/(365*params.bec.spd)); % years per second
params.bec.mpercm    = 0.01;                     % meters per cm
params.bec.t0_kelvin = 273.16;                   % 0C in Kelvin
params.bec.rho0      = 1027.4;                   % Bousinesq reference density (kg/m3)
params.bec.Cp        = 3985;                     % Specific heat of seawater (J/kgK)

% Air-sea gas exchange constants 
% (others are stored directly in functions/*satu.m and functions/cschmidt_*.m)
params.bec.a    = 0.31 .* (2.7778e-6); % Piston velocity conversion coefficient (cm/hr to m/s)
params.bec.xn2o = 300.0 .* 1e-9;       % Atmopheric concentration (ppb)
params.bec.xn2  = 0.78084;             % Molar ratio of N2 in air

% Misc. Iron constants
params.bec.fe_scavenge_thres1 = 0.8e-3;             % Upper thres. for Fe scavenging (mmol/m^3)
params.bec.dust_fescav_scale  = 1.0e10;             % Dust scavenging scale factor (was 1e9 in CESM)
params.bec.fe_max_scale2      = 1200.0;             % Unitless scaling coeff.
params.bec.dust_to_Fe         = 0.035/55.847*1.0e6; % Dust conversion to Fe 

% Partioning of phytoplankton growth, grazing, and losses
params.bec.caco3_poc_min       = 0.4;   % Minimum proportionality between QCaCO3 and grazing losses to POC (mmol C/mmol CaCO3)
params.bec.spc_poc_fac         = 0.14;  % Small phyto grazing factor (1/mmolC)
params.bec.f_graze_sp_poc_lim  = 0.36;  % 
params.bec.f_photosp_CaCO3     = 0.4;   % Proportionality between small phyto production and CaCO3 production
params.bec.f_graze_CaCO3_remin = 0.33;  % Fraction of spCaCO3 grazing which is remin
params.bec.f_graze_si_remin    = 0.35;  % Fraction of diatom Si grazing which is remin

% N fixation 
params.bec.r_Nfix_photo = 1.25; % N fixation relative to C fixation (non-dim)

% Fixed stoichiometric ratios
params.bec.Q          = 0.137;   % N/C ratio (mmol/mmol) of phyto & zoo
params.bec.Qp_zoo_pom = 0.00855; % P/C ratio (mmol/mmol) zoo & pom
params.bec.Qfe_zoo    = 3.0e-6;  % Zooplankton fe/C ratio
params.bec.gQsi_0     = 0.137;   % Initial Si/C ratio
params.bec.gQsi_max   = 0.8;     % Max Si/C ratio
params.bec.gQsi_min   = 0.0429;  % Min Si/C ratio
params.bec.QCaCO3_max = 0.4;     % Max QCaCO3

% Redfield ratios
% * assumes OM is C117H297O85N16P
params.bec.Red_D_C_P       = 117.0;                                        % carbon:phosphorus
params.bec.Red_D_N_P       =  16.0;                                        % nitrogen:phosphorus
params.bec.Red_D_O2_P      = 170.0;                                        % oxygen:phosphorus
params.bec.Remin_D_O2_P    = 138.0;                                        % oxygen:phosphorus
params.bec.Red_P_C_P       = params.bec.Red_D_C_P;                         % carbon:phosphorus
params.bec.Red_D_C_N       = params.bec.Red_D_C_P/params.bec.Red_D_N_P;    % carbon:nitrogen
params.bec.Red_P_C_N       = params.bec.Red_D_C_N;                         % carbon:nitrogen
params.bec.Red_D_C_O2      = params.bec.Red_D_C_P/params.bec.Red_D_O2_P;   % carbon:oxygen for HNO3 uptake*
params.bec.Remin_D_C_O2    = params.bec.Red_D_C_P/params.bec.Remin_D_O2_P; % carbon:oxygen for NH3 uptake* 
params.bec.Red_P_C_O2      = params.bec.Red_D_C_O2;                        % carbon:oxygen
params.bec.Red_Fe_C        = 3.0e-6;                                       % iron:carbon
params.bec.Red_D_C_O2_diaz = params.bec.Red_D_C_P/150.0;                   % carbon:oxygen for diazotrophs
params.bec.Red_D_C_O2_NO2V = params.bec.Red_D_C_P/162.0;                   % carbon:oxygen for HNO2 uptake*

% Denitrification ratios
params.bec.denitrif_C_N           = params.bec.Red_D_C_P/136.0; % carbon:nitrogen of denitrification
params.bec.denitrif_NO3_C         = 276.0 / 117.0;              % carbon:nitrogen of denitrif1 
params.bec.denitrif_NO2_C         = 276.0 / 117.0;              % carbon:nitrogen of denitrif2
params.bec.denitrif_N2O_C         = 276.0 / 117.0;              % carbon:nitrogen of denitrif3
params.bec.denitrif_NO3_limit     = 5.0;                        % threshold for reducing water column denitrification
params.bec.sed_denitrif_NO3_limit = 5.0;                        % threshold for reducing sediment denitrification

% Loss term threshold parameters and Chl:C ratios 
params.bec.thres_z1          = 100.0; % Threshold = C_loss_thres for z shallower than this (m)
params.bec.thres_z2          = 150.0; % Threshold = 0 for z deeper than this (m)
params.bec.CaCO3_temp_thres1 = 6.0;   % Upper temp threshold for CaCO3 prod
params.bec.CaCO3_temp_thres2 = -2.0;  % Lower temp threshold
params.bec.CaCO3_sp_thres    = 2.5;   % Bloom condition thres (mmolC/m3)

% PAR
params.bec.f_qsw_par      = 0.45;  % Fraction of incoming shortwave assumed to be PAR
params.bec.PAR_thres_pChl = 1e-10; % Threshold for PAR used in computation fo pChl

% Temperature parameters
params.bec.Tref = 30;        % Reference temperature (C)
params.bec.Q_10 = 1.7;       % Factor for temperature dependence (non-dim)
params.bec.Q_10_diat = 1.55; % Factor for diatom temperature dependence (non-dim) 

% BEC Constants used below
params.bec.epsC            = 1e-8;    % small C concentration (mmol C / m3)
params.bec.epsN            = 1e-30;   % small N concentration (mmol N / m3)
params.bec.epsTinv         = 3.17e-8; % small inverse time scale (1/year) (1/sec) 
params.bec.epsnondim       = 1e-6;    % small non-dimensional number
params.bec.Q               = 0.137;   % N/C ratio of phyto & zoo
params.bec.gQsi_max        = 0.8;     % Max Si/C ratio
params.bec.gQsi_min        = 0.0429;  % Min Si/C ratio
params.bec.QCaCO3_max      = 0.4;     % Max QCaCO3 
params.bec.cks             = 9;       % Constant used in Fe quota modification
params.bec.cksi            = 5;       % Constant used in Si quota modification
params.bec.par_thresh_pChl = 1e-10;   % Threshold for PAR used in pChl computation
params.bec.f_prod_caco3    = 0.055;   % Fraction of sp production as CaCO3 

% Other BEC constants
params.bec.Fe_bioavail        = 1.0;                   %
params.bec.o2_min             = 1.0;                   %
params.bec.lowo2_remin_factor = 3.3;                   %
params.bec.o2_min_delta       = 2.0;                   %
params.bec.kappa_nitrif       = 0.06 * params.bec.dps; %
params.bec.nitrif_par_lim     = 1.0;                   %
params.bec.z_mort_0           = 0.1 * params.bec.dps;  %
params.bec.z_mort2_0          = 0.4 * params.bec.dps;  %
params.bec.labile_ratio       = 0.85;                  % Fraction of loss instantly remineralized
params.bec.POMbury            = 1.0;                   %
params.bec.ktfunc_soft        = 0.055;                 %

% ...more constants
params.bec.epsC      = 1.00e-8; % small C concentration (mmol C/m^3)
params.bec.epsTinv   = 3.17e-8; % small inverse time scale (1/year) (1/sec)
params.bec.epsnondim = 1.00e-6; % small non-dimensional number (non-dim)
params.bec.cks       = 9.0;     % constant used in Fe quota modification
params.bec.cksi      = 5.0;     % constant used in Si quota modification

% NitrOMZ parameters
params.bec.kao           = 0.0500 * params.bec.dps; % max ammox (1/s)
params.bec.kno           = 0.0500 * params.bec.dps; % max nitrox rate (1/s)
params.bec.koxic         = 0.08 * params.bec.dps;   % max aerobic remin (1/s)
params.bec.kden1         = 0.0160 * params.bec.dps; % max denitrif1 (1/s)
params.bec.kden2         = 0.008 * params.bec.dps;  % max denitrif2 (1/s)
params.bec.kden3         = 0.0496 * params.bec.dps; % max denitrif3 (1/s)
params.bec.kax           = 0.441 * params.bec.dps;  % max anammox (1/s)
params.bec.ko2_ao        = 0.333;                   % half saturation for o2 uptake during ammox (mmol m-3)
params.bec.knh4_ao       = 0.305;                   % half saturation for nh4 uptake during ammox (mmol m-3)
params.bec.ko2_no        = 0.778;                   % half saturation for o2 uptake during nitrox (mmol m-3)
params.bec.kno2_no       = 0.509;                   % half sautration for no2 uptake during nitrox (mmol m-3)
params.bec.kno3_den1     = 1.0;                     % half saturation for no3 uptake during denitrif1 (mmol m-3) 
params.bec.kno2_den2     = 0.01;                    % half saturation for no2 uptake during denitrif2 (mmol m-3) 
params.bec.kn2o_den3     = 0.159;                   % half saturation for n2o uptake during denitrif3 (mmol m-3)
params.bec.ko2_oxic      = 1.0;                     % half sautration for o2 uptake during aerobic remin (mmol m-3)
params.bec.knh4_ax       = 1.0;                     % half saturation for nh4 uptake during anammox (mmol m-3) 
params.bec.kno2_ax       = 1.0;                     % half saturation for no2 uptake during anammox (mmol m-3)
params.bec.ko2_den1      = 6.0;                     % O2 inhibition for denitrif1 (mmol m-3)
params.bec.ko2_den2      = 2.3;                     % O2 inhibition for denitrif2 (mmol m-3)
params.bec.ko2_den3      = 0.506;                   % O2 inhibition for denitrif2 (mmol m-3)
params.bec.ko2_ax        = 6.0;                     % O2 inhibition for anammox (mmol m-3)
params.bec.r_no2tonh4_ax = 1.00;                    % ratio of no2 to nh4 uptake during anammox  
params.bec.n2o_ji_a      = 0.3;                     % O2-dependent N2O production from ammox (mmol O2)
params.bec.n2o_ji_b      = 0.1;                     % background N2O production from ammox 

% Constants for particulate terms
params.bec.BSIbury           = 0.65;   % x1 default
params.bec.Fe_scavenge_rate0 = 2.5;    % test initial scavenging rate 3 times higher (previous default : 3.0)
params.bec.f_prod_sp_CaCO3   = 0.055;  % x1 default
params.bec.POC_diss          = 105.0;  % dissociation length scale
params.bec.POC_mass          = 12.01;  % molecular weight of POC
params.bec.POC_gamma         = 0;      % fraction of production to hard class
params.bec.SiO2_diss         = 300.0;  % dissociation length scale
params.bec.SiO2_gamma        = 0.03;   % fraction of production to hard class
params.bec.SiO2_mass         = 60.08;  % molecular weight of SiO2
params.bec.SiO2_rho          = 0.05;   % QA mass ratio
params.bec.CaCO3_diss        = 180.0;  % dissociation length scale
params.bec.CaCO3_gamma       = 0.30;   % fraction of production to hard class
params.bec.CaCO3_mass        = 100.09; % molecular weight of CaCO3
params.bec.CaCO3_rho         = 0.05;   % QA mass ratio
params.bec.dust_diss         = 200;    % dissociation length scale
params.bec.dust_gamma        = 0.97;   % fraction of production to hard class
params.bec.dust_mass         = 1e6;    % base units are already in kg
params.bec.dust_rho          = 0.05;   % QA mass ratio
params.bec.scalelen_z        = [100 300 600 900]; 
params.bec.scalelen_vals     = [1.0 3.5 7.5 8.5];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Phytoplankton %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Small phytoplankton parameters 
params.sp.PCref         = 5 .* params.bec.dps;    % Max growth (s-1)
params.sp.K_no3         = 0.02;                   % Half-saturation constant for NO3 uptake (mmol N)
params.sp.K_no2         = 0.03;                   % Half-saturation constant for NO2 uptake (mmol N)
params.sp.K_nh4         = 0.0025;                 % Half-saturation constant for NH4 uptake (mmol N)
params.sp.K_fe          = 0.025e-3;               % Half-saturation constant for Fe uptake (mmol Fe)
params.sp.K_po4         = 0.0075;                 % Half-saturation constant for PO4 uptake (mmol PO4)
params.sp.K_dop         = 0.22;                   % Half-saturation constant for DOP uptake (mmol DOP)
params.sp.K_sio3        = 0.0;                    % Half-saturation constant for SiO3 uptake (mmol SiO3);
params.sp.K_I           = 10;                     % Light half-saturation constant (W m-2)
params.sp.e_o           = 150/16;                 % Production of O2 (mmol O2 / mmol N)    
params.sp.m_l           = 0.1 .* params.bec.dps;  % Linear mortality (uM BN-1 s-1) 
params.sp.m_q           = 0.01 .* params.bec.dps; % Quadratic mortality (uM BN-1 s-1)
params.sp.temp_thresh   = -10;                    % Temperature threshold
params.sp.loss_thresh   = 0.02;                   %
params.sp.loss_thresh2  = 0.0;                    %
params.sp.Qp            = 0.00855;                % P/C ratio
params.sp.gQfe_0        = 33e-6;                  %
params.sp.alphaPI       = 0.4 .* params.bec.dps;  % Chl specific initial slow of P_I curve ((mmol C / m2) / (mg Chl W sec))
params.sp.gQfe_min      = 2.7e-6;                 % Initial and minimum Fe/C ratios
params.sp.thetaN_max    = 2.5;                    % Max thetaN (Chl/N) (mg Chl/mmol N)
params.sp.Nfixer        = false;                  % Perform N fixation?
params.sp.imp_calcifier = true;                   % Implicit CaCO3 production?
params.sp.exp_calcifier = false;                  % Explicit CaCO3 production?
params.sp.agg_rate_max  = 0.5 .* params.bec.dps;  % Maximum aggregation rate (1/s)
params.sp.agg_rate_min  = 0.01 .* params.bec.dps; % Minimum aggregation rate (1/s)
params.sp.z_grz         = 1.2;                    % Grazing coefficient (mmol C / m3)
params.sp.graze_zoo     = 0.3;                    % Grazing fraction to zoo biomass (%)
params.sp.graze_poc     = 0.0;                    % Grazing fraction to poc (%)
params.sp.graze_doc     = 0.06;                   % Grazing fraction to doc (%)
params.sp.loss_poc      = 0.0;                    % Fraction of loss routed to poc (%)
params.sp.z_umax_0      = 3.3 .* params.bec.dps;  % Max grazing rate (1/s)
params.sp.f_zoo_detr    = 0.1;                    % Fraction of zoo losses to detrital pool (%) 
params.sp.Si_ind        = false;                  % Source or sink of SiO3? 
params.sp.CaCO3_ind     = true;                   % Source or sink of CacO3?

% Diatom parameters
params.diat.PCref         = 5 .* params.bec.dps;    % Max growth (s-1)
params.diat.K_no3         = 0.6;                    % Half-saturation constant for NO3 uptake (mmol N)
params.diat.K_no2         = 0.24;                   % Half-saturation constant for NO2 uptake (mmol N)
params.diat.K_nh4         = 0.02;                   % Half-saturation constant for NH4 uptake (mmol N)
params.diat.K_fe          = 0.05e-3;                % Half-saturation constant for Fe uptake (mmol Fe)
params.diat.K_po4         = 0.06;                   % Half-saturation constant for PO4 uptake (mmol PO4)
params.diat.K_dop         = 0.6;                    % Half-saturation constant for DOP uptake (mmol DOP)
params.diat.K_sio3        = 0.6;                    % Half-saturation constant for SiO3 uptake (mmol SiO3);
params.diat.K_I           = 10;                     % Light half-saturation constant (W m-2)
params.diat.e_o           = 150/16;                 % Production of O2 (mmol O2 / mmol N)    
params.diat.m_l           = 0.1 .* params.bec.dps;  % Linear mortality (uM BN-1 s-1) 
params.diat.m_q           = 0.01 .* params.bec.dps; % Quadratic mortality (uM BN-1 s-1)
params.diat.temp_thresh   = -10;                    % Temperature threshold
params.diat.loss_thresh   = 0.02;                   %
params.diat.loss_thresh2  = 0.0;                    %
params.diat.Qp            = 0.00855;               % P/C ratio 
params.diat.gQfe_0        = 33e-6;                  %
params.diat.alphaPI       = 0.31 .* params.bec.dps; % Chl specific initial slow of P_I curve ((mmol C / m2) / (mg Chl W sec))
params.diat.gQfe_min      = 2.7e-6;                 % Initial and minimum Fe/C ratios
params.diat.thetaN_max    = 4.0;                    % Max thetaN (Chl/N) (mg Chl/mmol N)
params.diat.Nfixer        = false;                  % Perform N fixation?
params.diat.imp_calcifier = false;                  % Implicit CaCO3 production?
params.diat.exp_calcifier = false;                  % Explicit CaCO3 production?
params.diat.agg_rate_max  = 0.5 .* params.bec.dps;  % Maximum aggregation rate (1/s)
params.diat.agg_rate_min  = 0.01 .* params.bec.dps; % Minimum aggregation rate (1/s)
params.diat.z_grz         = 1.2;                    % Grazing coefficient (mmol C / m3)
params.diat.graze_zoo     = 0.25;                   % Grazing fraction to zoo biomass (%)
params.diat.graze_poc     = 0.4;                    % Grazing fraction to poc (%)
params.diat.graze_doc     = 0.06;                   % Grazing fraction to doc (%)
params.diat.loss_poc      = 0.0;                    % Fraction of loss routed to poc (%)
params.diat.z_umax_0      = 3.05 .* params.bec.dps; % Max grazing rate (1/s)
params.diat.f_zoo_detr    = 0.2;                    % Fraction of zoo losses to detrital pool (%) 
params.diat.Si_ind        = true;                   % Source or sink of SiO3? 
params.diat.CaCO3_ind     = false;                  % Source or sink of CacO3?

% Diazotroph parameters
params.diaz.PCref         = 2.5 .* params.bec.dps;  % Max growth (s-1)
params.diaz.K_no3         = 2.0;                    % Half-saturation constant for NO3 uptake (mmol N)
params.diaz.K_no2         = 0.8;                    % Half-saturation constant for NO2 uptake (mmol N)
params.diaz.K_nh4         = 0.0667;                 % Half-saturation constant for NH4 uptake (mmol N)
params.diaz.K_fe          = 0.025e-3;               % Half-saturation constant for Fe uptake (mmol Fe)
params.diaz.K_po4         = 0.015;                  % Half-saturation constant for PO4 uptake (mmol PO4)
params.diaz.K_dop         = 0.05;                   % Half-saturation constant for DOP uptake (mmol DOP)
params.diaz.K_sio3        = 0.0;                    % Half-saturation constant for SiO3 uptake (mmol SiO3)
params.diaz.K_I           = 10;                     % Light half-saturation constant (W m-2)
params.diaz.e_o           = 150/16;                 % Production of O2 (mmol O2 / mmol N)    
params.diaz.m_l           = 0.1 .* params.bec.dps;  % Linear mortality (uM BN-1 s-1) 
params.diaz.m_q           = 0.01 .* params.bec.dps; % Quadratic mortality (uM BN-1 s-1)
params.diaz.temp_thresh   = 16;                     % Temperature threshold
params.diaz.loss_thresh   = 0.02;                   %
params.diaz.loss_thresh2  = 0.0;                    %
params.diaz.Qp            = 0.00855;                % P/C ratio 
params.diaz.gQfe_0        = 66e-6;                  %
params.diaz.alphaPI       = 0.33 .* params.bec.dps; % Chl specific initial slow of P_I curve ((mmol C / m2) / (mg Chl W sec))
params.diaz.gQfe_min      = 6.0e-6;                 % Initial and minimum Fe/C ratios
params.diaz.thetaN_max    = 2.5;                    % Max thetaN (Chl/N) (mg Chl/mmol N)
params.diaz.Nfixer        = true;                   % Perform N fixation?
params.diaz.imp_calcifier = false;                  % Implicit CaCO3 production?
params.diaz.exp_calcifier = false;                  % Explicit CaCO3 production?
params.diaz.agg_rate_max  = 0.5 .* params.bec.dps;  % Maximum aggregation rate (1/s)
params.diaz.agg_rate_min  = 0.01 .* params.bec.dps; % Minimum aggregation rate (1/s)
params.diaz.z_grz         = 1.2;                    % Grazing coefficient (mmol C / m3)
params.diaz.graze_zoo     = 0.3;                    % Grazing fraction to zoo biomass (%)
params.diaz.graze_poc     = 0.1;                    % Grazing fraction to poc (%)
params.diaz.graze_doc     = 0.06;                   % Grazing fraction to doc (%)
params.diaz.loss_poc      = 0.0;                    % Fraction of loss routed to poc (%)
params.diaz.z_umax_0      = 3.05 .* params.bec.dps; % Max grazing rate (1/s)
params.diaz.f_zoo_detr    = 0.1;                    % Fraction of zoo losses to detrital pool (%) 
params.diaz.Si_ind        = false;                  % Source or sink of SiO3? 
params.diaz.CaCO3_ind     = false;                  % Source or sink of CacO3?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Zooplankton %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Large zoo
params.zoo.m_l        = 0.1 .* params.bec.dps; % Linear mortality (uM C-1 s-1) 
params.zoo.m_q        = 0.4 .* params.bec.dps; % Quadratic mortality (uM C-1 s-1)
params.zoo.Y          = 0.5;                   % Fraction of digestion
params.zoo.bmin       = 0.06;                  % Minimum biomass where losses stop 

