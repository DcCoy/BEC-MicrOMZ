%-----------------------------------------------------------------------
%  Initializes particulate terms
%-----------------------------------------------------------------------
%  Code is based on Ballast model from Armstrong et al. 2000
%
%  July 2002, length scale for excess POC and bSI modified by temperature
%  Value given here is at Tref of 30 deg. C, JKM
%
%   diss       dissolution length for soft subclass
%   gamma      fraction of production -> hard subclass
%   mass       mass of 1e6 base units in kg    (WAS: 1e9 base units in g)
%   rho        QA mass ratio of POC to this particle class
%
%   Base units:
%     POC:        mmol C      (WAS: nmol C)
%     P_CaCO3:    mmol CaCO3  (WAS: nmol CaCO3)
%     P_SiO2:     mmol SiO2   (WAS: nmol SiO2)
%     dust:       kg dust     (WAS: g dust)
%     P_iron:     mmol Fe     (WAS: nmol Fe)
%
%  Units of fluxes:
%     sflux_in:    incoming flux of soft subclass (base units/m^2/sec)
%     hflux_in:    incoming flux of hard subclass (base units/m^2/sec)
%     prod:        production term (base units/m^3/sec)
%     sflux_out:   outgoing flux of soft subclass (base units/m^2/sec)
%     hflux_out:   outgoing flux of hard subclass (base units/m^2/sec)
%     remin:       remineralization term (base units/m^3/sec)
%    NOTE: Area/volume units were cm^2 and cm^3%
%-----------------------------------------------------------------------

POC_diss      = params.bec.POC_diss;   % diss. length (m), modified by TEMP
POC_mass      = params.bec.POC_mass;   % molecular weight of POC
POC_gamma     = params.bec.POC_gamma;  % fraction of production to hard class

P_CaCO3_diss  = params.bec.CaCO3_diss;                          % diss. length (m)
P_CaCO3_gamma = params.bec.CaCO3_gamma;                         % prod frac -> hard subclass
P_CaCO3_mass  = params.bec.CaCO3_mass;                          % molecular weight of CaCO3
P_CaCO3_rho   = params.bec.CaCO3_rho * P_CaCO3_mass / POC_mass; % QA mass ratio for CaCO3

P_SiO2_diss   = params.bec.SiO2_diss;                         % diss. length (m), modified by TEMP
P_SiO2_gamma  = params.bec.SiO2_gamma;                        % prod frac -> hard subclass
P_SiO2_mass   = params.bec.SiO2_mass;                         % molecular weight of SiO2
P_SiO2_rho    = params.bec.SiO2_rho * P_SiO2_mass / POC_mass; % QA mass ratio for SiO2

dust_diss     = params.bec.dust_diss;                       % diss. length (m) (DL: changed from cm)
dust_gamma    = params.bec.dust_gamma;                      % prod frac -> hard subclass
dust_mass     = params.bec.dust_mass;                       % base units are already kg
dust_rho      = params.bec.dust_rho * dust_mass / POC_mass; % QA mass ratio for dust

P_iron_gamma  = 0;              % prod frac -> hard subclass

%-----------------------------------------------------------------------
%  Set incoming fluxes
%-----------------------------------------------------------------------

P_CaCO3_sflux_out = 0;
P_CaCO3_hflux_out = 0;
P_SiO2_sflux_out  = 0;
P_SiO2_hflux_out  = 0;

dust_sflux_out = (1 - dust_gamma) * frc.dust;
dust_hflux_out = dust_gamma * frc.dust;

P_iron_sflux_out = 0;
P_iron_hflux_out = 0;

%-----------------------------------------------------------------------
%  Hard POC is QA flux and soft POC is excess POC.
%
%  Note (MF):
%  These names are convenient given the mineral associated soft and
%  hard components, but rather confusing when reading about the
%  particulate/ballast models in Armstrong et al 2002, and Lima et
%  al 2014 (there is no such thing as POC hard, POC_gamma
%  doesnt exist here or is set to c0 in CESM BEC). Soft/Hard actually
%  applies only to the mineral associated components. Think about
%  POC_hflux as being the QA component of POC and POC_sflux the excess POC.
%-----------------------------------------------------------------------

POC_sflux_out = 0;
POC_hflux_out = 0;

%-----------------------------------------------------------------------
%  Compute initial QA(dust) POC flux deficit.
%-----------------------------------------------------------------------

QA_dust_def = dust_rho * (dust_sflux_out + dust_hflux_out);

