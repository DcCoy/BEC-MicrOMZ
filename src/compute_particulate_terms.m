% -------------------------------------------
% Compute particulate terms away from surface
% -------------------------------------------
% Code is based on Ballast model from Armstrong et al. 2000
% Terms:
% QA_dust_def        incoming deficit in the QA (dust) POC flux
% dust_remin         remineralization of dust (base units/m^3/sec)
% sed_denitrif       sedimentary denitrification (mmolN/m^3/s)
% other_remin        sedimentary remin not due to oxic or denitrification
% flux_oxidated      carbon oxidation rate in sediment
% fesedflux          sedimentary Fe inputs
% scalelength       
% TfuncS             temperature scaling from soft POM remin (right now just applied to ballast)
% Tfunc_soft         temperature scaling from soft POM remin (Laufkoetter 2017)
% DECAY_Hard         scaling factor for dissolution of Hard Ballast
% DECAY_HardDust     scaling factor for dissolution of Hard dust
% poc_diss_loc       diss. length used (m)
% sio2_diss_loc      diss. length varies spatially with O2 (m)
% caco3_diss_loc    
% dust_diss_loc     
% decay_POC_E        scaling factor for dissolution of excess POC
% decay_SiO2         scaling factor for dissolution of SiO2
% decay_CaCO3        scaling factor for dissolution of CaCO3
% decay_dust         scaling factor for dissolution of dust
% poc_prod_avail     POC production available for excess POC flux (mmol/m^3/s)
% new_QA_dust_def    outgoing deficit in the QA(dust) POC flux
% flux, flux_alt     temp variables used to update sinking flux
% dz_loc, dzr_loc    Hz and its inverse at a particular i,j,k location
%-----------------------------------------------------------------------
% Parameters from original model 
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

%-----------------------------------------------------------------------
% Get current grid cell thickness 
%-----------------------------------------------------------------------
dz_loc  = grid.Hz(k);
dzr_loc = (1 ./ dz_loc); 

%-----------------------------------------------------------------------
%  incoming fluxes are outgoing fluxes from previous level
%-----------------------------------------------------------------------
P_CaCO3_sflux_in = P_CaCO3_sflux_out;
P_CaCO3_hflux_in = P_CaCO3_hflux_out;
P_SiO2_sflux_in  = P_SiO2_sflux_out;
P_SiO2_hflux_in  = P_SiO2_hflux_out;
dust_sflux_in    = dust_sflux_out;
dust_hflux_in    = dust_hflux_out;
POC_sflux_in     = POC_sflux_out;
POC_hflux_in     = POC_hflux_out;
P_iron_sflux_in  = P_iron_sflux_out;
P_iron_hflux_in  = P_iron_hflux_out;

%-----------------------------------------------------------------------
%  initialize loss to sediments = 0 and local copy of percent sed
%-----------------------------------------------------------------------
P_iron_sed_loss = 0;
POC_sed_loss = 0;
P_CaCO3_sed_loss = 0;
P_SiO2_sed_loss = 0;
dust_sed_loss = 0;
sed_denitrif = 0;
other_remin = 0;
fesedflux = 0;

%-----------------------------------------------------------------------
%  compute scalelength and decay factors
%-----------------------------------------------------------------------
if -grid.z_w(k) < params.bec.scalelen_z(1)
    scalelength = params.bec.scalelen_vals(1);
elseif -grid.z_w(k) >= params.bec.scalelen_z(end)
    scalelength = params.bec.scalelen_vals(end);
else
    for i = 2:length(params.bec.scalelen_vals)
        if -grid.z_w(k) < params.bec.scalelen_z(i)
            scalelength = params.bec.scalelen_vals(i-1) + ...
                (params.bec.scalelen_vals(i) - params.bec.scalelen_vals(i-1)) .* ...
                (-grid.z_w(k) - params.bec.scalelen_z(i-1)) ./ ...
                (params.bec.scalelen_z(i) - params.bec.scalelen_z(i-1));
        end
    end
end
DECAY_Hard     = exp(-grid.Hz(k) / 4.0e4);  
DECAY_HardDust = exp(-grid.Hz(k) / 1.2e5);  

%----------------------------------------------------------------------
% Get temperature dependent function(s) 
%-----------------------------------------------------------------------
TfuncS = params.bec.Q_10.^((phy.temp(k) - params.bec.Tref) / 10);
if opt.TDEP_REMIN
    Tfunc_soft = exp(params.bec.ktfunc_soft*(phy.temp(k) - params.bec.Tref));
end
          
%----------------------------------------------------------------------
% Initialize dissociation length scales 
%-----------------------------------------------------------------------
POC_diss_loc = POC_diss;
sio2_diss_loc = P_SiO2_diss;
caco3_diss_loc = P_CaCO3_diss;
dust_diss_loc = dust_diss;

%-----------------------------------------------------------------------
% increase POC diss length scale where O2 concentrations are low
%-----------------------------------------------------------------------
if opt.O2_REMIN
    if tr.o2(k) >= 5.0 & tr.o2(k) < 40.0
       POC_diss_loc = POC_diss*(1 + (params.bec.lowo2_remin_factor - 1) .* ...
           (40.0 - tr.o2(k)) / 35.0);
    elseif tr.o2(k) < 5.0
       POC_diss_loc = POC_diss * params.bec.lowo2_remin_factor;
    end
end

%-----------------------------------------------------------------------
%  apply scalelength factor to length scales
%-----------------------------------------------------------------------
if opt.TDEP_REMIN
    POC_diss_loc = scalelength * POC_diss_loc / Tfunc_soft;
    sio2_diss_loc = scalelength * sio2_diss_loc / Tfunc_soft;
    caco3_diss_loc = scalelength * caco3_diss_loc / Tfunc_soft;
    dust_diss_loc = scalelength * dust_diss_loc / Tfunc_soft;
else
    POC_diss_loc = scalelength * POC_diss_loc; 
    sio2_diss_loc = scalelength * sio2_diss_loc; 
    caco3_diss_loc = scalelength * caco3_diss_loc; 
    dust_diss_loc = scalelength * dust_diss_loc; 
end

%-----------------------------------------------------------------------
%  apply temperature dependence to sio2_diss length scale
%-----------------------------------------------------------------------
sio2_diss_loc = sio2_diss_loc / TfuncS;

%-----------------------------------------------------------------------
%  decay_POC_E and decay_SiO2 set locally, modified by O2
%-----------------------------------------------------------------------
decay_POC_E = exp(-dz_loc / POC_diss_loc);
decay_SiO2  = exp(-dz_loc / sio2_diss_loc);
decay_CaCO3 = exp(-dz_loc / caco3_diss_loc);
decay_dust  = exp(-dz_loc / dust_diss_loc);

%-----------------------------------------------------------------------
%  Set outgoing fluxes for non-iron pools.
%  The outoing fluxes for ballast materials are from the
%  solution of the coresponding continuous ODE across the model
%  level. The ODE has a constant source term and linear decay.
%  It is assumed that there is no sub-surface dust production.
%-----------------------------------------------------------------------
P_CaCO3_sflux_out = P_CaCO3_sflux_in * decay_CaCO3 + ...
    P_CaCO3_prod * ((1 - P_CaCO3_gamma) * (1 - decay_CaCO3) * caco3_diss_loc);
P_CaCO3_hflux_out = P_CaCO3_hflux_in * DECAY_Hard + ...
    P_CaCO3_prod * (P_CaCO3_gamma * dz_loc);
P_SiO2_sflux_out = P_SiO2_sflux_in * decay_SiO2 + ...
    P_SiO2_prod * ((1 - P_SiO2_gamma) * (1 - decay_SiO2) * sio2_diss_loc);
P_SiO2_hflux_out = P_SiO2_hflux_in * DECAY_Hard + ...
    P_SiO2_prod * (P_SiO2_gamma * dz_loc);
dust_sflux_out = dust_sflux_in * decay_dust;
dust_hflux_out = dust_hflux_in * DECAY_HardDust;

%-----------------------------------------------------------------------
%  Compute how much POC_PROD is available for deficit reduction
%  and excess POC flux after subtracting off fraction of non-dust
%  ballast production from net POC_PROD.
%-----------------------------------------------------------------------
poc_prod_avail = poc_prod - (P_CaCO3_rho * P_CaCO3_prod) - (P_SiO2_rho * P_SiO2_prod);

%-----------------------------------------------------------------------
%  Check for POC production bounds violations
%-----------------------------------------------------------------------
if (poc_prod_avail < 0)
    disp('subroutine compute_particulate_terms: non_dust ballast production exceeds POC production');
    disp(['k = ',num2str(k)]);
    disp(['poc_prod_avail: ',num2str(poc_prod_avail)]);
    disp(['poc_prod: ',num2str(poc_prod)]);
    disp(['P_CaCO3_rho*P_CaCO3_prod: ',num2str(P_CaCO3_rho * P_CaCO3_prod)])
    disp(['P_SiO2_rho * P_SiO2_prod: ',num2str(P_SiO2_rho * P_SiO2_prod)]);
    kill
end

%-----------------------------------------------------------------------
%  Compute 1st approximation to new QA_dust_def, the QA_dust
%  deficit leaving the cell (implicit). In the case of explicit 
%  sinking, new_QA_dust_def equals dust flux leaving the layer.
%  Ignore poc_prod_avail at this stage.
%-----------------------------------------------------------------------
if QA_dust_def > 0
    new_QA_dust_def = QA_dust_def * ...
        (dust_sflux_out + dust_hflux_out) ./ ... 
        (dust_sflux_in + dust_hflux_in);
else
    new_QA_dust_def = 0;
end

%-----------------------------------------------------------------------
%  Use poc_prod_avail to reduce new_QA_dust_def.
%-----------------------------------------------------------------------
if new_QA_dust_def > 0
    new_QA_dust_def = new_QA_dust_def - poc_prod_avail * dz_loc;
    if new_QA_dust_def < 0
       poc_prod_avail = -new_QA_dust_def * dzr_loc;
       new_QA_dust_def = 0;
    else
       poc_prod_avail = 0;
    end
end
QA_dust_def = new_QA_dust_def;

%-----------------------------------------------------------------------
%  Compute outgoing POC fluxes. QA POC flux is computing using
%  ballast fluxes and new_QA_dust_def. If no QA POC flux came in
%  and no production occured, then no QA POC flux goes out. This
%  shortcut is present to avoid roundoff cancellation errors from
%  the dust_rho * dust_flux_out - QA_dust_def computation.
%  Any poc_prod_avail still remaining goes into excess POC flux.
%-----------------------------------------------------------------------
if POC_hflux_in == 0 & poc_prod == 0
    POC_hflux_out = 0;
else
    POC_hflux_out = P_CaCO3_rho * ...
        (P_CaCO3_sflux_out + P_CaCO3_hflux_out) + ...
        P_SiO2_rho * ...
        (P_SiO2_sflux_out + P_SiO2_hflux_out) + ...
        dust_rho * ...
        (dust_sflux_out + dust_hflux_out) - ...
        new_QA_dust_def;
    POC_hflux_out = max(POC_hflux_out, 0);
end

POC_sflux_out = POC_sflux_in * decay_POC_E + ...
    poc_prod_avail * ((1 - decay_POC_E) * POC_diss);

%-----------------------------------------------------------------------
%  Compute remineralization terms. It is assumed that there is no
%  sub-surface dust production.
%-----------------------------------------------------------------------
P_CaCO3_remin = P_CaCO3_prod + ...
    ((P_CaCO3_sflux_in - P_CaCO3_sflux_out) + ...
    (P_CaCO3_hflux_in - P_CaCO3_hflux_out)) * dzr_loc;
P_SiO2_remin = P_SiO2_prod + ...
    ((P_SiO2_sflux_in - P_SiO2_sflux_out) + ...
    (P_SiO2_hflux_in - P_SiO2_hflux_out)) * dzr_loc;
poc_remin = poc_prod + ...
    ((POC_sflux_in - POC_sflux_out) + ...
    (POC_hflux_in - POC_hflux_out)) * dzr_loc;
dust_remin = ...
    ((dust_sflux_in - dust_sflux_out) + ...
    (dust_hflux_in - dust_hflux_out)) * dzr_loc;

%-----------------------------------------------------------------------
%  Compute iron remineralization and flux out.
%-----------------------------------------------------------------------
if (POC_sflux_in+POC_hflux_in) == 0
    P_iron_remin = (poc_remin * params.bec.Red_Fe_C);
else
    P_iron_remin = (poc_remin * ...
        (P_iron_sflux_in + P_iron_hflux_in) / (POC_sflux_in + POC_hflux_in));
end

% Why is this hard-coded?
P_iron_remin = P_iron_remin + P_iron_sflux_in * 1.5e-5;

P_iron_sflux_out = P_iron_sflux_in + dz_loc * ...
    ((1 - P_iron_gamma) * P_iron_prod - P_iron_remin);

if (P_iron_sflux_out < 0)
    P_iron_sflux_out = 0;
    P_iron_remin = P_iron_sflux_in * dzr_loc + (1 - P_iron_gamma) * P_iron_prod;
end

%-----------------------------------------------------------------------
%  Compute iron release from dust remin/dissolution
%
%  dust remin gDust = 0.035 / 55.847 * 1.0e9 = 626712.0 nmolFe
%                      gFe     molFe     nmolFe
%  Also add in Fe source from sediments if applicable to this cell.
%-----------------------------------------------------------------------
P_iron_remin = P_iron_remin + dust_remin * params.bec.dust_to_Fe;
P_iron_hflux_out = P_iron_hflux_in;

%-----------------------------------------------------------------------
%  Bottom Sediments Cell?
%  If so compute sedimentary burial and denitrification N losses.
%  Using empirical relations from Bohlen et al., 2012 (doi:10.1029/2011GB004198) for Sed Denitrification
%  other_remin estimates organic matter remineralized in the sediments
%      by the processes other than oxic remin and denitrification (SO4 and CO2,
%      etc..)
%      based on Soetaert et al., 1996, varies between 10% and 50%
%      0.4_r8 is a coefficient with units mmolC/cm2/yr sinking flux,
%      other_remin is 50% above this high flux value,
%      In special case where bottom O2 has been depleted to < 1.0 uM,
%               all sedimentary remin is due to DENITRIFICATION + other_remin
%  POC burial from Dunne et al. 2007 (doi:10.1029/2006GB002907), maximum of 80% burial efficiency imposed
%  Bsi preservation in sediments based on
%     Ragueneau et al. 2000 (doi:10.1016/S0921-8181(00)00052-7)
%  Calcite is preserved in sediments above the lysocline, dissolves below.
%       Here a constant depth is used for lysocline.
%-----------------------------------------------------------------------

if opt.SEDIMENTS & k == length(grid.z_r)
    flux = POC_sflux_out+POC_hflux_out;  % mmol C/m^2/s
    if flux > 0
        flux_alt = flux.*params.bec.spd; % convert to mmol C/m^2/day
        POC_sed_loss = flux * min(0.8, params.bec.POMbury .* ...
            (0.013 + 0.53 * flux_alt*flux_alt / (7.0 + flux_alt).^2));
        if POC_sed_loss > 0 & tr.o2(k) > 0 & tr.no3(k) > 0 & tr.no3(k) < params.bec.sed_denitrif_NO3_limit & (tr.o2(k) - tr.no3(k)) > 0
            sed_denitrif = dzr_loc*min(...
                flux.*(0.06+0.19*0.99.^(tr.o2(k)-tr.no3(k))),...
                flux-(POC_sed_loss*params.bec.denitrif_C_N*(tr.no3(k)/(tr.no3(k)+params.bec.sed_denitrif_NO3_limit))));
        else
            sed_denitrif = 0;
        end
        flux_alt = flux *1e-4 * params.bec.spd * 365; % convert to mmol C/cm^2/year
        other_remin = dzr_loc.*min(...
            min(0.1 + flux_alt,0.5)*(flux - POC_sed_loss),...
            (flux-POC_sed_loss-(sed_denitrif*dz_loc*params.bec.denitrif_C_N)));

        % if bottom water O2 is depleted, assume all remin is denitrif + other               
        if tr.o2(k) < 1
            other_remin = dzr_loc * ...
                (flux - POC_sed_loss - (sed_denitrif * dz_loc * params.bec.denitrif_C_N));
        end

        % fesedflux parametrization from Dale 2015
        % carbon oxidation rate in mmol m-2 s-1 
        flux_oxidated = flux - POC_sed_loss;

        % convert from mmol m-2 s-1 to mmol m-2 d-1
        flux_oxidated = flux_oxidated * 86400;

        % constant 170 in umol m-2 d-1, flux_oxidated in mmol m-2
        % day-1, o2 in muM
        if (flux_oxidated > 0)
            fesedflux = 170.0 * tanh(flux_oxidated/tr.o2(k));
        end
        % convert from umol m-2 d-1 to mmol m-2 s-1
        fesedflux = fesedflux*0.001/86400*1.0;
    end  % flux > 0

    flux = P_SiO2_sflux_out+P_SiO2_hflux_out;
    flux_alt = flux*params.bec.spd; % convert to mmol/m^2/day
    % first compute burial efficiency, then compute loss to sediments
    if (flux_alt > 1)
       P_SiO2_sed_loss = 0.3 * flux_alt - 0.06;
    else
       P_SiO2_sed_loss = 0.04;
    end
    P_SiO2_sed_loss = flux * params.bec.BSIbury * P_SiO2_sed_loss;

    flux = 0; % set flux back to zero
    if (-grid.z_w(k) < 3300.0)
       flux = P_CaCO3_sflux_out + P_CaCO3_hflux_out;
       P_CaCO3_sed_loss = flux;
    end

    %----------------------------------------------------------------------------------
    %  Update sinking fluxes and remin fluxes, accounting for sediments.
    %  flux used to hold sinking fluxes before update.
    %----------------------------------------------------------------------------------

    flux = P_CaCO3_sflux_out + P_CaCO3_hflux_out;
    if (flux > 0)
        P_CaCO3_remin = P_CaCO3_remin + ...
            ((flux - P_CaCO3_sed_loss) * dzr_loc);
    end

    flux = P_SiO2_sflux_out + P_SiO2_hflux_out;
    if (flux > 0)
        P_SiO2_remin = P_SiO2_remin + ((flux - P_SiO2_sed_loss) * dzr_loc);
    end

    flux = POC_sflux_out + POC_hflux_out;
    if (flux > 0)
        poc_remin = poc_remin + ((flux - POC_sed_loss) * dzr_loc);
    end

    %-----------------------------------------------------------------------
    %   Remove all Piron and dust that hits bottom, sedimentary Fe source 
    %        accounted for by fesedflux elsewhere.
    %-----------------------------------------------------------------------

    flux = (P_iron_sflux_out + P_iron_hflux_out);
    if (flux > 0)
       P_iron_sed_loss = flux;
    end
    dust_sed_loss = dust_sflux_out + dust_hflux_out;
 
    %-----------------------------------------------------------------------
    %   Bottom layer: set all outgoing fluxes to 0.0
    %-----------------------------------------------------------------------
    P_CaCO3_sflux_out = 0;
    P_CaCO3_hflux_out = 0;

    P_SiO2_sflux_out = 0;
    P_SiO2_hflux_out = 0;

    dust_sflux_out = 0;
    dust_hflux_out = 0;

    POC_sflux_out = 0;
    POC_hflux_out = 0;

    P_iron_sflux_out = 0;
    P_iron_hflux_out = 0;

    P_iron_remin = P_iron_remin + (fesedflux * dzr_loc);
end

