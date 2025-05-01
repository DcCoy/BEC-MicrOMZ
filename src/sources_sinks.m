function [sms diag tr] =  sources_sinks(tr,phy,grid,opt,params,frc,diag_out); 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specifies the biogeochemical and ecosystem sources and sinks 
% The model is based off of the Biogeochemical Elemental Cycling (BEC) model from
% Moore et al., 2004 as implemented in ROMS-BEC via Deutsch et al., 2020
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For safety, reduces zero and negative variables to small value
vars = fields(tr);
epsn = 1e-24;
for i = 1:length(vars)
    for k = 1:length(grid.z_r)
        tr.(vars{i})(k) = max(epsn,tr.(vars{i})(k));
    end
end

% Set interior grid points (starting at surface)
for k = 1:length(grid.z_r)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Various k == 1 initializations and air-sea flux 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %----------------------------------------------------------
    % Initialize surface PAR and particulate terms 
    %----------------------------------------------------------
    if k == 1
        % Get par at surface based on forcing
        % NOTE: frc.swrad is already in the correct units
        % No need to multiply by rho0*Cp
        par_out = max(0, params.bec.f_qsw_par * frc.swrad);
        
        % Get initial particulate pools at surface cell
        % Mostly sets their pools to 0 at surface cell, but also feeds in dust forcing (if any)
        initialize_particulate_terms
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute PAR related quantities 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %----------------------------------------------------------
    % PAR exponentially attenuates at depth
    %----------------------------------------------------------
    % in = out, i.e. from layer above (k-1)
    par_in = par_out; 
    % Attentuation parameter, increased if phyto pool is large (shading)
    kpardz = max(tr.sp_chl(k) + tr.diat_chl(k) + tr.diaz_chl(k),0.02);
    if kpardz < 0.13224
        kpardz = 0.0919 * (kpardz^0.3536) * grid.Hz(k);
    else
        kpardz = 0.1131 * (kpardz^0.4562) * grid.Hz(k);
    end
    % Calculate PAR out (bottom of cell)
    par_out = par_in * exp(-kpardz);
    % Get PAR averaged over layer for BGC calculations
    par_lay = par_in * (1 - exp(-kpardz)) / kpardz;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Temperature dependent growth 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %----------------------------------------------------------
    % Tfunc is 0.2 --> 1.0 between 0 and 30 degC
    % Tfunc_diat is higher at lower temperatures (0.27 --> 1)
    %----------------------------------------------------------
    Tfunc = (params.bec.Q_10)^(((phy.temp(k) + params.bec.t0_kelvin) - ...
            (params.bec.Tref + params.bec.t0_kelvin)) / 10); 
    Tfunc_diat = (params.bec.Q_10_diat)^(((phy.temp(k) + params.bec.t0_kelvin) - ...
            (params.bec.Tref + params.bec.t0_kelvin)) / 10);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% Phytoplankton %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Begin loop over phytoplankton
    for i = 1:length(opt.auto_ind);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get stoichiometry of phytoplankton
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %----------------------------------------------------------
        % Determine local elemental ratios for new growth (new biomass)
        % Modify the initial Fe/C ratios under low ambient iron conditions
        % Modify the initial Si/C ratio under low ambient Si conditions
        %----------------------------------------------------------
        % Get local Chl:C ratio of autotroph
        thetaC(i) = tr.([opt.auto_ind{i},'_chl'])(k) / (tr.(opt.auto_ind{i})(k) + params.bec.epsC); 
        % Get local Fe:C ratio of autotroph
        Qfe(i)    = tr.([opt.auto_ind{i},'_fe'])(k) / (tr.(opt.auto_ind{i})(k) + params.bec.epsC);
        % If autotroph requires SiO3 ...
        if params.(opt.auto_ind{i}).K_sio3 > 0
            % Get local Si:C ratio of autotroph, restrict to maximum ratio to inhibit excess SiO3 uptake
            Qsi(i) = min((tr.sio3(k) / tr.(opt.auto_ind{i})(k) + params.bec.epsC), params.bec.gQsi_max); 
        else
            Qsi(i) = 0; % Set to 0 if autotroph does not require SiO3
        end
        % Get initial Fe:C quota
        gQfe(i) = params.(opt.auto_ind{i}).gQfe_0;
        % Modify quota if local Fe concentration is below autotroph-specific threshold
        if tr.fe(k) < (params.bec.cks * params.(opt.auto_ind{i}).K_fe);
            % Update quota, restrict to minimum to allow growth at low Fe 
            gQfe(i) = max((gQfe(i) * tr.fe(k) / (params.bec.cks * params.(opt.auto_ind{i}).K_fe)), ...
                params.(opt.auto_ind{i}).gQfe_min);
        end    
        % If autotroph requires SiO3 ...
        if params.(opt.auto_ind{i}).K_sio3 > 0
            % Get initial Si:C quota
            gQsi(i) = params.bec.gQsi_0;
            % Modify quota if local Fe concentration is below autotroph-specific threshold
            if tr.fe(k) < (params.bec.cksi * params.(opt.auto_ind{i}).K_sio3) 
                % If local Fe is available ...
                if tr.fe(k) > 0 
                    % ... and if local SiO3 is greater than autotroph-specific threshold
                    if tr.sio3(k) > (params.bec.cksi * params.(opt.auto_ind{i}).K_sio3) 
                        % Update quota, restrict to minimum to allow growth at low SiO3
                        gQsi(i) = min((gQsi(i) * params.bec.cksi * params.(opt.auto_ind{i}).K_fe / tr.fe(k)), ...
                            params.bec.gQsi_max);
                    end
                end
            end
            % If no local Fe is available, set quota to maximum
            if tr.fe(k) == 0
                gQsi(i) = gQsi_max;
            end
            % If local SiO3 is high, allow for increased SiO3 quota 
            if tr.sio3(k) > (params.bec.cksi * params.(opt.auto_ind{i}).K_sio3);
                gQsi(i) = max((gQsi(i) * tr.sio3(k) / (params.bec.cksi * params.(opt.auto_ind{i}).K_sio3)), ...
                    params.bec.gQsi_min);
            end
        end

        % Get CaCO3 quota for coccolithophores
        if params.(opt.auto_ind{i}).CaCO3_ind
            QCaCO3(i) = min(tr.([opt.auto_ind{i},'_caco3'])(k) / (tr.(opt.auto_ind{i})(k) + params.bec.epsC), ...
                params.bec.QCaCO3_max);
        else
            QCaCO3(i) = 0;
        end
    end % end phytoplankton loop

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Depth dependent mortality thresholds 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %----------------------------------------------------------
    % Calculate the loss threshold interpolation factor
    %----------------------------------------------------------
    % Prevent extinction at the surface
    if -grid.z_r(k) > params.bec.thres_z1
        if -grid.z_r(k) < params.bec.thres_z2
            f_loss_thres = (params.bec.thres_z2 + grid.z_r(k))/(params.bec.thres_z2 - params.bec.thres_z1);
        else
            f_loss_thres = 0;
        end
    % Allow for extinction deeper in the water column
    else
        f_loss_thres = 1;
    end

    % Begin loop over phytoplankton
    for i = 1:length(opt.auto_ind) 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Autotroph Pprime
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %----------------------------------------------------------
        % Compute Pprime for all autotrophs, used for loss terms
        % Pprime is zero when biomass is at or below threshold
        % When biomass is high, Pprime =~ carbon biomass
        % As Pprime goes to zero, phytoplankton losses are also zero,
        % meaning they cannot go extinct
        %----------------------------------------------------------
        if phy.temp(k) < params.(opt.auto_ind{i}).temp_thresh;
            C_loss_thresh = f_loss_thres * params.(opt.auto_ind{i}).loss_thresh2; 
        else
            C_loss_thresh = f_loss_thres * params.(opt.auto_ind{i}).loss_thresh;
        end
        Pprime(i) = max(0,tr.(opt.auto_ind{i})(k) - C_loss_thresh);
    end % end autotroph loop    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Nutrient limitation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %----------------------------------------------------------
    % Here, find the limiting nutrient between:
    % N, P, Fe, Silicate (diatoms only), and light
    % Limitation takes the form of Michaelis-Menten kinetics 
    %----------------------------------------------------------
    for i = 1:length(opt.auto_ind);
        % Individual DIN limitation
        % Coded such that VNO3 + VNO2 + VNH4 = value between [0,1]
        VNO3(i) = (tr.no3(k) / params.(opt.auto_ind{i}).K_no3 / (1 + ...
            (tr.no3(k) / params.(opt.auto_ind{i}).K_no3) + ... 
            (tr.no2(k) / params.(opt.auto_ind{i}).K_no2) + ... 
            (tr.nh4(k) / params.(opt.auto_ind{i}).K_nh4)));
        VNO2(i) = (tr.no2(k) / params.(opt.auto_ind{i}).K_no2 / (1 + ...
            (tr.no3(k) / params.(opt.auto_ind{i}).K_no3) + ... 
            (tr.no2(k) / params.(opt.auto_ind{i}).K_no2) + ... 
            (tr.nh4(k) / params.(opt.auto_ind{i}).K_nh4)));
        VNH4(i) = (tr.nh4(k) / params.(opt.auto_ind{i}).K_nh4 / (1 + ...
            (tr.no3(k) / params.(opt.auto_ind{i}).K_no3) + ... 
            (tr.no2(k) / params.(opt.auto_ind{i}).K_no2) + ... 
            (tr.nh4(k) / params.(opt.auto_ind{i}).K_nh4)));
        % Get total N limitation (ignored for diaz, who can fix N)
        % Value between [0,1]
        if params.(opt.auto_ind{i}).Nfixer
            VNtot(i) = 1;
        else
            VNtot(i) = VNO3(i) + VNO2(i) + VNH4(i);
        end
        
        % Iron
        % Value between [0,1]
        VFe(i) = tr.fe(k) / (tr.fe(k) + params.(opt.auto_ind{i}).K_fe); 

        % Find the limiting nutrient between Fe, DIN 
        f_nut = min(VNtot(i), VFe(i));

        % Phosphorous limitation 
        % Unlike nitrogen, autotrophs can use DOP for P quota
        % Coded such that VPO4 + VDOP = value between [0,1]
        VPO4(i) = (tr.po4(k) / params.(opt.auto_ind{i}).K_po4 / (1 + ...
            (tr.po4(k) / params.(opt.auto_ind{i}).K_po4) + ...
            (tr.dop(k) / params.(opt.auto_ind{i}).K_dop)));
        VDOP(i) = (tr.dop(k) / params.(opt.auto_ind{i}).K_dop / (1 + ...
            (tr.po4(k) / params.(opt.auto_ind{i}).K_po4) + ...
            (tr.dop(k) / params.(opt.auto_ind{i}).K_dop)));
        VPtot(i) = VPO4(i) + VDOP(i);

        % Update limiting nutrient between Fe/N and P
        f_nut = min(f_nut, VPtot(i));

        % Silicate limitation (for diatoms only)
        if params.(opt.auto_ind{i}).K_sio3 > 0
            VSiO3(i) = tr.sio3(k) / (tr.sio3(k) + params.(opt.auto_ind{i}).K_sio3);
            % Update limiting nutrient between Fe/N/P and SiO3
            f_nut = min(f_nut, VSiO3(i));
        else
            VSiO3(i) = 0;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Apply nutrient + light limitation 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %----------------------------------------------------------
        % PCref is the maximum growth rate
        % PCmax is PCref augmented by limiting functions:
        % ... nutrient (f_nut)
        % ... temperature (Tfun)
        % ... and light (light_lim)
        %----------------------------------------------------------
        
        % Apply nutrient limititation, modified by temperature
        if ~strcmp(opt.auto_ind{i},'diat');
            PCmax = params.(opt.auto_ind{i}).PCref * f_nut * Tfunc;
        else
            PCmax = params.(opt.auto_ind{i}).PCref * f_nut * Tfunc_diat;
        end
        % Prevent growth below temperature threshold
        if phy.temp(k) < params.(opt.auto_ind{i}).temp_thresh;
            PCmax = 0;
        end

        % Apply light limitation
        light_lim = 1 - exp((-1 * params.(opt.auto_ind{i}).alphaPI * thetaC(i) * par_lay) / (PCmax + params.bec.epsTinv)); 

        % Get max rate after applying all limiting terms
        PCphoto = PCmax * light_lim;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute biomass growth
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %----------------------------------------------------------
        % PCphoto is the max rate after limitation
        % Multiply by local biomass to get growth over timestep
        %----------------------------------------------------------

        % Get growth after limitation terms
        photoC(i) = PCphoto * tr.(opt.auto_ind{i})(k);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Nutrient uptake rates
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %----------------------------------------------------------
        % Calculate nutrient uptake based on C growth rates
        %----------------------------------------------------------
        % DIN uptake 
        NO3_V(i) = 0; % NO3 contribution to N limitation 
        NO2_V(i) = 0; % NO2 contribution to N limitation
        NH4_V(i) = 0; % NH4 contribution to N limitation
        VNC(i)   = 0; % 
        % If DIN is available ...
        if VNtot(i) > 0
            % Uptake DIN based on individual contributions to N limitation
            % ... and the growth rate (converted from C to N)
            if tr.no3(k) > 0
                NO3_V(i) = (VNO3(i) / VNtot(i)) * photoC(i) * params.bec.Q;
            end
            if tr.no3(k) > 0
                NO2_V(i) = (VNO2(i) / VNtot(i)) * photoC(i) * params.bec.Q; 
            end
            if tr.nh4(k) > 0
                NH4_V(i) = (VNH4(i) / VNtot(i)) * photoC(i) * params.bec.Q;
            end
            % Get rate (C N-1 s-1) 
            VNC(i) = PCphoto * params.bec.Q;
        end

        % Phosphorous
        PO4_V(i) = 0; % PO4 contribution to P limitation 
        DOP_V(i) = 0; % DOP contribution to P limitation
        % If P is available ...
        if VPtot(i) > 0
            % Uptake P based on individual contributions to P limitation
            % ... and the growth rate (converted from C to P) 
            PO4_V(i) = (VPO4(i) / VPtot(i)) * photoC(i) * params.(opt.auto_ind{i}).Qp;
            DOP_V(i) = (VDOP(i) / VPtot(i)) * photoC(i) * params.(opt.auto_ind{i}).Qp;
        end

        % Iron
        % Uptake based on iron quota of autotrophs (modified above at low Fe)
        photoFe(i) = photoC(i) * gQfe(i);

        % Silicate
        % Uptake based on SiO3 quota of autotrophs (modified above at low SiO3)
        if params.(opt.auto_ind{i}).K_sio3 > 0
            photoSi(i) = photoC(i) * gQsi(i);
        else
            photoSi(i) = 0;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Photoadaptation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %----------------------------------------------------------
        %  calculate pChl, (used in photoadapt., GD98)
        %  2.3 = max value of thetaN (Chl/N ratio) (mg Chl/mmol N)
        %  GD 98 Chl. synth. term
        %----------------------------------------------------------
        work1 = params.(opt.auto_ind{i}).alphaPI * thetaC(i) * par_lay; 
        if work1 > params.bec.par_thresh_pChl;
            pChl = params.(opt.auto_ind{i}).thetaN_max * PCphoto / work1;
            photoacc(i) = (pChl * VNC(i) / thetaC(i)) * tr.([opt.auto_ind{i},'_chl'])(k);
        else
            photoacc(i) = 0;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Implicit CaCO3 production 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %----------------------------------------------------------
        %  Parameterized as function of small phyto production:
        %  - decrease CaCO3 as function of nutrient limitation
        %  - decrease CaCO3 prod at low temperatures
        %  - increase CaCO3 prod under bloom conditions
        %  Maximum calcification rate is 40% of primary production.
        %  explicit CaCO3 Production, parametrized as function of photoC of PFT
        %  - think about: temp dependence, f_nut dependence, CaCO3:Corg
        %  - ratio (right now set to one)
        %----------------------------------------------------------
        caco3_prod(i) = 0;
        % Only take up CaCO3 if autotroph is a calcifier
        if params.(opt.auto_ind{i}).imp_calcifier
            % Get fraction of CaCO3 production rate from C growth rate 
            caco3_prod(i) = params.bec.f_prod_caco3 * photoC(i);
            % Limit production based on nutrient limitation
            caco3_prod(i) = caco3_prod(i) * f_nut;
            % Decrease production at cold temperatures
            if phy.temp(k) < params.bec.CaCO3_temp_thres1 
                caco3_prod(i) = caco3_prod(i) * max((phy.temp(k) - params.bec.CaCO3_temp_thres2),0) / ...
                    (params.bec.CaCO3_temp_thres1 - params.bec.CaCO3_temp_thres2);
            end
            % At high biomass (bloom), increase production 
            if tr.(opt.auto_ind{i}) > params.bec.CaCO3_sp_thres
                caco3_prod(i) = min((caco3_prod(i) * tr.(opt.auto_ind{i}) / params.bec.CaCO3_sp_thres),...
                    (params.bec.f_photosp_CaCO3 * photoC(i))); 
            end
        elseif params.(opt.auto_ind{i}).exp_calcifier % off for all autotrophs
            caco3_prod(i) = QCaCO3(i) * photoC(i);
        else
            caco3_prod(i) = 0;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Autotroph mortality / aggregation 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %----------------------------------------------------------
        % Mortality of the phytoplankton groups is based on a linear 
        % mortality (m_l, 1/s) and a quadratic mortality (m_q, 1/s) 
        % But is restricted to between agg_rate_max and min
        %----------------------------------------------------------
        % Linear mortality losses are reduced at colder temperatures
        if strcmp(opt.auto_ind{i},'diat') 
            auto_loss(i) = params.(opt.auto_ind{i}).m_l * Pprime(i) * Tfunc_diat;
        else
            auto_loss(i) = params.(opt.auto_ind{i}).m_l * Pprime(i) * Tfunc;
        end
        % Agg collects quadratic mortality, but is restricted to min/max
        auto_agg(i) = min(params.(opt.auto_ind{i}).agg_rate_max * Pprime(i),...
            params.(opt.auto_ind{i}).m_q * Pprime(i) * Pprime(i));
        auto_agg(i) = max(params.(opt.auto_ind{i}).agg_rate_min * Pprime(i),...
            auto_agg(i));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Autotroph grazing 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ---------------------------------------------------------------------------------
        % Grazing by zooplankton takes the simple form of:
        %    total_grazing = mu(zoo) * B(zoo) *  B(prey) / ( B(prey) + K )
        % Grazing losses to a specific phyto is calculated as a proportion of their biomass
        % relative to the total biomass of phytoplankton, such that:
        % diatom_grazing = (B(diat) / B(phyto)) * total_grazing
        % ---------------------------------------------------------------------------------
        % Retrieve the biomass of all prey
        work1 = 0;
        work1 = work1 + sum(Pprime);
        if opt.EXPLICIT_MICROBES
            % Also include chemoautotrophs as zoo prey
            for j = 1:length(opt.chemo_ind)
                work1 = work1 + (max(tr.(opt.chemo_ind{j})(k) - params.(opt.chemo_ind{j}).bmin,0));
            end
        end

        % Limit grazing at low O2
        z_o2lim = min(1,max(0, (exp(-tr.o2(k) / 10))));
        z_umax = params.(opt.auto_ind{i}).z_umax_0 * Tfunc * (1 - z_o2lim);
        if work1 > 0
            % Get total grazing and route it based on prey fraction
            auto_graze(i) = (Pprime(i) / work1) * ...
				z_umax * tr.zoo(k) * (work1 / (work1 + params.(opt.auto_ind{i}).z_grz));    
        else
            auto_graze(i) = 0;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % N fixation (diazotrophs only) 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %-----------------------------------------------------------------------
        %  Get N fixation by diazotrophs based on C fixation,
        %  Diazotrophs fix more than they need, then 20% is excreted
        %-----------------------------------------------------------------------
        if params.(opt.auto_ind{i}).Nfixer
            work1 = photoC(i) * params.bec.Q;
            Nfix(i) = (work1 * params.bec.r_Nfix_photo) - ...
                NO3_V(i) - NO2_V(i) - NH4_V(i);
            Nexcrete(i) = Nfix(i) + NO3_V(i) + NO2_V(i) + NH4_V(i) - work1;
        else
            Nfix(i)     = 0;
            Nexcrete(i) = 0;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Route grazing and loss to Zoo, POC, DIC, and DOC 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %---------------------------------------------------------%
        % all aggregation goes to POC
        % currently assumes that 33% of grazed caco3 is remineralized
        % if autotrophs(sp_ind)%graze_zoo ever changes, coefficients on routing grazed sp must change!
        % min.%C routed to POC from grazing for ballast requirements = 0.4 * QCaCO3
        % min.%C routed from sp_loss = 0.59 * QCaCO3, or P_CaCO3%rho
        % if autotrophs(diat_ind)%graze_zoo is changed, coeff.s for poc,doc and dic must change!
        %-----------------------------------------------------------------------
        % Route grazing to zooplankton biomass
        auto_graze_zoo(i) = params.(opt.auto_ind{i}).graze_zoo * auto_graze(i);

        % Route grazing to POC
        if params.(opt.auto_ind{i}).imp_calcifier
            % Include CaCO3 production from calcifiers
            auto_graze_poc(i) = auto_graze(i) * ...
                max((params.bec.caco3_poc_min * QCaCO3(i)),...
                min(params.bec.spc_poc_fac * ...
                max(1.0, Pprime(i)),params.bec.f_graze_sp_poc_lim));
        else
            auto_graze_poc(i) = params.(opt.auto_ind{i}).graze_poc * auto_graze(i);    
        end

        % Route grazing to DOC
        auto_graze_doc(i) = params.(opt.auto_ind{i}).graze_doc * auto_graze(i);

        % Route razing to DIC (remaining fraction not routed to zoo, poc, doc)
        auto_graze_dic(i) = auto_graze(i) - ...
            (auto_graze_zoo(i) + auto_graze_poc(i) + auto_graze_doc(i)); 

        % Route mortality to POC
        if params.(opt.auto_ind{i}).imp_calcifier
            % For small phyto, base fraction on local CaCO3:C ratio
            auto_loss_poc(i) = QCaCO3(i) * auto_loss(i); 
        else
            % Otherwise, set it to 'loss_poc' fraction
            auto_loss_poc(i) = params.(opt.auto_ind{i}).loss_poc * auto_loss(i);
        end

        % Route mortality to DOC
        auto_loss_doc(i) = (1 - params.bec.labile_ratio) * (auto_loss(i) - auto_loss_poc(i)); 
        % Route mortality to DIC
        auto_loss_dic(i) = params.bec.labile_ratio * (auto_loss(i) - auto_loss_poc(i)); 
    end % end autotroph loop

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Grazing of chemoautotrophs via zoo 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %-----------------------------------------------------------
    % Similar to autotroph grazing, but get fraction out of loop
    %-----------------------------------------------------------
    if opt.EXPLICIT_MICROBES
		if opt.CHEMO_GRAZING
			% Retrieve 'consumable' biomass of zoo prey (autotrophs + chemoautotrophs)
			work1 = sum(Pprime);
			for i = 1:length(opt.chemo_ind)
				work1 = work1 + max(tr.(opt.chemo_ind{i})(k) - params.(opt.chemo_ind{i}).bmin,0);
			end

			% Get grazing based on fraction of total prey
			if opt.CHEMO_GRAZING_O2
				z_o2lim = min(1,max(0, (exp(-tr.o2(k)/10))));
			else
				z_o2lim = 0;
			end
			z_umax = 0.5 * params.sp.z_umax_0 * Tfunc * (1 - z_o2lim);
			for i = 1:length(opt.chemo_ind)
				if work1 > 0
					% Get 'consumable' biomass of chemo
					this_prime = max(tr.(opt.chemo_ind{i})(k) - params.(opt.chemo_ind{i}).bmin,0);
					% Calculate grazing based on local fraction of prey
					% ... and assuming grazing rate is 50% of small phytoplankton
					chemo_graze(i) = (this_prime / work1) * ...
						(z_umax * tr.zoo(k)) * (work1 / (work1 + params.(opt.chemo_ind{i}).kzoo));
				else
					chemo_graze(i) = 0;
				end
			end
		else
			for i = 1:length(opt.chemo_ind)
				chemo_graze(i) = 0;
			end
		end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Grazing of heterotrophs via szoo 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %-----------------------------------------------------------
    % Similar to other grazing, but with small zoo 
    % and no O2 limitation
    %-----------------------------------------------------------
    if opt.EXPLICIT_MICROBES
		if opt.HETERO_GRAZING
			% Retrieve 'consumable' biomass of heterotrophs
			work1 = 0;
			for i = 1:length(opt.hetero_ind)
				work1 = work1 + max(tr.(opt.hetero_ind{i})(k) - params.(opt.hetero_ind{i}).bmin,0);
			end
			% Get maximum grazing rate (modified by temperature)
			if opt.HETERO_GRAZING_O2
				z_o2lim = min(1,max(0,(exp(-tr.o2(k) / 10))));
			else
				z_o2lim = 0;  
			end
			szoo_umax = params.szoo.u_max * Tfunc * (1 - z_o2lim);
			for i = 1:length(opt.hetero_ind)
				if work1 > 0
					% Get 'consumable' biomass of hetero
					this_prime = max(tr.(opt.hetero_ind{i})(k) - params.(opt.hetero_ind{i}).bmin,0);
					% Get grazing based on local fraction of total prey
					hetero_graze(i) = (this_prime / work1) * ...
						(szoo_umax * tr.szoo(k)) * (work1 / (work1 + params.(opt.hetero_ind{i}).kszoo)); 
				else
					hetero_graze(i) = 0;
				end
			end
		else
			for i = 1:length(opt.hetero_ind)
				hetero_graze(i) = 0;
			end
		end
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Route bacterial grazing to zoo, poc, dic, doc, don, dop
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %-----------------------------------------------------------
    % Assumes same fractioning as in small phytoplankton
    % Also collects 'spill over' effect of consuming nutrient rich bacteria    
    % Because zooplankton and microbes can have different elemental stoichiometries,
    % we need to account for the spill over effect of zooplankton consuming N-rich,
    % P-rich or Fe-rich foods. If consuming rich foods (high quality), then the
    % zooplankton will excrete the excess N, P and Fe --> DOM. If consuming poor quality
    % food (C:N of prey > C:N of Zoo), then these numbers will be negative, and the
    % zooplankton will need to consume extra prey.
    % However, supplementary feeding is yet to be coded, so:
    % C:N, C:P and C:Fe of Microbes <= C:N, C:P and C:Fe of Zoo
    %-----------------------------------------------------------
    if opt.EXPLICIT_MICROBES
        % Chemoautotrophs grazing to zoo, poc, doc, dic
        for i = 1:length(opt.chemo_ind)
            chemo_graze_zoo(i) = params.sp.graze_zoo * chemo_graze(i); % grazing --> zoo (biomass)
            chemo_graze_poc(i) = params.sp.graze_poc * chemo_graze(i); % grazing --> POC (sloppy feeding)
            chemo_graze_doc(i) = params.sp.graze_doc * chemo_graze(i); % grazing --> DOC (excretion of DOC)
            chemo_graze_dic(i) = chemo_graze(i) - ...                   % grazing --> DIC (excretion of DIC)
                (chemo_graze_zoo(i) + chemo_graze_poc(i) + chemo_graze_doc(i)); 
            chemo_graze_don(i) = (chemo_graze_doc(i) / params.(opt.chemo_ind{i}).CN); % via C:N of bacteria
            chemo_graze_dop(i) = (chemo_graze_doc(i) / params.(opt.chemo_ind{i}).CP); % via C:P of bacteria
            chemo_graze_din(i) = ... % NH4 production due to N-rich bacteria (> C:N of zoo)
                (chemo_graze_zoo(i) / params.(opt.chemo_ind{i}).CN) - ...
                (chemo_graze_zoo(i) * params.bec.Q) + ...
                (chemo_graze_poc(i) / params.(opt.chemo_ind{i}).CN) - ...
                (chemo_graze_poc(i) * params.bec.Q); 
            chemo_graze_dip(i) = ... % PO4 production due to P-rich bacteria (> C:P of zoo)
                (chemo_graze_zoo(i) / params.(opt.chemo_ind{i}).CP) - ...
                (chemo_graze_zoo(i) * params.bec.Qp_zoo_pom) + ...
                (chemo_graze_poc(i) / params.(opt.chemo_ind{i}).CP) - ...
                (chemo_graze_poc(i) * params.bec.Qp_zoo_pom); 
            chemo_graze_fe(i) = ... % dissolved Fe production due to Fe-rich bacteria (> C:Fe of zoo)
                (chemo_graze_zoo(i) / params.(opt.chemo_ind{i}).CFe) - ...
                (chemo_graze_zoo(i) * params.bec.Qfe_zoo);
        end
        % Heterotrophs grazing to szoo, poc, doc, dic
		for i = 1:length(opt.hetero_ind);
			hetero_graze_szoo(i) = params.sp.graze_zoo * hetero_graze(i); % grazing --> zoo (biomass)
			hetero_graze_poc(i)  = params.sp.graze_poc * hetero_graze(i); % grazing --> POC (sloppy feeding)
			hetero_graze_doc(i)  = params.sp.graze_doc * hetero_graze(i); % grazing --> DOC (excretion of DOC)
			hetero_graze_dic(i)  = hetero_graze(i) - ...                   % grazing --> DIC (excretion of DIC)
				(hetero_graze_szoo(i) + hetero_graze_poc(i) + hetero_graze_doc(i)); 
			hetero_graze_don(i) = (hetero_graze_doc(i) / params.(opt.hetero_ind{i}).CN); % via C:N of bacteria
			hetero_graze_dop(i) = (hetero_graze_doc(i) / params.(opt.hetero_ind{i}).CP); % via C:N of bacteria
			hetero_graze_din(i) = ... % NH4 production due to N-rich bacteria (> C:N of zoo)
				(hetero_graze_szoo(i) / params.(opt.hetero_ind{i}).CN) - ...
				(hetero_graze_szoo(i) * params.bec.Q) + ...
				(hetero_graze_poc(i) / params.(opt.hetero_ind{i}).CN) - ...
				(hetero_graze_poc(i) * params.bec.Q);  
			hetero_graze_dip(i) = ... % PO4 production due to P-rich bacteria (> C:P of zoo)
				(hetero_graze_szoo(i) / params.(opt.hetero_ind{i}).CP) - ...
				(hetero_graze_szoo(i) * params.bec.Qp_zoo_pom) + ...
				(hetero_graze_poc(i) / params.(opt.hetero_ind{i}).CP) - ...
				(hetero_graze_poc(i) * params.bec.Qp_zoo_pom); 
			hetero_graze_fe(i) = ... % dissolved Fe production due to Fe-rich bacteria (> C:Fe of zoo)
				(hetero_graze_szoo(i) / params.(opt.hetero_ind{i}).CFe) - ...
				(hetero_graze_szoo(i) * params.bec.Qfe_zoo);
		end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Losses of zooplankton 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %-----------------------------------------------------------------------% 
    % Losses of zooplankton to detritus (i.e., POC) are scaled by what they eat. 
    % If they eat more diatoms, then more of their losses are routed to POC
    % 
    % Then, calculate losses of zooplankton via quadratic closure, and 
    % route these losses to DOC, DIC, and POC
    %--------------------------------------------------------------
    % Get fractional factor for routing of zoo losses, based on food supply
    % More material is routed to large detrital pool when diatoms eaten
    work1 = 0;
    work2 = 0;
	% Get total fractions of autotroph grazing routed to POC 
    for i = 1:length(opt.auto_ind);
        work1 = work1 + params.(opt.auto_ind{i}).f_zoo_detr * ...
            (auto_graze(i) + (params.bec.epsC * params.bec.epsTinv));
        work2 = work2 + ...
            (auto_graze(i) + (params.bec.epsC * params.bec.epsTinv));
    end
	% Get total fractions of chemoautotroph grazing routed to POC 
	% Fractions are based off small phytoplankton type
	for i = 1:length(opt.chemo_ind)
		work1 = work1 + params.sp.f_zoo_detr * ...
			(chemo_graze(i) + (params.bec.epsC * params.bec.epsTinv));
        work2 = work2 + ...
            (chemo_graze(i) + (params.bec.epsC * params.bec.epsTinv));
	end
    f_zoo_detr_loc = work1 / work2;

    % Get concentration beneath which losses do not occur 
    C_loss_thres = f_loss_thres * params.zoo.bmin;  
    Zprime = max(tr.zoo(k) - C_loss_thres,0);

    % Calculate loss of zooplankton (both quadratic and linear mortality)
    zoo_loss = ((params.zoo.m_q * Zprime^1.5) + (params.zoo.m_l * Zprime)) * Tfunc;

    % Route zooplankton losses to DOC, DIC and POC
    zoo_loss_doc = zoo_loss * (1 - params.bec.labile_ratio) * (1 - f_zoo_detr_loc);
    zoo_loss_dic = zoo_loss * (1 - f_zoo_detr_loc) * params.bec.labile_ratio ;
    zoo_loss_poc = zoo_loss * f_zoo_detr_loc;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Losses of small zooplankton
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %-----------------------------------------------------------------------% 
    % Same as above for zooplankton
    %-----------------------------------------------------------------------% 
    if opt.EXPLICIT_MICROBES

		% Get fractional factor for routing of zoo losses, based on food supply
		work1 = 0;
		work2 = 0;
		% Get total fractions of heterootroph grazing routed to POC
		% Fractions are based off small phytoplankton type
		for i = 1:length(opt.hetero_ind)
			work1 = work1 + params.sp.f_zoo_detr * ...
				(hetero_graze(i) + (params.bec.epsC * params.bec.epsTinv));
			work2 = work2 + ...
				(hetero_graze(i) + (params.bec.epsC * params.bec.epsTinv));
		end
		f_szoo_detr_loc = work1 / work2;

		% Calculate the concentration beneath which losses do not occur
		C_loss_thres = f_loss_thres * params.szoo.bmin;
		sZprime = max(tr.szoo(k) - C_loss_thres, 0);

		% Calculate loss of small zooplankton (both quadratic and linear mortality)
		szoo_loss = ((params.szoo.m_q * sZprime^1.5) + (params.szoo.m_l * sZprime)) * Tfunc;

		% Route zooplankton losses to DOC, DIC and POC
		szoo_loss_doc = szoo_loss * (1 - params.bec.labile_ratio) * (1 - f_szoo_detr_loc);
		szoo_loss_dic = szoo_loss * (0 + params.bec.labile_ratio) * (1 - f_szoo_detr_loc);
		szoo_loss_poc = szoo_loss * f_szoo_detr_loc;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mortality of bacteria 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %-----------------------------------------------------------------------% 
    % Calculate total mortality (linear + quadratic)
    % Then, route to doc, dic, poc, don, dop
    %-----------------------------------------------------------------------% 
    if opt.EXPLICIT_MICROBES
        % Chemoautotrophs
        for i = 1:length(opt.chemo_ind);
            % Calculate the concentration beneath which losses do not occur 
            this_prime = max(tr.(opt.chemo_ind{i})(k) - params.(opt.chemo_ind{i}).bmin,0);
            
            % Calculate losses (both quadratic and linear)
            chemo_loss(i) = (...
                (params.(opt.chemo_ind{i}).m_l * this_prime) + ...
                (params.(opt.chemo_ind{i}).m_q * (this_prime * this_prime))) * Tfunc;

            % Route losses to DOC, POC, DIC
            chemo_loss_doc(i) = chemo_loss(i) * ... % Mortality --> DOC
                params.(opt.chemo_ind{i}).mortdoc;
            chemo_loss_poc(i) = chemo_loss(i) * ... % Mortality --> POC
                params.(opt.chemo_ind{i}).mortpoc;
            chemo_loss_dic(i) = chemo_loss(i) * ... % Mortality --> DIC
                params.(opt.chemo_ind{i}).mortdic;

            % Route losses to DON, DOP
            % 1. Get total C production (chemo_loss), convert to N (or P) via C:N (or C:P) of bacteria
            % 2. Subtract fraction of N (or P) production routed to DIN (or DIP)
            % 3. Subtract fraction of N (or P) production routed to PON (or POP)
            % 4. Remainder is N (or P) routed to DON (or DOP)
			% NOTE: This is different to the partitioning via auto losses to DON/DOP
			% This is because P and Z have the same stoichiometry, so you can simply
			% Convert all DOC production to DON/DOP via BEC ratios
            chemo_loss_don(i) = ...
                (chemo_loss(i)     / params.(opt.chemo_ind{i}).CN) - (... % total N from biomass mortality 
                (chemo_loss_dic(i) / params.(opt.chemo_ind{i}).CN) +  ... % subtract DIN production 
                (chemo_loss_poc(i) * params.bec.Q));                      % subtract PON production 
            chemo_loss_dop(i) = ...
                (chemo_loss(i)     / params.(opt.chemo_ind{i}).CP) - (... % total P from biomass mortality
                (chemo_loss_dic(i) / params.(opt.chemo_ind{i}).CP) +  ... % subtract DIP production 
                (chemo_loss_poc(i) * params.bec.Qp_zoo_pom));             % subtract POP production         
        end

        % Heterotrophs
        for i = 1:length(opt.hetero_ind);
            % Calculate the concentration beneath which losses do not occur 
            this_prime = max(tr.(opt.hetero_ind{i})(k) - params.(opt.hetero_ind{i}).bmin,0);
            
            % Calculate losses (both quadratic and linear)
            hetero_loss(i) = (...
                (params.(opt.hetero_ind{i}).m_l * this_prime) + ...
                (params.(opt.hetero_ind{i}).m_q * (this_prime * this_prime))) * Tfunc;

            % Route losses to DOC, POC, DIC
            hetero_loss_doc(i) = hetero_loss(i) * ... % Mortality --> DOC
                (params.(opt.hetero_ind{i}).mortdoc);
            hetero_loss_poc(i) = hetero_loss(i) * ... % Mortality --> POC
                params.(opt.hetero_ind{i}).mortpoc;
            hetero_loss_dic(i) = hetero_loss(i) * ... % Mortality --> DIC
                params.(opt.hetero_ind{i}).mortdic;

            % Route losses to DON, DOP
            % 1. Get total C production (chemo_loss), convert to N (or P) via C:N (or C:P) of bacteria
            % 2. Subtract fraction of N (or P) production routed to DIN (or DIP)
            % 3. Subtract fraction of N (or P) production routed to PON (or POP)
            % 4. Remainder is N (or P) routed to DON (or DOP)
			% NOTE: This is different to the partitioning via auto losses to DON/DOP
			% This is because P and Z have the same stoichiometry, so you can simply
			% Convert all DOC production to DON/DOP via BEC ratios
            hetero_loss_don(i) = ...
                (hetero_loss(i)     / params.(opt.hetero_ind{i}).CN) - ( ... % total N from biomass mortality
                (hetero_loss_dic(i) / params.(opt.hetero_ind{i}).CN) +   ... % subtract DIN production 
                (hetero_loss_poc(i) * params.bec.Q));                        % subtract PON production 
            hetero_loss_dop(i) = ...
                (hetero_loss(i)     / params.(opt.hetero_ind{i}).CP) - ( ... % total P from biomass mortality
                (hetero_loss_dic(i) / params.(opt.hetero_ind{i}).CP) +   ... % subtract DIP production 
                (hetero_loss_poc(i) * params.bec.Qp_zoo_pom));               % subtract POP production 
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Collect POM production 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % POC prod
    poc_prod = ...
		+ zoo_loss_poc ...        % zoo mortality routed to POC
		+ sum(auto_graze_poc) ... % auto grazing routed to POC
		+ sum(auto_agg) ...       % auto aggregation routed to POC
		+ sum(auto_loss_poc);     % auto mortality routed to POC
    if opt.EXPLICIT_MICROBES
        % Include bacterial production
        poc_prod = poc_prod ...
			+ sum(chemo_loss_poc) ...   % chemo mortality routed to POC
			+ sum(chemo_graze_poc) ...  % chemo grazing routed to POC
			+ sum(hetero_graze_poc) ... % hetero grazing routed to POC
			+ sum(hetero_loss_poc) ...  % hetero mortality routed to POC
			+ szoo_loss_poc;            % small zoo mortality routed to POC
    end

    % CaCO3 prod    
    % In BEC, 33% (f_graze_CaCO3_remin) of CaCO3 is remineralized when phyto are grazed
    % Anything left over is routed to particulate pool
    P_CaCO3_prod = ((1 - params.bec.f_graze_CaCO3_remin) * auto_graze(strcmp('sp',opt.auto_ind)) + ...
        auto_loss(strcmp('sp',opt.auto_ind)) + ...
        auto_agg(strcmp('sp',opt.auto_ind))) * QCaCO3(strcmp('sp',opt.auto_ind)); 

    % SiO2 prod
    % In BEC, 35% (f_graze_si_remin) of SiO2 is remineralized when diatoms are grazed
    % Anything left over is routed to particulate pool
    P_SiO2_prod = 0;
    for i = 1:length(opt.auto_ind)
        if params.(opt.auto_ind{i}).K_sio3 > 0
            P_SiO2_prod = P_SiO2_prod + ...
                Qsi(i) * ((1 - params.bec.f_graze_si_remin) * auto_graze(i) + ...
                auto_agg(i) + params.(opt.auto_ind{i}).loss_poc * auto_loss(i));    
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Iron scavenging and particulate iron production 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %-----------------------------------------------------------------------
    % Compute in terms of loss per year per unit iron (%/year/fe)
    % Scale by sinking POMx10 + Dust + bSi + CaCO3 flux
    %-----------------------------------------------------------------------
    % Get initial scavenging rate
    Fe_scavenge_rate = params.bec.Fe_scavenge_rate0;

    % Modify based on sinking pools
    Fe_scavenge_rate = Fe_scavenge_rate * (...
        (POC_sflux_out + POC_hflux_out) * params.bec.POC_mass * 10 + ...
        (P_CaCO3_sflux_out + P_CaCO3_hflux_out) * params.bec.CaCO3_mass + ...
        (P_SiO2_sflux_out + P_SiO2_hflux_out) * params.bec.SiO2_mass + ...
        (dust_sflux_out + dust_hflux_out) * params.bec.dust_fescav_scale); 

    % Increase scavenging at higher iron
    if tr.fe(k) > params.bec.fe_scavenge_thres1 
        Fe_scavenge_rate = Fe_scavenge_rate + ...
            (tr.fe(k) - params.bec.fe_scavenge_thres1) * params.bec.fe_max_scale2; 
    end

    % Convert to loss-per-second
    Fe_scavenge = tr.fe(k) * Fe_scavenge_rate * params.bec.yps; 

    % Now compute total particulate iron prod (including iron scavenging by particles)
    P_iron_prod = Fe_scavenge;
    if opt.EXPLICIT_MICROBES
        % ... add chemo losses to Fe production
        for i = 1:length(opt.chemo_ind)
            P_iron_prod = P_iron_prod + ...
                (chemo_loss_poc(i) / params.(opt.chemo_ind{i}).CFe) + ...
                (chemo_graze_poc(i) / params.(opt.chemo_ind{i}).CFe); 
        end
		% ... add hetero losses to Fe production
        for i = 1:length(opt.hetero_ind);
            P_iron_prod = P_iron_prod + ...
                (hetero_loss_poc(i) / params.(opt.hetero_ind{i}).CFe) + ...
				(hetero_graze_poc(i) / params.(opt.hetero_ind{i}).CFe);
		end
        % ... add zoo losses to Fe production 
		P_iron_prod = P_iron_prod + ...
			((zoo_loss_poc + szoo_loss_poc) * params.bec.Qfe_zoo);
    end
    % ... add phytoplankton mortality, grazing, aggregation to Fe production
    for i = 1:length(opt.auto_ind)
        P_iron_prod = P_iron_prod + ...
            Qfe(i) * (auto_agg(i) + auto_graze_poc(i) + auto_loss_poc(i));            
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sinking of POC and its remineralization 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    compute_particulate_terms

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Collect POM production 
	% NOTE: No longer tracking dofe, all routed to fe
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DOC
    doc_prod = ...
		+ zoo_loss_doc ...        % zoo mortality routed directly to DOC
		+ sum(auto_loss_doc) ...  % auto mortality routed directly to DOC
		+ sum(auto_graze_doc) ... % auto grazing routed directly to DOC
        + poc_remin;              % z-level 'remineralization' of POC to DOC
    if opt.EXPLICIT_MICROBES
        % ... include bacterial production
        doc_prod = doc_prod ...
			+ sum(chemo_loss_doc) ...   % chemo mortality routed directly to DOC
		    + sum(chemo_graze_doc) ...  % chemo grazing routed directly to DOC
		    + sum(hetero_loss_doc) ...  % hetero mortality routed directly to DOC
			+ sum(hetero_graze_doc) ... % hetero grazing routed directly to DOC
			+ szoo_loss_doc;            % small zoo mortality routed directly to DOC
    end
    
    % DON    
    don_prod = ...
		+ (params.bec.Q * zoo_loss_doc) ...        % same as for DOC, but converted to DON
		+ (params.bec.Q * sum(auto_loss_doc)) ...  % same as for DOC, but converted to DON
        + (params.bec.Q * sum(auto_graze_doc)) ... % same as for DOC, but converted to DON
	    + (params.bec.Q * poc_remin);              % same as for DOC, but converted to DON
    if opt.EXPLICIT_MICROBES
        % ... include bacterial and small zoo production
        don_prod = don_prod ... 
			+ sum(chemo_loss_don) ...          % DOC converted to DON (considers different C:N)
			+ sum(chemo_graze_don) ...         % DOC converted to DON (considers different C:N)
			+ sum(hetero_loss_don) ...         % DOC converted to DON (considers different C:N)
			+ sum(hetero_graze_don) ...        % DOC converted to DON (considers different C:N)
			+ (params.bec.Q * szoo_loss_doc);  % same as for DOC, but converted to DON 
    end

    % DOP
    dop_prod = ...
		+ (params.bec.Qp_zoo_pom * zoo_loss_doc) ...        % same as for DOC, but converted to DOP
		+ (params.bec.Qp_zoo_pom * sum(auto_loss_doc)) ...  % same as for DOC, but converted to DON
        + (params.bec.Qp_zoo_pom * sum(auto_graze_doc)) ... % same as for DOC, but converted to DON
		+ (params.bec.Qp_zoo_pom * poc_remin);              % same as for DOC, but converted to DOP
    if opt.EXPLICIT_MICROBES
        % ... include bacterial and small zoo production
        dop_prod = dop_prod ...
            + sum(chemo_loss_dop) ...                   % DOC converted to DOP (considers different C:N)
			+ sum(chemo_graze_dop) ...                  % DOC converted to DOP (considers different C:N)
			+ sum(hetero_loss_dop) ...                  % DOC converted to DOP (considers different C:N)
			+ sum(hetero_graze_dop) ...                 % DOC converted to DOP (considers different C:N)
			+ (params.bec.Qp_zoo_pom * szoo_loss_doc);  % same as for DOC, but converted to DOP 
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Chemoautotroph metabolism (anammox + nitrif) 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Chemoautotroph stoichiometry is simple
    % They require a certain concentration of oxidants and reductants
    % to build biomass in a fixed ratio (C_5 H_7 O_2.6 N_1)
    % See (functions/get_stoichiometry.m)
    if opt.EXPLICIT_MICROBES
		if tr.o2(k) > 0
			if opt.KELLY_N2O_YIELD
				[yaoa_n2o_nh4,yaoa_n2o_hyb,yaoa_no2] = n2o_yield_Kelly(tr.o2(k),params);
            else
				[yaoa_n2o_nh4,yaoa_n2o_hyb,yaoa_no2] = n2o_yield_Ji(tr.o2(k),params);
            end
        else
            yaoa_n2o_nh4 = 0;
			yaoa_n2o_hyb = 0;
            yaoa_no2 = 0;
        end

        % Calculate growth rate (1 / s)
        for i = 1:length(opt.chemo_ind);
            this_oxy = tr.(params.(opt.chemo_ind{i}).oxy)(k); % uM oxidant 
            this_red = tr.(params.(opt.chemo_ind{i}).red)(k); % uM reductant
            % Compute oxidant and reductant uptake
            if this_oxy > 0 & this_red > 0
                if strcmp(params.(opt.chemo_ind{i}).oxyup,'diffusive');
                    % Diffusive rate (m3/mmol BC s) * (mmol oxidant / m3) =(mol oxy / mol BC s)
                    oxy_up = params.(opt.chemo_ind{i}).pcoef * this_oxy; 
                elseif strcmp(params.(opt.chemo_ind{i}).oxyup,'michaelis');
                    % VMax (mol oxy / mol BC s) * (%) = (mol oxy / mol BC s)
                    oxy_up = params.(opt.chemo_ind{i}).Vmax_oxy * ...
                        (this_oxy ./ (this_oxy + params.(opt.chemo_ind{i}).K_oxy));
                end
                % Vmax (mol red / mol BC s) * (%) = (mol red / mol BC s)
                red_up = params.(opt.chemo_ind{i}).Vmax_red * ...
                    (this_red ./ (this_red + params.(opt.chemo_ind{i}).K_red));
            else
                oxy_up = 0;
                red_up = 0;
            end
            % Get Liebig's 'Law of the minimum' 
            % (mol red / mol BC s) * (mol BC / mol red) = (1 / s)
            % (mol oxy / mol BC s) * (mol BC / mol oxy) = (1 / s)
            chemo_mu = max(0, (min(...
                (oxy_up / params.(opt.chemo_ind{i}).y_oxy),...
                (red_up / params.(opt.chemo_ind{i}).y_red))))*Tfunc;
            % Get growth rate
            % (1 / s) * (mmol BC / m3) = (mmol BC / m3 s)
            chemo_growth(i) = chemo_mu * tr.(opt.chemo_ind{i})(k); 
        end
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Heterotroph metabolism (aerobic + denitrif remin) 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Heterotroph stoichiometry is more complicated
    % Here, when estimating the yields (functions/get_stoichiometry), 
    % we assume  the stoichiometry of DOM is similar to 
    % Anderson & Sarmiento, 1995: [C_117 H_297 O_85 N_16] 
    % NOTE: You can change this in options.m
    % However, the ratios of local DOM can differ from these ratios
    % Later on, we'll calculate the local ratios, and any excess
    % or lack of N/P will be taken up from the inorganic pool  

    if opt.EXPLICIT_MICROBES
        doc_tot = tr.doc(k) + tr.docr(k);
        for i = 1:length(opt.hetero_ind)
            % Get aerobic and anaerobic oxidants
            this_oxy_aer = tr.(params.(opt.hetero_ind{i}).oxy_aer)(k);
            this_oxy_ana = tr.(params.(opt.hetero_ind{i}).oxy_ana)(k);
            % Get limitation for DOC and DOCr uptake (reducants)
            doc_lim(i)  = tr.doc(k)  / (tr.doc(k)  + params.(opt.hetero_ind{i}).K_doc); 
            docr_lim(i) = tr.docr(k) / (tr.docr(k) + params.(opt.hetero_ind{i}).K_docr); 
            % Get both aerobic and anaerobic metabolisms
            % Aerobic metabolism
            if this_oxy_aer > 0 & doc_tot > 0 
                if strcmp(params.(opt.hetero_ind{i}).oxyup_aer,'diffusive');
                    % Diffusive rate (m3/mmol BC s) * (mmol oxidant / m3) = mmol oxy / mmol BC s
                    oxy_up_aer = params.(opt.hetero_ind{i}).pcoef * this_oxy_aer; 
                elseif strcmp(params.(opt.hetero_ind{i}).oxyup_aer,'michaelis');
                    % Vmax (uM oxy / uM BC s) * (%) = uM oxy / BC s
                    oxy_up_aer = params.(opt.hetero_ind{i}).Vmax_oxy_aer * ...
                        (this_oxy_aer / (this_oxy_aer + params.(opt.hetero_ind{i}).K_oxy_aer));
                end
                % Vmax (uM red / BC s) * (% + %) = uM red / BC s
                red_up_aer = params.(opt.hetero_ind{i}).Vmax_red_aer * (doc_lim(i) + docr_lim(i));
            else
                oxy_up_aer = 0;
                red_up_aer = 0;
            end
            % The equation below solves for Biomass growth
            % (uM oxy / uM BC s) * (uM B / uM oxy) = (1 / s)
            % (uM red / uM BC s) * (uM B / uM red) = (1 / s)
            % The 'loser' determines the substrate (oxidant or reductant) limiting growth rate
            mu_aerobic = max(0, (min(...
                (oxy_up_aer / params.(opt.hetero_ind{i}).y_oxy_aer), ...
                (red_up_aer / params.(opt.hetero_ind{i}).y_red_aer))) * Tfunc);

            % Anaerobic metabolism
            if this_oxy_ana > 0 & doc_tot > 0
                if strcmp(params.(opt.hetero_ind{i}).oxyup_ana,'diffusive');
                    % Diffusive rate (m3/mmol BC s) * (mmol oxidant / m3) = mmol oxy / mmol BC s 
                    oxy_up_ana = params.(opt.hetero_ind{i}).pana * this_oxy_ana;
                elseif strcmp(params.(opt.hetero_ind{i}).oxyup_ana,'michaelis');
                    % Vmax (uM oxy / uM BC s) * (%) = uM oxy / BC s 
                    oxy_up_ana = params.(opt.hetero_ind{i}).Vmax_oxy_ana * ...
                        (this_oxy_ana / (this_oxy_ana + params.(opt.hetero_ind{i}).K_oxy_ana));
                end
                % Vmax (mol red / mol BC s) * (%) = mol oxy / mol BC s 
                red_up_ana = params.(opt.hetero_ind{i}).Vmax_red_ana * (doc_lim(i) + docr_lim(i));
            else
                oxy_up_ana = 0;
                red_up_ana = 0;
            end
            % (uM red / s) * (1 / uM red) = (1 / s)
            % (uM oxy / s) * (1 / uM oxy) = (1 / s)
            mu_anaerobic = max(0, (min(...
                (oxy_up_ana / params.(opt.hetero_ind{i}).y_oxy_ana), ...
                (red_up_ana / params.(opt.hetero_ind{i}).y_red_ana))) * Tfunc);

            % Get the winning metabolism
            if opt.FACULTATIVE_MICROBES
                if mu_anaerobic > mu_aerobic
                    facultative(i) = 1; % anaerobic wins
                else
                    facultative(i) = 0; % aerobic wins
                end 
                hetero_mu = max(mu_anaerobic,mu_aerobic);
                if i == 1 % aer
                    hetero_mu      = mu_aerobic; % override for obligate aerobe 'aer'
                    facultative(i) = 0;
                end
            else
                if i == 1 % aer
                    facultative(i) = 0;
                    hetero_mu      = mu_aerobic; % override for obligate aerobe 'aer'
                else
                    facultative(i) = 1; % force anaerobic
                    hetero_mu = mu_anaerobic;
                end
            end
            % (1 / s) * (uM BC) = uM BC / s
            hetero_growth(i) = max(0, (hetero_mu * tr.(opt.hetero_ind{i})(k)));
        end    
    end
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute remineralization of DOM into C,N,P pools 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %-----------------------------------------------------------------------
    % We treat the remineralisation of C,N,P differently from Fe.
    % N and P are remineralised in their ratio with C, since heterotrophs get N and P from DOM.
    % Fe, however, is assimilated as dFe by bacteria, not consumed as part of DOM. No DOFe pool.
    % For N and P, first get the elemental ratios of DOM in the environment.
    %-----------------------------------------------------------------------
    % Get the elemental ratios of DOM in the environment
    % DOM (C N_x P_y Fe_z)
    qDON  = (tr.don(k) + params.bec.epsN)  / (tr.doc(k) + params.bec.epsN);  % N_x
    qDOP  = (tr.dop(k) + params.bec.epsN)  / (tr.doc(k) + params.bec.epsN);  % P_y
    qDONr = (tr.donr(k) + params.bec.epsN) / (tr.docr(k) + params.bec.epsN); % Nr_x
    qDOPr = (tr.dopr(k) + params.bec.epsN) / (tr.docr(k) + params.bec.epsN); % Pr_y

    % Calculate DOC consumption based on growth and yields
    % Total growth is partitioned between DOC and DOCr remineralization
    doc_remin  = 0;
    docr_remin = 0;
    for i = 1:length(opt.hetero_ind);
        doc_remin = doc_remin + hetero_growth(i) * (...
            (params.(opt.hetero_ind{i}).y_red_aer * (1 - facultative(i))) + ...    % aerobic or
            (params.(opt.hetero_ind{i}).y_red_ana * (0 + facultative(i)))) * (...  % anaerobic
            (doc_lim(i) + params.bec.epsN) / (doc_lim(i) + docr_lim(i) + params.bec.epsN)); 
        docr_remin = docr_remin + hetero_growth(i) * (...
            (params.(opt.hetero_ind{i}).y_red_aer * (1 - facultative(i))) + ...    % aerobic or
            (params.(opt.hetero_ind{i}).y_red_ana * (0 + facultative(i)))) * (...  % anaerobic
            (docr_lim(i) + params.bec.epsN) / (doc_lim(i) + docr_lim(i) + params.bec.epsN));
    end

    % Apply the elemental ratios
    % This assumes that microbes indiscriminately munch on DOM that is available
    % i.e., that they consume whatever DOM is available in the stoich that is available
    % This approach should also ensure that DON and DOP never become negative
    don_remin  = max(0,(doc_remin * qDON));   % Amount of DON consumed
    dop_remin  = max(0,(doc_remin * qDOP));   % Amount of DOP consumed
    donr_remin = max(0,(docr_remin * qDONr)); % Amount of DONr consumed
    dopr_remin = max(0,(docr_remin * qDOPr)); % Amount of DOPr consumed
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute sources and sinks for BGC tracers 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% ------- %
    % Nitrate %
	% ------- %
    sms.no3(k,1) = 0 ...
        + (chemo_growth(opt.nob_ind) * params.nob.e_no3) ...                                 % nob NO3 excretion
        + (chemo_growth(opt.aox_ind) * params.aox.e_no3) ...                                 % aox NO3 excretion
        - (hetero_growth(opt.nar_ind) * facultative(opt.nar_ind) * params.nar.y_oxy_ana) ... % nar NO3 yield
        - (hetero_growth(opt.nai_ind) * facultative(opt.nai_ind) * params.nai.y_oxy_ana) ... % nai NO3 yield
        - (hetero_growth(opt.nao_ind) * facultative(opt.nao_ind) * params.nao.y_oxy_ana) ... % nao NO3 yield
        - (sed_denitrif) ...                                                                 % NO3 --> N2 from sediments
        - (sum(NO3_V));                                                                      % NO3 --> autotroph biomass
	if opt.DNRA1
		sms.no3(k,1) = sms.no3(k,1) ...
			- (hetero_growth(opt.dnra1_ind) * facultative(opt.dnra1_ind) * params.dnra1.y_oxy_ana); % DNRA1 NO3 yield
	end

	% ------- %
    % Nitrite %
	% ------- %
    sms.no2(k,1) = 0 ...
        + (hetero_growth(opt.nar_ind) * facultative(opt.nar_ind) * params.nar.e_no2_ana) ... % nar NO2 excretion
        + (chemo_growth(opt.aoa_ind) * params.aoa.e_no2 * yaoa_no2) ...                      % aoa NO2 excretion * NO2 fraction 
		- (chemo_growth(opt.aoa_ind) * params.aoa.e_no2 * 0.5 * yaoa_n2o_hyb) ...            % hybrid N2O production (uses NO2)
        - (chemo_growth(opt.nob_ind) * params.nob.y_red) ...                                 % nob NO2 yield
        - (chemo_growth(opt.aox_ind) * params.aox.y_oxy) ...                                 % aox NO2 yield 
        - (hetero_growth(opt.nir_ind) * facultative(opt.nir_ind) * params.nir.y_oxy_ana) ... % nir NO2 yield
        - (hetero_growth(opt.nio_ind) * facultative(opt.nio_ind) * params.nio.y_oxy_ana) ... % nio NO2 yield
        - (sum(NO2_V));                                                                      % NO2 --> autotroph biomass
	if opt.DNRA2
		sms.no2(k,1) = sms.no2(k,1) ...
			- (hetero_growth(opt.dnra2_ind) * facultative(opt.dnra2_ind) * params.dnra2.y_oxy_ana); % DNRA2 NO2 yield
	end

	% ------------- %
    % Nitrous oxide %
	% ------------- %
    sms.n2o(k,1) = 0 ...
        + (chemo_growth(opt.aoa_ind) * params.aoa.e_no2 * 0.5 * yaoa_n2o_hyb) ...            % aoa NO2 exretion * N2O fraction
        + (chemo_growth(opt.aoa_ind) * params.aoa.e_no2 * 0.5 * yaoa_n2o_nh4) ...            % aoa NO2 exretion * N2O fraction
        + (hetero_growth(opt.nai_ind) * facultative(opt.nai_ind) * params.nai.e_n2o_ana) ... % nai N2O excretion
        + (hetero_growth(opt.nir_ind) * facultative(opt.nir_ind) * params.nir.e_n2o_ana) ... % nir N2O excretion 
        - (hetero_growth(opt.nos_ind) * facultative(opt.nos_ind) * params.nos.y_oxy_ana);    % nos N2O yield

	% ---------- %
    % Dinitrogen %
	% ---------- %
    sms.n2(k,1) = 0 ...
        + (chemo_growth(opt.aox_ind) * params.aox.e_n2) ...                                 % aox N2 excretion
        + (hetero_growth(opt.nao_ind) * facultative(opt.nao_ind) * params.nao.e_n2_ana) ... % nao N2 excretion
        + (hetero_growth(opt.nio_ind) * facultative(opt.nio_ind) * params.nio.e_n2_ana) ... % nio N2 excretion
        + (hetero_growth(opt.nos_ind) * facultative(opt.nos_ind) * params.nos.e_n2_ana) ... % nos N2 excretion
        + (0.5 * sed_denitrif);                                                             % NO3 --> N2 from sediments 
            
	% -------- %
    % Ammonium % 
	% -------- %
    sms.nh4(k,1) = 0 ...
        + (params.bec.Q * zoo_loss_dic) ...        % zoo mortality to NH4
		+ (params.bec.Q * szoo_loss_dic) ...       % szoo mortality to NH4
        + (params.bec.Q * sum(auto_loss_dic)) ...  % autotroph mortality to NH4
        + (params.bec.Q * sum(auto_graze_dic)) ... % autotroph grazing to NH4
        - (sum(NH4_V));                            % NH4 --> autotroph biomass 
	% ... add chemoautotroph consumption
	sms.nh4(k,1) = sms.nh4(k,1) ...
        - (chemo_growth(opt.aoa_ind) * params.aoa.y_red) ... % aoa NH4 yield
        - (chemo_growth(opt.nob_ind) * params.nob.y_nh4) ... % nob NH4 yield
        - (chemo_growth(opt.aox_ind) * params.aox.y_red);    % aox NH4 yield   
	% ... add heterotroph excretion
	for i = 1:length(opt.hetero_ind)
		sms.nh4(k,1) = sms.nh4(k,1) ...
			+ (hetero_growth(i) * (0 + facultative(i)) * params.(opt.hetero_ind{i}).e_nh4_ana) ...  % anaerobic NH4 excretion
			+ (hetero_growth(i) * (1 - facultative(i)) * params.(opt.hetero_ind{i}).e_nh4_aer);     % aerobic NH4 excretion
	end
    % ... add chemoautotroph losses, grazing to NH4
    for i = 1:length(opt.chemo_ind);
        sms.nh4(k,1) = sms.nh4(k,1) ...
            + (chemo_loss_dic(i) / params.(opt.chemo_ind{i}).CN) ...  % chemo mortality to NH4 
            + (chemo_graze_dic(i) / params.(opt.chemo_ind{i}).CN) ... % chemo grazing to NH4
            + (chemo_graze_din(i));                                   % chemo grazing to NH4 (N-rich) 
    end
    % ... add heterotroph losses, grazing to NH4
    for i = 1:length(opt.hetero_ind);
        sms.nh4(k,1) = sms.nh4(k,1) ...
            + (hetero_loss_dic(i) / params.(opt.hetero_ind{i}).CN)...   % hetero mortality to NH4 
			+ (hetero_graze_dic(i) / params.(opt.hetero_ind{i}).CN) ... % hetero grazing to NH4
			+ (hetero_graze_din(i));                                    % hetero grazing to NH4 (N-rich) 
    end
    % ... add excretion from autotrophs
    for i = 1:length(opt.auto_ind);
        if params.(opt.auto_ind{i}).Nfixer
            sms.nh4(k,1) = sms.nh4(k,1) + Nexcrete(i); % auto excretion to NH4
        end
    end

	% --------- %
    % Phosphate %
	% --------- %
	sms.po4(k,1) = 0 ...
		+ (params.bec.Qp_zoo_pom * zoo_loss_dic) ...  % zoo mortality to PO4
		+ (params.bec.Qp_zoo_pom * szoo_loss_dic) ... % small zoo mortality to PO4
		- (sum(PO4_V));	                            % PO4 --> autotroph biomass
	% ... include uptake, losses, and excretion via chemos
	for i = 1:length(opt.chemo_ind);
		sms.po4(k,1) = sms.po4(k,1) ...
			- (chemo_growth(i) * params.(opt.chemo_ind{i}).y_po4) ... % chemo uptake of PO4
			+ (chemo_loss_dic(i) / params.(opt.chemo_ind{i}).CP) ...  % chemo mortality to PO4
			+ (chemo_graze_dic(i) / params.(opt.chemo_ind{i}).CP) ... % chemo grazing to PO4
			+ (chemo_graze_dip(i));                                   % chemo grazing to PO4 (P-rich)
	end
	% ... add heterotroph growth, losses
	for i = 1:length(opt.hetero_ind)
		sms.po4(k,1) = sms.po4(k,1) ...
			+ (hetero_loss_dic(i) / params.(opt.hetero_ind{i}).CP) ...                             % hetero mortality to PO4
			+ (hetero_graze_dic(i) / params.(opt.hetero_ind{i}).CP) ...                            % hetero grazing to PO4
			+ (hetero_growth(i) * (0 + facultative(i)) * params.(opt.hetero_ind{i}).e_po4_ana) ... % anaerobic growth to PO4
			+ (hetero_growth(i) * (1 - facultative(i)) * params.(opt.hetero_ind{i}).e_po4_aer) ... % aerobic growth to PO4
			+ (hetero_graze_dip(i));                                                               % hetero grazing to PO4 (P-rich)
	end
	% ... add excretion from autotrophs
	for i = 1:length(opt.auto_ind);
		sms.po4(k,1) = sms.po4(k,1) ...
			+ (params.(opt.auto_ind{i}).Qp * auto_loss_dic(i)) ... % auto mortality to PO4 (DIC)
			+ (params.(opt.auto_ind{i}).Qp * auto_graze_dic(i));   % auto grazing to PO4 (DIC)
	end

	% ---- %
    % Iron % 
	% ---- %
	sms.fe(k,1) = 0 ...
		+ (P_iron_remin) ...                       % particulate remineralization to Fe
		+ (params.bec.Qfe_zoo * zoo_loss_dic) ...  % zoo mortality to Fe (DIC)
		+ (params.bec.Qfe_zoo * zoo_loss_doc) ...  % zoo mortality to Fe (DOC)
		+ (params.bec.Qfe_zoo * szoo_loss_dic) ... % small zoo mortality to Fe (DIC)
		+ (params.bec.Qfe_zoo * szoo_loss_doc) ... % small zoo mortality to Fe (DOC)
		- (sum(photoFe)) ...                       % Fe -- > autotroph biomass
		- (Fe_scavenge);                           % Fe scavenging by particles
	% ... include uptake, losses, and excretion via chemos
	for i = 1:length(opt.chemo_ind);
		sms.fe(k,1) = sms.fe(k,1) ...
			- (chemo_growth(i) / params.(opt.chemo_ind{i}).CFe) ...    % chemo uptake of Fe
			+ (chemo_loss_dic(i) / params.(opt.chemo_ind{i}).CFe) ...  % chemo mortality to Fe (DIC) 
			+ (chemo_graze_dic(i) / params.(opt.chemo_ind{i}).CFe) ... % chemo grazing to Fe (DIC)
			+ (chemo_loss_doc(i) / params.(opt.chemo_ind{i}).CFe) ...  % chemo mortality to Fe (DOC) 
			+ (chemo_graze_doc(i) / params.(opt.chemo_ind{i}).CFe) ... % chemo grazing to Fe (DOC)
			+ (chemo_graze_fe(i));                                     % chemo grazing to Fe (Fe-rich) 
	end
	% ... add heterootroph losses and excretion to Fe
	for i = 1:length(opt.hetero_ind)
		sms.fe(k,1) = sms.fe(k,1) ...
			- (hetero_growth(i) / params.(opt.hetero_ind{i}).CFe) ...    % hetero uptake of Fe
			+ (hetero_loss_dic(i) / params.(opt.hetero_ind{i}).CFe) ...  % hetero mortality to Fe (DIC)
			+ (hetero_graze_dic(i) / params.(opt.hetero_ind{i}).CFe) ... % hetero grazing to Fe (DIC)
			+ (hetero_loss_doc(i) / params.(opt.hetero_ind{i}).CFe) ...  % hetero mortality to Fe (DIC)
			+ (hetero_graze_doc(i) / params.(opt.hetero_ind{i}).CFe) ... % hetero grazing to Fe (DIC)
			+ (hetero_graze_fe(i));                                      % hetero grazing to Fe (Fe-rich)
	end
	% ... add excretion from autotrophs
	for i = 1:length(opt.auto_ind);
		sms.fe(k,1) = sms.fe(k,1) ...
			+ (Qfe(i) * auto_loss_dic(i)) ...           % zoo mortality to Fe (DIC)
			+ (Qfe(i) * auto_graze_dic(i)) ...          % zoo grazing to Fe (DIC)
			+ (Qfe(i) * auto_loss_doc(i)) ...           % zoo loss to Fe (DOC)
			+ (Qfe(i) * auto_graze_doc(i)) ...          % zoo grazing to Fe (DIC)
			+ (Qfe(i) * auto_graze_zoo(i)) ...          % difference between Fe uptake from zoo growth 
			- (params.bec.Qfe_zoo * auto_graze_zoo(i)); % and any altered Fe stoichiometry of autotrophs
	end    

	% ---- %
    % SiO3 %
	% ---- %
	sms.sio3(k,1) = P_SiO2_remin; % particulate remineralization to SiO3
	% ... add uptake and release via phytoplankton
	for i = 1:length(opt.auto_ind)
		if params.(opt.auto_ind{i}).K_sio3 > 0
			sms.sio3(k,1) = sms.sio3(k,1) ...
				- photoSi(i) ...                                                     % SiO3 --> autotroph biomass
				+ Qsi(i) * (params.bec.f_graze_si_remin * auto_graze(i)) ...         % auto grazing to SiO3
				+ Qsi(i) * ((1 - params.(opt.auto_ind{i}).loss_poc) * auto_loss(i)); % auto mortality to SiO3
		end
	end

	% ------------------------ %
    % Dissolved organic matter %
	% ------------------------ %
    sms.doc(k,1)  = (doc_prod * (1 - params.bec.docrefract)) - doc_remin;               % production - remineralization
    sms.docr(k,1) = (doc_prod * (0 + params.bec.docrefract)) - docr_remin;              % production - remineralization
    sms.don(k,1)  = (don_prod * (1 - params.bec.donrefract)) - don_remin;               % production - remineralization
    sms.donr(k,1) = (don_prod * (0 + params.bec.donrefract)) - donr_remin;              % production - remineralization
    sms.dop(k,1)  = (dop_prod * (1 - params.bec.doprefract)) - dop_remin - sum(DOP_V);  % production - remineralization
    sms.dopr(k,1) = (dop_prod * (0 + params.bec.doprefract)) - dopr_remin;              % production - remineralization

	% -------------------------- %
    % Dissolved inorganic carbon %
	% -------------------------- %
    sms.dic(k,1) = 0 ...
         + sum(auto_loss_dic) ...    % auto mortality to DIC
         + sum(auto_graze_dic) ...   % auto grazing to DIC
         - sum(photoC) ...           % auto uptake of DIC
         + P_CaCO3_remin ...         % CacO3 remin to DIC
         + zoo_loss_dic ...          % zoo mortality to DIC
		 + szoo_loss_dic ...         % szoo mortality to DIC
		 + sum(hetero_graze_dic)...  % hetero grazing to DIC
         + sum(chemo_loss_dic) ...   % chemo mortality to DIC
         + sum(chemo_graze_dic) ...  % chemo grazing to DIC
         + sum(hetero_loss_dic);     % hetero mortality to DIC
    % ... add production from heterotroph growth 
    for i = 1:length(opt.hetero_ind)
        sms.dic(k,1) = sms.dic(k,1) ...
            + (hetero_growth(i) * (1 - facultative(i)) * (params.(opt.hetero_ind{i}).e_dic_aer)) ... % DIC excretion (if aerobic)
            + (hetero_growth(i) * (0 + facultative(i)) * (params.(opt.hetero_ind{i}).e_dic_ana));    % DIC excretion (if anaerobic)
    end
	% ... add consumption from chemoautotrophic growth
	for i = 1:length(opt.chemo_ind)
		sms.dic(k,1) = sms.dic(k,1) ...
			- (chemo_growth(i) * params.(opt.chemo_ind{i}).y_dic);
	end
    % ... add autotroph production and consumption
    for i = 1:length(opt.auto_ind);
        if params.(opt.auto_ind{i}).CaCO3_ind
            sms.dic(k,1) = sms.dic(k,1) ...
                + params.bec.f_graze_CaCO3_remin * auto_graze(i) * QCaCO3(i) ... % production due to CaCO3 remin
                - caco3_prod(i);                                                 % consumption due to CaCO3 prod
        end
    end

	% ---------- %
    % Alkalinity %
	% ---------- %
    sms.alk(k,1) = 0 ...
        + sms.nh4(k,1) ...        % NH4 contribution (+)
        + (2 * P_CaCO3_remin) ... % CaCO3 contribution (+)
        - sms.no3(k,1) ...        % NO3 contribution (-)
        - sms.no2(k,1);           % NO2 contribution (-)
    % ... add autotroph production and consumption
    for i = 1:length(opt.auto_ind);
        if params.(opt.auto_ind{i}).CaCO3_ind
            sms.alk(k,1) = sms.alk(k,1) ...
                + 2 * (params.bec.f_graze_CaCO3_remin * auto_graze(i) * QCaCO3(i)) ... % CaCO3 remin (+)
                - 2 * (caco3_prod(i));                                                 % CaCO3 prod (-)
        end
    end

	% ------ %
    % Oxygen %
	% ------ %
    % Gather production during autotroph biomass growth
    o2_production = 0;
    for i = 1:length(opt.auto_ind)
        if params.(opt.auto_ind{i}).Nfixer 
            if photoC(i) > 0
                o2_production = o2_production ...
                    + photoC(i) * (...
                    ((NO3_V(i) / (NO3_V(i) + NO2_V(i) + NH4_V(i) + Nfix(i))) / params.bec.Red_D_C_O2) + ...      % photosynthesis with NO3
                    ((NO2_V(i) / (NO3_V(i) + NO2_V(i) + NH4_V(i) + Nfix(i))) / params.bec.Red_D_C_O2_NO2V) + ... % photosynthesis with NO2
                    ((NH4_V(i) / (NO3_V(i) + NO2_V(i) + NH4_V(i) + Nfix(i))) / params.bec.Red_D_C_O2) + ...      % photosynthesis with NH4
                    ((Nfix(i)  / (NO3_V(i) + NO2_V(i) + NH4_V(i) + Nfix(i))) / params.bec.Red_D_C_O2_diaz) ...   % photosynthesis with N2
                    );
            end
        else
            if photoC(i) > 0
                o2_production = o2_production ...
                    + photoC(i) * (...
                    ((NO3_V(i) / (NO3_V(i) + NO2_V(i) + NH4_V(i))) / params.bec.Red_D_C_O2) + ...      % photosynthesis with NO3
                    ((NO2_V(i) / (NO3_V(i) + NO2_V(i) + NH4_V(i))) / params.bec.Red_D_C_O2_NO2V) + ... % photosynthesis with NO2
                    ((NH4_V(i) / (NO3_V(i) + NO2_V(i) + NH4_V(i))) / params.bec.Red_D_C_O2)...         % photosynthesis with NH4
                    );
            end
        end
    end
    % Gather consumption due to remineralization
	% Losses to DIC follow the approach in Paulmier (2009), 3.1.1, to calculate O2 required for instant remin
	% OM requires 138 mmol O2/m3
	% Bacterial B requires 56.25 mmol O2/m3
    o2_consumption = 0 ...
        + (sum(auto_loss_dic) / params.bec.Remin_D_C_O2) ...    % O2 loss via auto mortality
        + (sum(auto_graze_dic) / params.bec.Remin_D_C_O2) ...   % O2 loss via auto grazing
        + (sum(chemo_loss_dic) / (55/56.25)) ...                % O2 loss via chemo mortality
        + (sum(hetero_loss_dic) / (55/56.25)) ...               % O2 loss via hetero mortality
        + (sum(chemo_graze_dic) / (55/56.25)) ...               % O2 loss via chemo grazing
		+ (sum(hetero_graze_dic) / (55/56.25)) ...              % O2 loss via hetero grazing
        + (zoo_loss_dic / params.bec.Remin_D_C_O2) ...          % O2 loss via zoo mortality
		+ (szoo_loss_dic / params.bec.Remin_D_C_O2);            % O2 loss via szoo mortality
    % ... add consumption via chemo uptake
	o2_consumption = o2_consumption ...
		+ (chemo_growth(opt.aoa_ind) * params.aoa.y_oxy) ... % aoa O2 yield 
		+ (chemo_growth(opt.nob_ind) * params.nob.y_oxy);    % nob O2 yield 
    % ... add consumption via heterotroph uptake
    for i = 1:length(opt.hetero_ind)
        o2_consumption = o2_consumption ...
            + (hetero_growth(i) * (1 - facultative(i)) * params.(opt.hetero_ind{i}).y_oxy_aer); % hetero O2 yield (if aerobic) 
    end

    % ... get difference
    sms.o2(k,1) = o2_production - o2_consumption;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute sources and sinks for autotrophs and their cellular components 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:length(opt.auto_ind);
        work1 = auto_graze(i) + auto_loss(i) + auto_agg(i);                      % loss terms
        % BIO carbon
        sms.(opt.auto_ind{i})(k,1) = photoC(i) - work1;                          % gain - loss (carbon)
        % BIO chlA
        sms.([opt.auto_ind{i},'_chl'])(k,1) = photoacc(i) - (thetaC(i) * work1); % gain - loss (ChlA)
        % BIO fe
        sms.([opt.auto_ind{i},'_fe'])(k,1) = photoFe(i) - (Qfe(i) * work1);      % gain - loss (Fe)
        % BIO si
        if params.(opt.auto_ind{i}).Si_ind
            sms.([opt.auto_ind{i},'_si'])(k,1) = photoSi(i) - (Qsi(i) * work1);  % gain - loss (Si)
        else
            sms.([opt.auto_ind{i},'_si'])(k,1) = 0;                                       
        end
        % BIO caco3
        if params.(opt.auto_ind{i}).CaCO3_ind
            sms.([opt.auto_ind{i},'_caco3'])(k,1) = caco3_prod(i) - (QCaCO3(i) * work1); % gain - loss (CaCO3)
        else
            sms.([opt.auto_ind{i},'_caco3'])(k,1) = 0;
        end
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute sources and sinks for chemoautotrophs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    for i = 1:length(opt.chemo_ind)
        sms.(opt.chemo_ind{i})(k,1) = 0 ...
            + (chemo_growth(i)) ... % chemo growth
            - (chemo_graze(i)) ...  % chemo grazing
            - (chemo_loss(i));      % chemo mortality
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute sources and sinks for heterotrophs 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    for i = 1:length(opt.hetero_ind)
        sms.(opt.hetero_ind{i})(k,1) = 0 ...
            + (hetero_growth(i)) ... % hetero growth
            - (hetero_loss(i)) ...   % hetero mortality
			- (hetero_graze(i));     % hetero grazing
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute sources and sinks for zooplankton 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sms.zoo(k,1)  = 0 ...
        + (sum(auto_graze_zoo)) ...  % zoo growth via auto grazing
        + (sum(chemo_graze_zoo)) ... % zoo growth via chemo grazing
        - (zoo_loss);                % zoo mortality
	sms.szoo(k,1) = 0 ...;
		+ (sum(hetero_graze_szoo)) ... % szoo growth via hetero grazing
		- (szoo_loss);                 % szoo mortality

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Add any diagnostics below 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save diagnostics?
    if diag_out == 1;
        % N-Cycle reactions
        diag.ammox(k,1)      = chemo_growth(opt.aoa_ind) * params.aoa.y_red;                                 % aoa NH4 uptake rate (uM N/s)
        diag.nitrox(k,1)     = chemo_growth(opt.nob_ind) * params.nob.y_red;                                 % nob NO2 uptake rate (uM N/s)
        diag.anammox(k,1)    = chemo_growth(opt.aox_ind) * params.aox.y_red;                                 % aox NH4 uptake rate (uM N/s)
		diag.o2_remin(k,1)   = hetero_growth(opt.aer_ind) * params.aer.y_oxy_aer;                            % aer O2 uptake rate  (uM O2/s)
        diag.denitrif1(k,1)  = hetero_growth(opt.nar_ind) * facultative(opt.nar_ind) * params.nar.y_oxy_ana; % nar NO3 uptake rate (uM N/s)
        diag.denitrif2(k,1)  = hetero_growth(opt.nir_ind) * facultative(opt.nir_ind) * params.nir.y_oxy_ana; % nir NO2 uptake rate (uM N/s)
        diag.denitrif3(k,1)  = hetero_growth(opt.nos_ind) * facultative(opt.nos_ind) * params.nos.y_oxy_ana; % nos N2O uptake rate (uM N2/s)
        diag.denitrif4(k,1)  = hetero_growth(opt.nai_ind) * facultative(opt.nai_ind) * params.nai.y_oxy_ana; % nai NO3 uptake rate (uM N/s)
        diag.denitrif5(k,1)  = hetero_growth(opt.nio_ind) * facultative(opt.nio_ind) * params.nio.y_oxy_ana; % nio NO2 uptake rate (uM N/s)
        diag.denitrif6(k,1)  = hetero_growth(opt.nao_ind) * facultative(opt.nao_ind) * params.nao.y_oxy_ana; % nao NO3 uptake rate (uM N/s)
		if opt.DNRA1
			diag.dnra_no3(k,1) = hetero_growth(opt.dnra1_ind) * facultative(opt.dnra1_ind) * params.dnra1.y_oxy_ana; % dnra1 NO3 uptake rate (uM N/s)
		end
		if opt.DNRA2
			diag.dnra_no2(k,1) = hetero_growth(opt.dnra2_ind) * facultative(opt.dnra2_ind) * params.dnra2.y_oxy_ana; % dnra1 NO3 uptake rate (uM N/s)
		end

        % Reactants 
        diag.o2ammox(k,1)    = chemo_growth(opt.aoa_ind) * params.aoa.y_oxy; % aoa O2 uptake rate (uM O2/s) 
        diag.o2nitrox(k,1)   = chemo_growth(opt.nob_ind) * params.nob.y_oxy; % nob O2 uptake rate (uM O2/s)
        diag.nh4nitrox(k,1)  = chemo_growth(opt.nob_ind) * params.nob.y_nh4; % nob NH4 uptake rate (uM N/s)
        diag.no2anammox(k,1) = chemo_growth(opt.aox_ind) * params.aox.y_oxy; % aox NO2 uptake rate (uM N/s)

        % Products
        diag.no2ammox(k,1)   = chemo_growth(opt.aoa_ind) * params.aoa.e_no2 * yaoa_no2;           % aoa NO2 excretion rate (uM N/s)
        diag.hydroxln2o(k,1) = chemo_growth(opt.aoa_ind) * params.aoa.e_no2 * 0.5 * yaoa_n2o_nh4; % aoa N2O excretion rate (uM N/s)
        diag.hybridn2o(k,1)  = chemo_growth(opt.aoa_ind) * params.aoa.e_no2 * 0.5 * yaoa_n2o_hyb; % aoa N2O excretion rate (uM N/s)
        diag.no3nitrox(k,1)  = chemo_growth(opt.nob_ind) * params.nob.e_no3;                      % nob NO3 excretion rate (uM N/s)
        diag.no3anammox(k,1) = chemo_growth(opt.aox_ind) * params.aox.e_no3;                      % aox NO3 excretion rate (uM N/s)
        diag.n2anammox(k,1)  = chemo_growth(opt.aox_ind) * params.aox.e_n2;                       % aox N2 excretion rate (uM N2/s)

        % Fixed N uptake
        diag.spNO3uptake(k,1)   = NO3_V(opt.sp_ind);   % small phyto NO3 uptake rate (uM N/s)
        diag.spNO2uptake(k,1)   = NO2_V(opt.sp_ind);   % small phyto NO2 uptake rate (uM N/s)
        diag.spNH4uptake(k,1)   = NH4_V(opt.sp_ind);   % small phyto NH4 uptake rate (uM N/s)
        diag.diatNO3uptake(k,1) = NO3_V(opt.diat_ind); % diatom NO3 uptake rate (uM N/s)
        diag.diatNO2uptake(k,1) = NO2_V(opt.diat_ind); % diatom NO2 uptake rate (uM N/s)
        diag.diatNH4uptake(k,1) = NH4_V(opt.diat_ind); % diatom NH4 uptake rate (uM N/s)
        diag.diazNO3uptake(k,1) = NO3_V(opt.diaz_ind); % diazotroph NO3 uptake rate (uM N/s)
        diag.diazNO2uptake(k,1) = NO2_V(opt.diaz_ind); % diazotroph NO2 uptake rate (uM N/s)
        diag.diazNH4uptake(k,1) = NH4_V(opt.diaz_ind); % diazotroph NH4 uptake rate (uM N/s)

        % Organic matter
        diag.poc_prod(k,1)    = poc_prod;                    % POC production (uM C/s)
        diag.poc_flux_in(k,1) = POC_sflux_in + POC_hflux_in; % POC flux into cell (uM C/s)
        diag.poc_remin(k,1)   = poc_remin;                   % POC remineralization (uM C/s)
        diag.doc_prod(k,1)    = doc_prod;                    % DOC production (uM C/s)
        diag.doc_remin(k,1)   = doc_remin;                   % DOC remineralization (uM C/s)
        diag.docr_remin(k,1)  = docr_remin;                  % DOCr remineralization (uM C/s)
        diag.don_prod(k,1)    = don_prod;                    % DON production (uM N/s)        
        diag.don_remin(k,1)   = don_remin;                   % DON remineralization (uM N/s)
        diag.donr_remin(k,1)  = donr_remin;                  % DONr remineralization (uM N/s)
        diag.dop_prod(k,1)    = dop_prod;                    % DOP production (uM P/s)       
        diag.dop_remin(k,1)   = dop_remin;                   % DOP remineralization (uM P/s)
        diag.dopr_remin(k,1)  = dopr_remin;                  % DOPr remineralization (uM P/s)

        % PAR
        diag.par_out(k,1) = par_out; % PAR, bottom of cell (W/m2)
        diag.par_in(k,1)  = par_in;  % PAR, top of cell (W/m2)
        diag.par_lay(k,1) = par_lay; % PAR, cell average (W/m2)

        % Oxygen
        diag.o2_production(k,1)  = o2_production;  % Gross O2 production (uM O2/s)
        diag.o2_consumption(k,1) = o2_consumption; % Gross O2 consumption (uM O2/s)

        % Chemoautotrophs
        for i = 1:length(opt.chemo_ind);
            diag.([opt.chemo_ind{i},'_growth'])(k,1) = chemo_growth(i); % chemo growth
            diag.([opt.chemo_ind{i},'_graze'])(k,1)  = chemo_graze(i);  % chemo grazing 
            diag.([opt.chemo_ind{i},'_loss'])(k,1)   = chemo_loss(i);   % chemo mortality
        end
		diag.chemo_growth(k,1) = sum(chemo_growth);
		diag.chemo_graze(k,1)  = sum(chemo_graze);
		diag.chemo_loss(k,1)   = sum(chemo_loss);

        % Heterotrophs    
        for i = 1:length(opt.hetero_ind);
            diag.([opt.hetero_ind{i},'_growth'])(k,1) = hetero_growth(i); % hetero growth
            diag.([opt.hetero_ind{i},'_loss'])(k,1)   = hetero_loss(i);   % hetero mortality
			diag.([opt.hetero_ind{i},'_graze'])(k,1)  = hetero_graze(i);  % hetero grazing
        end
		diag.hetero_growth(k,1) = sum(hetero_growth);
		diag.hetero_graze(k,1)  = sum(hetero_graze);
		diag.hetero_loss(k,1)   = sum(hetero_loss);

		% Others
		diag.zoo_loss_doc(k,1)     = zoo_loss_doc;          % zoo mortality to doc
		diag.auto_loss_doc(k,1)    = sum(auto_loss_doc);    % auto mortality to doc
		diag.auto_graze_doc(k,1)   = sum(auto_graze_doc);   % auto grazing to doc
		diag.chemo_loss_doc(k,1)   = sum(chemo_loss_doc);   % chemo mortality to doc
		diag.chemo_graze_doc(k,1)  = sum(chemo_graze_doc);  % chemo grazing to doc
		diag.hetero_loss_doc(k,1)  = sum(hetero_loss_doc);  % hetero mortality to doc
		diag.hetero_graze_doc(k,1) = sum(hetero_graze_doc); % hetero grazing to doc
		diag.szoo_loss_doc(k,1)    = szoo_loss_doc;         % small zoo mortality to doc
		diag.sp_graze(k,1)         = auto_graze(1);         % small phyto grazing to zoo bio
		diag.diat_graze(k,1)       = auto_graze(2);         % diatom grazing to zoo bio
		diag.diaz_graze(k,1)       = auto_graze(3);         % diaz grazing to zoo bio
    else
        % Blank structure
        diag = struct;
    end
end % end k-loop

