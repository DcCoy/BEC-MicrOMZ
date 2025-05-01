# BEC-MicrOMZ 
BEC-MicrOMZ is a biogeochemical ocean model which is based on the Biogeochemical Elemental Cycling (BEC) model as represented in UCLA ROMS (Moore et al., 2004; Shchepetkin & McWilliams, 2005; Deutsch et al., 2020).

Here, we expand the BEC module to explicitly represent the growth, mortality, and grazing of chemoautotrophic and heterotrophic bacteria functional type populations (herein termed 'MicrOMZ', see Zakem et al., 2019, 2022).

By modeling the metabolisms of these bacterial populations, we are able to resolve the key nitrogen cycle reactions of nitrification, denitrification, and anammox using Tillman resource competition arguments, rather than parameterizing the reactions as a function of seawater chemistry only (e.g., 'NitrOMZ' from Bianchi et al., 2022).

The hybrid BEC-MicrOMZ model is imbedded here in a 1-D advection-diffusion-reaction model representative of a typical open-ocean water column overlying an oxygen-minimum-zone (OMZ). 

Initial conditions, lateral restoring, and surface forcing are extracted from UCLA ROMS Eastern Tropical South Pacific (ETSP) simulations (McCoy et al., 2023).

Contact Daniel McCoy (dmccoy@carnegiescience.edu) for assistance in configuring BEC-MicrOMZ for other oceanographic regions/regimes.
    
## Table of Contents

- [Updates](#updates)
- [Getting started](#getting-started)
- [Code structure](#code-structure)
- [Support](#support)
- [How to cite](#how-to-cite)

Requires MATLAB 2013 or above.

## Updates
* 03/08/2024 -- First commit of BEC-MicrOMZ 
* 03/18/2024 -- First 'working' version finished
* 03/19/2024 -- Simplified air-sea flux module
* 03/20/2024 -- Added restoring of Fe and SiO3 to initial conditions
* 03/22/2024 -- Added routine to calculate growth yields
* 03/27/2024 -- Corrected stoichiometry of OM to match BEC 
* 04/04/2024 -- Corrected bug in chemoautotroph growth 
* 04/25/2024 -- Revived air-sea flux module, added additional switches to fix yOM for anaerobes 
* 04/29/2024 -- Fixed bugs with air-sea flux, changed equation for N2O yield from NH4 oxidation 
* 05/01/2025 -- Upload latest version 

## Getting started
#### Set run options via 'runscripts/options.m'
    Set model settings, such as:
        REGION                -- 'ETSP' or 'ETNP', to load initial conditions and forcing from UCLA ROMS simulations
        RESTART               -- Switch to start from a previous run, or from initial conditions 
        ADVECTION             -- Currently set to FTCS with UPWIND (add more methods later?) 
        DEPTHVAR_WUP          -- Switch to use constant or depth-dependent upwelling velocity
        DEPTHVAR_KV           -- Switch to use constant or depth-dependent diffusion profile
        CONSTANT_TSTEP        -- Switch to use constant or variable time-stepping (not coded yet)
        RESTORING             -- Switch to turn restoring on for all variables
        TAUZVAR               -- Switch to use a depth-dependent restoring time-scale
        HIST_VERBOSE          -- Prompts a message at each output timestep 
        REST_VERBOSE          -- Prompts a message when loading a restart file
        TDEP_REMIN            -- Switch to use temperature dependent implicit remineralization
        O2_REMIN              -- Switch to increase POC remineralization length scale at low O2 (5 -- 40)
        N2_FULL               -- Switch to use full N2 concentrations, or only excess N2 (produced locally)
        SEDIMENTS             -- Switch to use the sediment module in 'compute_particulate_terms'
        GASEXCHANGE           -- Switch to 'turn on' surface cell (advection, sms) and allow for gas exchange

    Set MicrOMZ model settings, such as:
        EXPLICIT_MICROBES     -- Switch to use MicrOMZ modeul (explicit heterotrophs and chemoautotrophs)
        FACULTATIVE_MICROBES  -- Switch to allow hetertrophs to switch between aerobic and anaerobic metabolisms
        HETERO_GRAZING        -- Switch to include small zooplankton predator that consumes heterotrophs
        HETERO_GRAZING_O2     -- Switch to limit heterotroph grazing at low O2
        CHEMO_GRAZING         -- Switch to allow chemoautotrophs to be consumed by zooplankton 
        CHEMO_GRAZING_O2      -- Switch to limit chemoautotroph grazing at low O2
        O2_PULSE              -- Switch to add eddy-fluxes of O2 in OMZ at distinct periods
        KELLY_N2O_YIELD       -- Switch to use combined N2O yields from Kelly et al., 2024 (AOA)
        DNRA1                 -- Switch to include DNRA types (starting at NO3)
        DNRA2                 -- Switch to include DNRA types (starting at NO2)
        SINSABAUGH            -- Switch to calculate heterotroph OM yield using approach of Sinsabaugh (2013) 
        SET_AER_F             -- Switch to set 'f' for obligate aerobic bacteria (aer) to literature values 
        FACULTATIVE_PENALTY   -- Switch to assign a penalty to facultative anaerobes during 'f' calculation 
        RAPID_EXCLUSION       -- Switch to override linear and quadratic mortality to encourage competitive exclusion
        HIGH_UMAX             -- Switch to raise default maximum growth rates (Buchanan et al., 2024)

    The user can also toggles settings for:
        Data sources          -- For upwelling, restoring profiles, restart files, initial conditions, and forcing  
        Restoring switches    -- Turn on/off which tracers to restore to farfield/initial profiles
        The vertical grid     -- top, bottom, and dz
        Time-stepping         -- dt, and rate of history output
        Tracers in the model  -- 28 for BEC, +12 if EXPLICIT_MICROBES is true
        Stoichiometry         -- Set C:H:O:N of organic matter and bacterial biomass if EXPLICIT_MICROBES is true

#### Run the model 
    Open up Matlab
    Change directory to 'runscripts'
    run_model

## Code structure 
#### runscripts/  
    Folder where BEC-MicrOMZ runscripts and run settings files are stored
        options.m             -- Script to toggle main model settings
        run_model.m           -- Call this to run the model ('>> run_model')
                   
#### functions/
    Folder where minor model functions are stored
        cschmidt_n2.m (n2o/o2)         -- Computes Schmidth number of various dissolved gas tracers
        n2satu.m (n2o/o2)              -- Computes atmospheric saturation concentrations of dissolved gas tracers
        extract_forcing.m              -- Script to extract forcing (archived on GitHub, but only used by D.McCoy)
        extract_initial_conditions.m   -- Script to extract initial conditions (archived on GitHub, but only used by D.McCoy)
        extract_ROMS_output            -- Script to extract year 50 solution from McCoy et al., 2023 (used only by D.McCoy)
        wind_stress                    -- Calculates wind stress based on forcing
        calc_f                         -- Calculates the fraction of electrons using Gibbs free energy arguments 
        get_f_from_Gibbs.m             -- Routine that calls 'calc_f' if opt.CALC_F is true
        yield_sinsabaugh.m             -- Function to calculate biomass yield on organic matter via Sinsabaugh et al., 2013
        diffusion.m                    -- Diffusion function
        upwind.m                       -- Vertical advection function

#### src
    Folder where core model functions are stored
        upwind.m                       -- 1D advection function
        diffusion.m                    -- 1D diffusion function
        sources_sinks.m                -- BEC-MicrOMZ BGC module
        timestepping.m                 -- Main module where timestepping takes place 
        airsea_flux                    -- Air-sea flux module for dissolved gases (O2, N2, N2O) <--opt.GASEXCHANGE 
        restoring.m                    -- Script to restore tracers to farfield concentrations <--opt.RESTORING
        compute_particulate_terms.m    -- Computes interior particulate pools based on Armstrong et al., 2000
        initialize_particulate_terms.m -- Initializes surface particulate pools
        initialize_forcing.m           -- Load surface forcing from ROMS output into 'frc' structure
        initialize_model.m             -- Get timestepping, parameters, vertical grid, and physical environment based on 'options.m'
        initialize_restoring.m         -- (WORK-IN-PROGRESS) Used to process restoring profiles and timescales
        initialize_timestepping.m      -- Get timestepping arrays based on options.m
        initialize_tracers.m           -- Script to load tracers from restart, or from initial conditions

#### plotting/ 
    Folder where some useful plotting functions for BEC-MicrOMZ output are stored

#### parameters/ 
    Folder where physical and biogeochemistry parameters live 
        BEC_overrides.m        -- Add any changes to default BEC parameters here
        MicrOMZ_overrides.m    -- Add any changes to default MicrOMZ parameters here
        get_BEC_params.m       -- Default BEC parameters
        get_MicrOMZ_params.m   -- Default MicrOMZ parameters 
        get_Vmax.m             -- Calculate substrate uptake based on max growth and yields 
        get_diffusive_uptake.m -- Overrides pcoef for bacteria if opt.CALC_PCOEF is true
        get_f_from_Gibbs.m     -- Calculate yields based on thermodynamics of redox-reactions
        get_growth_yields.m    -- Calculates growth yields based on fraction 'f' of electrons used in biomass synthesis
        get_physical_params.m  -- Default physical environment parameters
        get_stoichiometry.m    -- Gets C:H:O:N ratios, normalized to C, based on options.m
        physical_overrides.m   -- Add any changes to default physical parameters here
        
#### data/
    Region-specific forcing and validation data
        compilation_ETNP_gridded.mat -- ETNP observed tracers, for validation 
        compilation_ETSP_gridded.mat -- ETSP observed tracers, for validation
        comprates_ETSP.mat           -- ETSP observed N transformation rates, for validation
        farfield_ETNP_gridded.mat    -- ETNP lateral restoring tracer concentrations 
        farfield_ETNP_gridded.mat    -- ETSP lateral restoring tracer concentrations
        Tau_restoring.mat            -- Region-specific restoring time-scales
        vertical_CESM.mat            -- Region-specific vertical upwelling data
        ETSP_ROMS_forcing.mat        -- ROMS annual surface forcing extracted via 'functions/extract_forcing.m'
        ETSP_ROMS_restart.mat        -- ROMS initial conditions extracted via 'functions/extract_initial_conditions.m'
        ETSP_ROMS_output.mat         -- ROMS year 50 solution at profile location, for validation purposes (McCoy et al., 2023)

#### output/
    Folder to store output (and restart files)
      
## Support
Contact Daniel McCoy at the Carnegie Institution for Science (dmccoy@carnegiescience.edu) 

