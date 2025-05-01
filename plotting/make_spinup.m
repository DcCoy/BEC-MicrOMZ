% List fname
clear all
fdir  = ['../output/'];
fname = ['MicrOMZ_exclude_control_50years.mat'];

% Make directory
cmd = ['mkdir plots/',fname(1:end-4)]; system(cmd);

% OPTIONS
% Tracer switches
make_bgc    = true;
make_bio    = true;
make_chemo  = true;
make_hetero = true;
make_OM     = true;
make_custom = false;
% Diagnostic switches
make_chemo_rates  = true;
make_hetero_rates = true;
make_diag        = false;
make_diag_custom = false; 

% Parallel?
parallel_pool = true;
if parallel_pool
    delete(gcp('nocreate'));
    parpool(12);
end

% ------------------------------- %
% Tracers
% ------------------------------- %
% BGC
if make_bgc
    outdir = [fname(1:end-4),'/bgc'];
    vars   = {'o2_avg','nh4_avg','no3_avg','no2_avg','n2o_avg','n2_avg'};
    tits   = {'O$_2$','NH$^{+}_4$','NO$^{-}_3$','NO$^{-}_2$','N$_2$O','N$_2$'};
    if parallel_pool
        plot_tracer_spinup_parallel([fdir,fname],outdir,vars,tits);
    else
        plot_tracer_spinup([fdir,fname],outdir,vars,tits);
    end
end

% BIO
if make_bio
    outdir = [fname(1:end-4),'/bio'];
    vars   = {'sp_avg','diat_avg','diaz_avg','zoo_avg','szoo_avg'};
    tits   = {'Small phyto','Diatoms','Diazotrophs','Large zoo','Small zoo'};
    if parallel_pool
        plot_tracer_spinup_parallel([fdir,fname],outdir,vars,tits);
    else
        plot_tracer_spinup([fdir,fname],outdir,vars,tits);
    end
end

% CHEMOAUTOTROPHS
if make_chemo
    outdir = [fname(1:end-4),'/chemo'];    
    vars   = {'aoa_avg','nob_avg','aox_avg'};
    tits   = {'AOA','NOB','AOX'};
    if parallel_pool
        plot_tracer_spinup_parallel([fdir,fname],outdir,vars,tits);
    else
        plot_tracer_spinup([fdir,fname],outdir,vars,tits);
    end
end

% HETEROTROPHS
if make_hetero
    outdir = [fname(1:end-4),'/hetero'];
    vars = {'aer_avg','nar_avg','nai_avg','nao_avg','nir_avg','nio_avg','nos_avg'};
    tits = {'AER','NAR','NAI','NAO','NIR','NIO','NOS'};
    if parallel_pool
        plot_tracer_spinup_parallel([fdir,fname],outdir,vars,tits);
    else
        plot_tracer_spinup([fdir,fname],outdir,vars,tits);
    end
end

% Organic Matter
if make_OM
    outdir = [fname(1:end-4),'/OM'];
    vars = {'doc_avg','docr_avg','don_avg','donr_avg','dop_avg','dopr_avg'};
    tits = {'DOC','DOCr','DON','DONr','DOP','DOPr'};
    if parallel_pool
        plot_tracer_spinup_parallel([fdir,fname],outdir,vars,tits);
    else
        plot_tracer_spinup([fdir,fname],outdir,vars,tits);
    end
end

% Custom
if make_custom
    outdir = [fname(1:end-4),'/custom'];
    %vars = {'o2','o2_adv','o2_dfz','o2_sms'};
    %tits = {'O$_2$','O$_2$-adv','O$_2$-dfz','O$_2$-sms'};
    vars = {'fe_avg','po4_avg','sio3_avg'};
    tits = {'Fe','PO$^{3-}_4$','SiO$^{2-}_3$'};
    if parallel_pool
        plot_tracer_spinup_parallel([fdir,fname],outdir,vars,tits);
    else
        plot_tracer_spinup([fdir,fname],outdir,vars,tits);
    end
end

% ------------------------------- %
% Diagnostics 
% ------------------------------- %
% Chemo rates
if make_chemo_rates
    outdir = [fname(1:end-4),'/chemo_rates'];
    vars = {'ammox_avg','nitrox_avg','anammox_avg'};
    tits = {'AMMOX','NITROX','ANAMMOX'};
    if parallel_pool
        plot_diag_spinup_parallel([fdir,fname],outdir,vars,tits);
    else
        plot_diag_spinup([fdir,fname],outdir,vars,tits);
    end
end

% Hetero rates
if make_hetero_rates
    outdir = [fname(1:end-4),'/hetero_rates'];
    vars = {'denitrif1_avg','denitrif2_avg','denitrif3_avg','denitrif4_avg','denitrif5_avg','denitrif6_avg'};
    tits = {'DENITRIF1','DENITRIF2','DENITRIF3','DENITRIF4','DENITRIF5','DENITRIF6'};
    if parallel_pool
        plot_diag_spinup_parallel([fdir,fname],outdir,vars,tits);
    else
        plot_diag_spinup([fdir,fname],outdir,vars,tits);
    end
end

% General diagnostics
if make_diag
    outdir = [fname(1:end-4),'/diag'];
    vars = {'poc_prod_avg','poc_remin_avg','poc_flux_in_avg','doc_prod_avg','doc_remin_avg'};
    tits = {'POC prod','POC remin','POC flux','DOC production','DOC remin'};
    if parallel_pool
        plot_diag_spinup_parallel([fdir,fname],outdir,vars,tits);
    else
        plot_diag_spinup([fdir,fname],outdir,vars,tits);
    end
end

% Custom diagnostics
if make_diag_custom
    outdir = [fname(1:end-4),'/diag_custom'];
    vars   = {'nos_growth_avg','nos_graze_avg','nos_loss_avg','nos_adv_avg','nos_dfz_avg','nos_sms_avg'};
    tits   = {'NOS growth','NOS grazing','NOS mortality','NOS advection','NOS diffusion','NOS SMS'};
    if parallel_pool
        plot_diag_spinup_parallel([fdir,fname],outdir,vars,tits);
    else
        plot_diag_spinup([fdir,fname],outdir,vars,tits);
    end
end

% ------------------------------- %
% Finish 
% ------------------------------- %
if parallel_pool
   delete(gcp('nocreate'));
end
   
