% List fname
clear all
fname = ['../output/MicrOMZ_facultative-10-Jun-2024.mat'];
fname = ['../output/MicrOMZ_facultative_low_K_doc-11-Jun-2024.mat'];
fname = ['../output/MicrOMZ_facultative_low_Km-13-Jun-2024.mat'];

% OPTIONS
% Tracer switches
make_bgc    = true;
make_bio    = true;
make_chemo  = true;
make_hetero = true;
make_OM     = true;
% Diagnostic switches
make_ncycle      = true;
make_diag        = true;
make_diag_custom = true; 

% ------------------------------- %
% Tracers
% ------------------------------- %
% BGC
if make_bgc
   outdir = ['final_bgc'];
   vars   = {'o2','nh4','no3','no2','n2o','n2'};
   tits   = {'O$_2$','NH$^{+}_4$','NO$^{-}_3$','NO$^{-}_2$','N$_2$O','N$_2$'};
   plot_tracer(fname,outdir,vars,tits);
end

% BIO
if make_bio
   outdir = ['final_bio'];
   vars   = {'sp','diat','diaz','zoo','szoo'};
   tits   = {'Small phyto','Diatoms','Diazotrophs','Large zoo','Small zoo'};
   plot_tracer(fname,outdir,vars,tits);
end

% CHEMOAUTOTROPHS
if make_chemo
   outdir = ['final_chemo'];   
   vars   = {'aoa','nob','aox'};
   tits   = {'AOA','NOB','AOX'};
   plot_tracer(fname,outdir,vars,tits);
end

% HETEROTROPHS
if make_hetero
   outdir = ['final_hetero'];
   vars = {'aer','nar','nai','nao','nir','nio','nos'};
   tits = {'AER','NAR','NAI','NAO','NIR','NIO','NOS'};
   plot_tracer(fname,outdir,vars,tits);
end

% Organic Matter
if make_OM
   outdir = ['final_OM'];
   vars = {'doc','docr','don','donr','dop','dopr'};
   tits = {'DOC','DOCr','DON','DONr','DOP','DOPr'};
   plot_tracer(fname,outdir,vars,tits);
end

% ------------------------------- %
% Diagnostics 
% ------------------------------- %
% Ncycle
if make_ncycle
   outdir = ['final_ncycle'];
   vars = {'ammox','nitrox','anammox','denitrif1','denitrif2','denitrif3','denitrif4','denitrif5','denitrif6'};
   tits = {'AMMOX','NITROX','ANAMMOX','DENITRIF1','DENITRIF2','DENITRIF3','DENITRIF4','DENITRIF5','DENITRIF6'};
   plot_diag(fname,outdir,vars,tits);
end

% General diagnostics
if make_diag
   outdir = ['final_diag'];
   vars = {'poc_remin','doc_prod','o2_production','o2_consumption'};
   tits = {'POC remin','DOC production','O$_2$ production','O$_2$ consumption'};
   plot_diag(fname,outdir,vars,tits);
end

% Custom diagnostics
if make_diag_custom
   outdir = ['final_diag_custom'];
   vars   = {'aoa_growth','aoa_graze','aoa_loss'};
   tits   = {'AOA growth','AOA grazing','AOA mortality'};
   plot_diag(fname,outdir,vars,tits);
end
