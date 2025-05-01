function params = get_growth_yields(params,opt)
% Calculates biomass and product yields based on 
% fraction (f) of electrons used for biomass synthesis
% and stoichiometry of OM and bacterial biomass
% NOTE: OM is already normalized to 1 mol C

% Extract stoichiometry
OM = params.stoich.OM;
BM = params.stoich.BM;

% Get ratios
% Organic matter
dOM = params.stoich.dOM;
cOM = params.stoich.cOM;
hOM = params.stoich.hOM;
oOM = params.stoich.oOM;
nOM = params.stoich.nOM;
pOM = params.stoich.pOM;

% Bacteria
dBM = params.stoich.dBM;
cBM = params.stoich.cBM;
hBM = params.stoich.hBM;
oBM = params.stoich.oBM;
nBM = params.stoich.nBM;
pBM = params.stoich.pBM;

% Use Sinsabaugh (2013) approach to OM yield of heterotrophs?
if opt.SINSABAUGH
    % Parameters used in approach
    y_max  = 0.6;           % maximum possible growth efficiency
    B_CN   = params.aer.CN; % C:N of bacterial biomass
    OM_CN  = cOM/nOM;       % C:N of organic matter
    K_CN   = 0.5;           % C:N half-saturation coefficient (Sinsabaugh & Follstad Shah 2012)
    EEA_CN = 1.123;         % Relative rate of enzymatic processing of complex C and complex N molecules
                            % into simple precursors for biosynthesis (Sinsabaugh & Follstad Shah 2012)
end

% ---------------------- %
% Heterotrophic bacteria %
% ---------------------- %

% ---------------------- %
% Obligate aerobic (AER) %
% ---------------------- %
% Reactants: OM + O2
% Products:  NH4 + DIC + BM
if opt.SINSABAUGH
    yOM_aer = yield_sinsabaugh(y_max,B_CN,OM_CN,K_CN,EEA_CN);
    f       = yOM_aer * (dBM/dOM);
else
    f = params.aer.f;
    yOM_aer  = (dBM/f)*(1/dOM);
end
yO2_aer  = (dBM/f)*((1-f)/4);
eNH4_aer = (dBM/f)*((nOM/dOM)-((nBM*f)/dBM)); 
eDIC_aer = (dBM/f)*((cOM/dOM)-((cBM*f)/dBM));
ePO4_aer = (dBM/f)*((pOM/dOM)-((pBM*f)/dBM));

% ---------------------------- %
% NO3 reduction to NO2 (DEN 1) %
% ---------------------------- %
% Reactants: OM + NO3
% Products:  NH4 + NO2 + DIC + BM
f = params.nar.f;
yOM_nar  = (dBM/f)*(1/dOM);
yNO3_nar = (dBM/f)*((1-f)/2);
eNO2_nar = (dBM/f)*((1-f)/2);
eNH4_nar = (dBM/f)*((nOM/dOM)-((nBM*f)/dBM));
eDIC_nar = (dBM/f)*((cOM/dOM)-((cBM*f)/dBM)); 
ePO4_nar = (dBM/f)*((pOM/dOM)-((pBM*f)/dBM));
% Facultative aerboic
f = params.nar.f_fac;
yOM_nar_fac  = (dBM/f)*(1/dOM);
yO2_nar_fac  = (dBM/f)*((1-f)/4);
eNH4_nar_fac = (dBM/f)*((nOM/dOM)-((nBM*f)/dBM)); 
eDIC_nar_fac = (dBM/f)*((cOM/dOM)-((cBM*f)/dBM));
ePO4_nar_fac = (dBM/f)*((pOM/dOM)-((pBM*f)/dBM));

% ---------------------------- %
% NO2 reduction to N2O (DEN 2) %
% ---------------------------- %
% Reactants: OM + NO2
% Products:  NH4 + N2O + DIC + BM
f = params.nir.f;
yOM_nir  = (dBM/f)*(1/dOM); 
yNO2_nir = (dBM/f)*((1-f)/2);
eN2O_nir = (dBM/f)*((1-f)/4);
eNH4_nir = (dBM/f)*((nOM/dOM)-((nBM*f)/dBM));
eDIC_nir = (dBM/f)*((cOM/dOM)-((cBM*f)/dBM));
ePO4_nir = (dBM/f)*((pOM/dOM)-((pBM*f)/dBM));
% Facultative aerboic
f = params.nir.f_fac;
yOM_nir_fac  = (dBM/f)*(1/dOM);
yO2_nir_fac  = (dBM/f)*((1-f)/4);
eNH4_nir_fac = (dBM/f)*((nOM/dOM)-((nBM*f)/dBM)); 
eDIC_nir_fac = (dBM/f)*((cOM/dOM)-((cBM*f)/dBM));
ePO4_nir_fac = (dBM/f)*((pOM/dOM)-((pBM*f)/dBM));

% --------------------------- %
% N2O reduction to N2 (DEN 3) %
% --------------------------- %
% Reactants: OM + N2O
% Products:  NH4 + N2 + DIC + BM
f = params.nos.f;
yOM_nos  = (dBM/f)*(1/dOM);
yN2O_nos = (dBM/f)*((1-f)/2);   
eN2_nos  = (dBM/f)*((1-f)/2); 
eNH4_nos = (dBM/f)*((nOM/dOM)-((nBM*f)/dBM));
eDIC_nos = (dBM/f)*((cOM/dOM)-((cBM*f)/dBM));
ePO4_nos = (dBM/f)*((pOM/dOM)-((pBM*f)/dBM));
% Facultative aerboic
f = params.nos.f_fac;
yOM_nos_fac  = (dBM/f)*(1/dOM);
yO2_nos_fac  = (dBM/f)*((1-f)/4);
eNH4_nos_fac = (dBM/f)*((nOM/dOM)-((nBM*f)/dBM)); 
eDIC_nos_fac = (dBM/f)*((cOM/dOM)-((cBM*f)/dBM));
ePO4_nos_fac = (dBM/f)*((pOM/dOM)-((pBM*f)/dBM));

% ---------------------------- %
% NO3 reduction to N2O (DEN 4) %
% ---------------------------- %
% Reactants: OM + NO3
% Products:  NH4 + N2O + DIC + BM
f = params.nai.f;
yOM_nai  = (dBM/f)*(1/dOM);
yNO3_nai = (dBM/f)*((1-f)/4);
eN2O_nai = (dBM/f)*((1-f)/8);
eNH4_nai = (dBM/f)*((nOM/dOM)-((nBM*f)/dBM));
eDIC_nai = (dBM/f)*((cOM/dOM)-((cBM*f)/dBM));
ePO4_nai = (dBM/f)*((pOM/dOM)-((pBM*f)/dBM));
% Facultative aerboic
f = params.nai.f_fac;
yOM_nai_fac  = (dBM/f)*(1/dOM);
yO2_nai_fac  = (dBM/f)*((1-f)/4);
eNH4_nai_fac = (dBM/f)*((nOM/dOM)-((nBM*f)/dBM)); 
eDIC_nai_fac = (dBM/f)*((cOM/dOM)-((cBM*f)/dBM));
ePO4_nai_fac = (dBM/f)*((pOM/dOM)-((pBM*f)/dBM));

% --------------------------- %
% NO2 reduction to N2 (DEN 5) %
% --------------------------- %
% Reactants: OM + NO2
% Products:  NH4 + N2 + DIC + BM
f = params.nio.f;
yOM_nio  = (dBM/f)*(1/dOM);
yNO2_nio = (dBM/f)*((1-f)/3);
eN2_nio  = (dBM/f)*((1-f)/6);
eNH4_nio = (dBM/f)*((nOM/dOM)-((nBM*f)/dBM));
eDIC_nio = (dBM/f)*((cOM/dOM)-((cBM*f)/dBM));
ePO4_nio = (dBM/f)*((pOM/dOM)-((pBM*f)/dBM));
% Facultative aerboic
f = params.nio.f_fac;
yOM_nio_fac  = (dBM/f)*(1/dOM);
yO2_nio_fac  = (dBM/f)*((1-f)/4);
eNH4_nio_fac = (dBM/f)*((nOM/dOM)-((nBM*f)/dBM)); 
eDIC_nio_fac = (dBM/f)*((cOM/dOM)-((cBM*f)/dBM));
ePO4_nio_fac = (dBM/f)*((pOM/dOM)-((pBM*f)/dBM));

% --------------------------- %
% NO3 reduction to N2 (DEN 6) %
% --------------------------- %
% Reactants: OM + NO3
% Products:  NH4 + N2 + DIC + BM
f = params.nao.f;
yOM_nao  = (dBM/f)*(1/dOM);
yNO3_nao = (dBM/f)*((1-f)/5);
eN2_nao  = (dBM/f)*((1-f)/10);
eNH4_nao = (dBM/f)*((nOM/dOM)-((nBM*f)/dBM));
eDIC_nao = (dBM/f)*((cOM/dOM)-((cBM*f)/dBM));
ePO4_nao = (dBM/f)*((pOM/dOM)-((pBM*f)/dBM));
% Facultative aerboic
f = params.nao.f_fac;
yOM_nao_fac  = (dBM/f)*(1/dOM);
yO2_nao_fac  = (dBM/f)*((1-f)/4);
eNH4_nao_fac = (dBM/f)*((nOM/dOM)-((nBM*f)/dBM)); 
eDIC_nao_fac = (dBM/f)*((cOM/dOM)-((cBM*f)/dBM));
ePO4_nao_fac = (dBM/f)*((pOM/dOM)-((pBM*f)/dBM));

% ----------------------------- %
% NO3 reduction to NH4 (DNRA 1) %
% ----------------------------- %
% Reactants: OM + NO3
% Products:  NH4 + BM + DIC
if opt.DNRA1
	f = params.dnra1.f;
	yOM_dnra1  = (dBM/f)*(1/dOM);
	yNO3_dnra1 = (dBM/f)*((1-f)/8); 
	eNH4_dnra1 = (dBM/f)*((nOM/dOM)-((nBM*f)/dBM)); % excrete extra NH4 in OM 
	eNH4_dnra1 = eNH4_dnra1 + (dBM/f)*((1-f)/8);    % excrete NH4 via reduction
	eDIC_dnra1 = (dBM/f)*((cOM/dOM)-((cBM*f)/dBM)); 
	ePO4_dnra1 = (dBM/f)*((pOM/dOM)-((pBM*f)/dBM));
	% Facultative aerboic
	f = params.dnra1.f_fac;
	yOM_dnra1_fac  = (dBM/f)*(1/dOM);
	yO2_dnra1_fac  = (dBM/f)*((1-f)/4);
	eNH4_dnra1_fac = (dBM/f)*((nOM/dOM)-((nBM*f)/dBM)); 
	eDIC_dnra1_fac = (dBM/f)*((cOM/dOM)-((cBM*f)/dBM));
	ePO4_dnra1_fac = (dBM/f)*((pOM/dOM)-((pBM*f)/dBM));
end

% ----------------------------- %
% NO2 reduction to NH4 (DNRA 2) %
% ----------------------------- %
% Reactants: OM + NO2
% Products:  NH4 + BM + DIC
if opt.DNRA2
	f = params.dnra2.f;
	yOM_dnra2  = (dBM/f)*(1/dOM);
	yNO2_dnra2 = (dBM/f)*((1-f)/6);
	eNH4_dnra2 = (dBM/f)*((nOM/dOM)-((nBM*f)/dBM)); % excrete extra NH4 in OM
	eNH4_dnra2 = eNH4_dnra2 + (dBM/f)*((1-f)/6);    % excrete NH4 via reduction 
	eDIC_dnra2 = (dBM/f)*((cOM/dOM)-((cBM*f)/dBM));
	ePO4_dnra2 = (dBM/f)*((pOM/dOM)-((pBM*f)/dBM));
	% Facultative aerboic
	f = params.dnra2.f_fac;
	yOM_dnra2_fac  = (dBM/f)*(1/dOM);
	yO2_dnra2_fac  = (dBM/f)*((1-f)/4);
	eNH4_dnra2_fac = (dBM/f)*((nOM/dOM)-((nBM*f)/dBM)); 
	eDIC_dnra2_fac = (dBM/f)*((cOM/dOM)-((cBM*f)/dBM));
	ePO4_dnra2_fac = (dBM/f)*((pOM/dOM)-((pBM*f)/dBM));
end

% ---------------------------------------
% Chemoautotrophic bacteria
% ---------------------------------------

% --------------------
% NH4 oxidation to NO2 
% --------------------
% Reactants: NH4 + O2 + DIC
% Products:  NO2 + BM
f = params.aoa.f;
yNH4_aoa = (dBM/f)*((1/6) + ((nBM*f)/dBM));
yO2_aoa  = (dBM/f)*((1-f)/4);
yDIC_aoa = (dBM/f)*((cBM*f)/dBM); 
yPO4_aoa = (dBM/f)*((pBM*f)/dBM);
eNO2_aoa = (dBM/f)*(1/6);

% --------------------
% NO2 oxidation to NO3
% --------------------
% Reactants: NO2 + NH4 + O2 + DIC
% Products:  NO3 + BM
f = params.nob.f;
yNO2_nob = (dBM/f)*(1/2);
yO2_nob  = (dBM/f)*((1-f)/4);
yDIC_nob = (dBM/f)*((cBM*f)/dBM); 
yNH4_nob = (dBM/f)*((nBM*f)/dBM); 
yPO4_nob = (dBM/f)*((pBM*f)/dBM);
eNO3_nob = (dBM/f)*(1/2);

% -------------------
% Anammox to NO3 + N2
% -------------------
% Reactants: NH4 + NO2 + DIC
% Products:  NO3 + N2 + BM
f = params.aox.f; 
x = params.aox.x;
yNH4_aox = (dBM/f)*((x/3)+((1-x)/8)+(f*nBM/dBM));
yNO2_aox = (dBM/f)*((1-f)/3);
yDIC_aox = (dBM/f)*((cBM*f)/dBM);
yPO4_aox = (dBM/f)*((pBM*f)/dBM);
eNO3_aox = (dBM/f)*((1-x)/8);
eN2_aox  = (dBM/f)*((x+1-f)/6);

% Heterotroph yields
% AER
% ... reactants
params.aer.y_red_aer = yOM_aer; % aerobic BM yield via DOC (mol DOC / mol BC)
params.aer.y_oxy_aer = yO2_aer; % aerobic BM yield via O2  (mol O2  / mol BC)
params.aer.y_red_ana = yOM_aer; % anaerobic BM yield via DOC (unused)
params.aer.y_oxy_ana = yO2_aer; % anaerboic BM yield via O2 (unused)
% ... products
params.aer.e_nh4_aer = eNH4_aer; % aerobic NH4 excretion (mol NH4 / mol BC)
params.aer.e_dic_aer = eDIC_aer; % aerobic DIC excretion (mol DIC / mol BC)
params.aer.e_po4_aer = ePO4_aer; % aerobic PO4 excretion (mol PO4 / mol BC)
params.aer.e_nh4_ana = eNH4_aer; % anaerboic NH4 excretion (unused)
params.aer.e_dic_ana = eDIC_aer; % anaerobic DIC excretion (unused)
params.aer.e_po4_ana = ePO4_aer; % anaerobic PO4 excretion (unused)

% NAR
% ... reactants
params.nar.y_red_aer = yOM_nar_fac;  % aerobic BM yield via DOC (facultative)
params.nar.y_oxy_aer = yO2_nar_fac;  % aerobic BM yield via O2  (facultative)
params.nar.y_red_ana = yOM_nar;      % anaerobic BM yield via DOC (mol DOC / mol BC)
params.nar.y_oxy_ana = yNO3_nar;     % anaerobic BM yield via NO3 (mol NO3 / mol BC)
% ... products
params.nar.e_nh4_aer = eNH4_nar_fac; % aerobic NH4 excretion (facultative)
params.nar.e_dic_aer = eDIC_nar_fac; % aerobic DIC excretion (facultative)
params.nar.e_po4_aer = ePO4_nar_fac; % aerobic PO4 excretion (facultative) 
params.nar.e_nh4_ana = eNH4_nar;     % anaerobic NH4 excretion (mol NH4 / mol BC)
params.nar.e_no2_ana = eNO2_nar;     % anaerboic NO2 excretion (mol NO2 / mol BC)
params.nar.e_dic_ana = eDIC_nar;     % anaerobic DIC excretion (mol DIC / mol BC)
params.nar.e_po4_ana = ePO4_nar;     % anaerobic PO4 excretion (mol PO4 / mol BC)

% NIR
% ... reactants
params.nir.y_red_aer = yOM_nir_fac;  % aerobic BM yield via DOC (facultative)
params.nir.y_oxy_aer = yO2_nir_fac;  % aerobic BM yield via O2  (facultative)
params.nir.y_red_ana = yOM_nir;      % anaerobic BM yield via DOC (mol DOC / mol BC)
params.nir.y_oxy_ana = yNO2_nir;     % anaerobic BM yield via NO2 (mol NO2 / mol BC)
% ... products
params.nir.e_nh4_aer = eNH4_nir_fac; % aerobic NH4 excretion (facultative)
params.nir.e_dic_aer = eDIC_nir_fac; % aerobic DIC excretion (facultative)
params.nir.e_po4_aer = ePO4_nir_fac; % aerobic PO4 excretion (facultative) 
params.nir.e_nh4_ana = eNH4_nir;     % anaerobic NH4 excretion (mol NH4 / mol BC)
params.nir.e_n2o_ana = eN2O_nir;     % anaerboic N2O excretion (mol N2O / mol BC)
params.nir.e_dic_ana = eDIC_nir;     % anaerobic DIC excretion (mol DIC / mol BC)
params.nir.e_po4_ana = ePO4_nir;     % anaerobic PO4 excretion (mol PO4 / mol BC)

% NOS
% ... reactants
params.nos.y_red_aer = yOM_nos_fac;  % aerobic BM yield via DOC (facultative)
params.nos.y_oxy_aer = yO2_nos_fac;  % aerobic BM yield via O2  (facultative)
params.nos.y_red_ana = yOM_nos;      % anaerobic BM yield via DOC (mol DOC / mol BC)
params.nos.y_oxy_ana = yN2O_nos;     % anaerobic BM yield via N2O (mol N2O / mol BC)
% ... products
params.nos.e_nh4_aer = eNH4_nos_fac; % aerobic NH4 excretion (facultative)
params.nos.e_dic_aer = eDIC_nos_fac; % aerobic DIC excretion (facultative)
params.nos.e_po4_aer = ePO4_nos_fac; % aerobic PO4 excretion (facultative) 
params.nos.e_nh4_ana = eNH4_nos;     % anaerobic NH4 excretion (mol NH4 / mol BC)
params.nos.e_n2_ana  = eN2_nos;      % anaerobic N2 excretion  (mol N2 / mol BC)
params.nos.e_dic_ana = eDIC_nos;     % anaerobic DIC excretion (mol DIC / mol BC)
params.nos.e_po4_ana = ePO4_nos;     % anaerobic PO4 excretion (mol PO4 / mol BC)

% NAI
% ... reactants
params.nai.y_red_aer = yOM_nai_fac;  % aerobic BM yield via DOC (facultative)
params.nai.y_oxy_aer = yO2_nai_fac;  % aerobic BM yield via O2  (facultative)
params.nai.y_red_ana = yOM_nai;      % anaerobic BM yield via DOC (mol DOC / mol BC)
params.nai.y_oxy_ana = yNO3_nai;     % anaerboic BM yield via NO3 (mol NO3 / mol BC)
% ... products
params.nai.e_nh4_aer = eNH4_nai_fac; % aerobic NH4 excretion (facultative)
params.nai.e_dic_aer = eDIC_nai_fac; % aerobic DIC excretion (facultative)
params.nai.e_po4_aer = ePO4_nai_fac; % aerobic PO4 excretion (facultative) 
params.nai.e_nh4_ana = eNH4_nai;     % anaerobic NH4 excretion (mol NH4 / mol BC)
params.nai.e_n2o_ana = eN2O_nai;     % anaerobic N2O excretion (mol N2O / mol BC) 
params.nai.e_dic_ana = eDIC_nai;     % anaerobic DIC excretion (mol DIC / mol BC)
params.nai.e_po4_ana = ePO4_nai;     % anaerobic PO4 excretion (mol PO4 / mol BC)

% NAO
% ... reactants
params.nao.y_red_aer = yOM_nao_fac;  % aerobic BM yield via DOC (facultative)
params.nao.y_oxy_aer = yO2_nao_fac;  % aerobic BM yield via O2  (facultative)
params.nao.y_red_ana = yOM_nao;      % anaerobic BM yield via DOC (mol DOC / mol BC)
params.nao.y_oxy_ana = yNO3_nao;     % anaerboic BM yield via NO3 (mol NO3 / mol BC)
% ... products
params.nao.e_nh4_aer = eNH4_nao_fac; % aerobic NH4 excretion (facultative)
params.nao.e_dic_aer = eDIC_nao_fac; % aerobic DIC excretion (facultative)
params.nao.e_po4_aer = ePO4_nao_fac; % aerobic PO4 excretion (facultative) 
params.nao.e_nh4_ana = eNH4_nao;     % anaerboic NH4 excretion (mol NH4 / mol BC)
params.nao.e_n2_ana  = eN2_nao;      % anaerobic N2 excretion  (mol N2  / mol BC)
params.nao.e_dic_ana = eDIC_nao;     % anaerobic DIC excretion (mol DIC / mol BC)
params.nao.e_po4_ana = ePO4_nao;     % anaerobic PO4 excretion (mol PO4 / mol BC)

% NIO
% ... reactants
params.nio.y_red_aer = yOM_nio_fac;  % aerobic BM yield via DOC (facultative)
params.nio.y_oxy_aer = yO2_nio_fac;  % aerobic BM yield via O2  (facultative)
params.nio.y_red_ana = yOM_nio;      % anaerobic BM yield via DOC (mol DOC / mol BC)
params.nio.y_oxy_ana = yNO2_nio;     % anaerobic BM yield via NO2 (mol NO2 / mol BC)
% ... products
params.nio.e_nh4_aer = eNH4_nio_fac; % aerobic NH4 excretion (facultative)
params.nio.e_dic_aer = eDIC_nio_fac; % aerobic DIC excretion (facultative)
params.nio.e_po4_aer = ePO4_nio_fac; % aerobic PO4 excretion (facultative) 
params.nio.e_nh4_ana = eNH4_nio;     % anaerobic NH4 excretion (mol NH4 / mol BC)
params.nio.e_n2_ana  = eN2_nio;      % anaerobic N2 excretion  (mol N2  / mol BC)
params.nio.e_dic_ana = eDIC_nio;     % anaerobic DIC excretion (mol DIC / mol BC)
params.nio.e_po4_ana = ePO4_nio;     % anaerobic PO4 excretion (mol PO4 / mol BC)

% DNRA 1
if opt.DNRA1
	% ... reactants
	params.dnra1.y_red_aer = yOM_dnra1_fac;  % aerobic BM yield via DOC (facultative)
	params.dnra1.y_oxy_aer = yO2_dnra1_fac;  % aerobic BM yield via O2  (facultative)
	params.dnra1.y_red_ana = yOM_dnra1;      % anaerobic BM yield via DOC (mol DOC / mol BC)
	params.dnra1.y_oxy_ana = yNO3_dnra1;     % anaerobic BM yield via NO2 (mol NO3 / mol BC)
	% ... products
	params.dnra1.e_nh4_aer = eNH4_dnra1_fac; % aerobic NH4 excretion (facultative)
	params.dnra1.e_dic_aer = eDIC_dnra1_fac; % aerobic DIC excretion (facultative)
	params.dnra1.e_po4_aer = ePO4_dnra1_fac; % aerobic PO4 excretion (facultative)
	params.dnra1.e_nh4_ana = eNH4_dnra1;     % anaerobic NH4 excretion (mol NH4 / mol BC)
	params.dnra1.e_dic_ana = eDIC_dnra1;     % anaerobic DIC excretion (mol DIC / mol BC)
	params.dnra1.e_po4_ana = ePO4_dnra1;     % anaerobic PO4 excretion (mol PO4 / mol BC)
end

% DNRA 2
if opt.DNRA2
	% ... reactants
	params.dnra2.y_red_aer = yOM_dnra2_fac;  % aerobic BM yield via DOC (facultative)
	params.dnra2.y_oxy_aer = yO2_dnra2_fac;  % aerobic BM yield via O2  (facultative)
	params.dnra2.y_red_ana = yOM_dnra2;      % anaerobic BM yield via DOC (mol DOC / mol BC)
	params.dnra2.y_oxy_ana = yNO2_dnra2;     % anaerobic BM yield via NO2 (mol NO3 / mol BC)
	% ... products
	params.dnra2.e_nh4_aer = eNH4_dnra2_fac; % aerobic NH4 excretion (facultative)
	params.dnra2.e_dic_aer = eDIC_dnra2_fac; % aerobic DIC excretion (facultative)
	params.dnra2.e_po4_aer = ePO4_dnra2_fac; % aerobic PO4 excretion (facultative)
	params.dnra2.e_nh4_ana = eNH4_dnra2;     % anaerobic NH4 excretion (mol NH4 / mol BC)
	params.dnra2.e_dic_ana = eDIC_dnra2;     % anaerobic DIC excretion (mol DIC / mol BC)
	params.dnra2.e_po4_ana = ePO4_dnra2;     % anaerobic PO4 excretion (mol PO4 / mol BC)
end

% Chemoautotroph yields
% AOA
% ... reactants
params.aoa.y_red = yNH4_aoa; % aerobic BM yield via NH4 (mol NH4 / mol BC)
params.aoa.y_oxy = yO2_aoa;  % aerobic BM yield via O2  (mol O2  / mol BC)
params.aoa.y_dic = yDIC_aoa; % aerobic BM yield via DIC (mol DIC / mol BC)
params.aoa.y_po4 = yPO4_aoa; % aerobic BM yield via PO4 (mol PO4 / mol BC)
% ... products
params.aoa.e_no2 = eNO2_aoa; % aerobic NO2 excretion (mol NO2 / mol BC) 

% NOB
% ... reactants
params.nob.y_red = yNO2_nob; % aerobic BM yield via NO2 (mol NO2 / mol BC)
params.nob.y_oxy = yO2_nob;  % aerobic BM yield via O2  (mol O2  / mol BC)
params.nob.y_nh4 = yNH4_nob; % aerobic BM yield via NH4 (mol NH4 / mol BC)
params.nob.y_dic = yDIC_nob; % aerobic BM yield via DIC (mol DIC / mol BC)
params.nob.y_po4 = yPO4_nob; % aerobic BM yield via PO4 (mol PO4 / mol BC)
% ... products
params.nob.e_no3 = eNO3_nob; % aerobic NO3 excretion (mol NO3 / mol BC)

% AOX
% ... reactants
params.aox.y_red = yNH4_aox; % anaerobic BM yield via NO2 (mol NO2 / mol BC)
params.aox.y_oxy = yNO2_aox; % anaerobic BM yield via NH4 (mol NH4 / mol BC)
params.aox.y_dic = yDIC_aox; % anaerobic BM yield via DIC (mol DIC / mol BC)
params.aox.y_po4 = yPO4_aox; % anaerobic BM yield via PO4 (mol PO4 / mol BC)
% ... products
params.aox.e_no3 = eNO3_aox; % anaerobic NO3 excretion (mol NO3 / mol BC)
params.aox.e_n2  = eN2_aox;  % anaerobic N2 excretion  (mol N2  / mol BC)

% ---------------------------------------
% Report values 
% ---------------------------------------
disp(' ')
disp('---------------')
disp('Chemoautotrophs')
disp('---------------')

disp(' ');disp('--- AOA ---');
eq = [num2str(yNH4_aoa),'(NH4) + ',...
      num2str(yDIC_aoa),'(DIC) + ',...
      num2str(yO2_aoa), '(O2) = ', ...
      '1(B) + ',...
      num2str(eNO2_aoa),'(NO2)'];
disp(eq);

disp(' ');disp('--- NOB ---');
eq = [num2str(yNO2_nob),'(NO2) + ',...
      num2str(yDIC_nob),'(DIC) + ',...
      num2str(yNH4_nob),'(NH4) + ',...
      num2str(yO2_nob), '(O2) = ',...
      '1(B) + ',...
      num2str(eNO3_nob),'(NO3)'];
disp(eq);

disp(' ');disp('--- AOX ---');
eq = [num2str(yNH4_aox),'(NH4) + ',...
      num2str(yDIC_aox),'(DIC) + ',...
      num2str(yNO2_aox),'(NO2) = ',...
      '1(B) + ',...
      num2str(eNO3_aox),'(NO3) + ',...
      num2str(eN2_aox), '(N2)']; 
disp(eq);

disp(' ');
disp('-------------')
disp('Heterotrophs')
disp('-------------')

% Aerobic metabolism
disp(' ');disp('--- AER ---');
eq = [num2str(yOM_aer), '(OM) + ',...
      num2str(yO2_aer), '(O2) = ',...
      '1(B) + ',...
      num2str(eDIC_aer),'(DIC) + ',...
      num2str(eNH4_aer),'(NH4)'];
disp(eq);

% Anaerobic metabolism
if opt.FACULTATIVE_MICROBES
	disp(' ');disp('--- NAR (FAC) ---');
	eq = [num2str(yOM_nar_fac), '(OM) + ',...
		  num2str(yO2_nar_fac), '(O2) = ',...
		  '1(B) + ',...
		  num2str(eDIC_nar_fac),'(DIC) + ',...
		  num2str(eNH4_nar_fac),'(NH4)'];
	disp(eq);
	disp(' ');disp('--- NIR (FAC) ---');
	eq = [num2str(yOM_nir_fac), '(OM) + ',...
		  num2str(yO2_nir_fac), '(O2) = ',...
		  '1(B) + ',...
		  num2str(eDIC_nir_fac),'(DIC) + ',...
		  num2str(eNH4_nir_fac),'(NH4)'];
	disp(eq);
	disp(' ');disp('--- NOS (FAC) ---');
	eq = [num2str(yOM_nos_fac), '(OM) + ',...
		  num2str(yO2_nos_fac), '(O2) = ',...
		  '1(B) + ',...
		  num2str(eDIC_nos_fac),'(DIC) + ',...
		  num2str(eNH4_nos_fac),'(NH4)'];
	disp(eq);
	disp(' ');disp('--- NAI (FAC) ---');
	eq = [num2str(yOM_nai_fac), '(OM) + ',...
		  num2str(yO2_nai_fac), '(O2) = ',...
		  '1(B) + ',...
		  num2str(eDIC_nai_fac),'(DIC) + ',...
		  num2str(eNH4_nai_fac),'(NH4)'];
	disp(eq);
	disp(' ');disp('--- NIO (FAC) ---');
	eq = [num2str(yOM_nio_fac), '(OM) + ',...
		  num2str(yO2_nio_fac), '(O2) = ',...
		  '1(B) + ',...
		  num2str(eDIC_nio_fac),'(DIC) + ',...
		  num2str(eNH4_nio_fac),'(NH4)'];
	disp(eq);
	disp(' ');disp('--- NAO (FAC) ---');
	eq = [num2str(yOM_nao_fac), '(OM) + ',...
		  num2str(yO2_nao_fac), '(O2) = ',...
		  '1(B) + ',...
		  num2str(eDIC_nao_fac),'(DIC) + ',...
		  num2str(eNH4_nao_fac),'(NH4)'];
	disp(eq);
	if opt.DNRA1
		disp(' ');disp('--- DNRA1 (FAC) ---');
		eq = [num2str(yOM_dnra1_fac), '(OM) + ',...
			  num2str(yO2_dnra1_fac), '(O2) = ',...
			  '1(B) + ',...
			  num2str(eDIC_dnra1_fac),'(DIC) + ',...
			  num2str(eNH4_dnra1_fac),'(NH4)'];
		disp(eq);
	end
	if opt.DNRA2
		disp(' ');disp('--- DNRA2 (FAC) ---');
		eq = [num2str(yOM_dnra2_fac), '(OM) + ',...
			  num2str(yO2_dnra2_fac), '(O2) = ',...
			  '1(B) + ',...
			  num2str(eDIC_dnra2_fac),'(DIC) + ',...
			  num2str(eNH4_dnra2_fac),'(NH4)'];
		disp(eq);
	end
end
disp(' ');disp('--- NAR (ANA) ---');
eq = [num2str(yOM_nar), '(OM) + ',...
	  num2str(yNO3_nar),'(NO3) = ',...
	  '1(B) + ',...
	  num2str(eDIC_nar),'(DIC) + ',...
	  num2str(eNH4_nar),'(NH4) + ',...
	  num2str(eNO2_nar),'(NO2)'];
disp(eq);
disp(' ');disp('--- NIR (ANA)---');
eq = [num2str(yOM_nir), '(OM) + ',...
	  num2str(yNO2_nir),'(NO2) = ',...
	  '1(B) + ',...
	  num2str(eDIC_nir),'(DIC) + ',...
	  num2str(eNH4_nir),'(NH4) + ',...
	  num2str(eN2O_nir),'(N2O)'];
disp(eq);
disp(' ');disp('--- NOS (ANA) ---');
eq = [num2str(yOM_nos), '(OM) + ',...
	  num2str(yN2O_nos),'(N2O) = ',...
	  '1(B) + ',...
	  num2str(eDIC_nos),'(DIC) + ',...
	  num2str(eNH4_nos),'(NH4) + ',...
	  num2str(eN2_nos), '(N2)'];
disp(eq);
disp(' ');disp('--- NAI (ANA) ---');
eq = [num2str(yOM_nai), '(OM) + ',...
	  num2str(yNO3_nai),'(NO3) = ',...
	  '1(B) + ',...
	  num2str(eDIC_nai),'(DIC) + ',...
	  num2str(eNH4_nai),'(NH4) + ',...
	  num2str(eN2O_nai),'(N2O)'];
disp(eq);
disp(' ');disp('--- NIO (ANA) ---');
eq = [num2str(yOM_nio), '(OM) + ',...
	  num2str(yNO2_nio),'(NO2) = ',...
	  '1(B) + ',...
	  num2str(eDIC_nio),'(DIC) + ',...
	  num2str(eNH4_nio),'(NH4) + ',...
	  num2str(eN2_nio), '(N2)'];
disp(eq);
disp(' ');disp('--- NAO (ANA) ---');
eq = [num2str(yOM_nao), '(OM) + ',...
	  num2str(yNO3_nao),'(NO3) = ',...
	  '1(B) + ',...
	  num2str(eDIC_nao),'(DIC) + ',...
	  num2str(eNH4_nao),'(NH4) + ',...
	  num2str(eN2_nao), '(N2)'];
disp(eq);
if opt.DNRA1
	disp(' ');disp('--- DNRA1 (ANA) ---');
	eq = [num2str(yOM_dnra1), '(OM) + ',...
		  num2str(yNO3_dnra1),'(NO3) = ',...
		  '1(B) + ',...
		  num2str(eDIC_dnra1),'(DIC) + ',...
		  num2str(eNH4_dnra1),'(NH4)'];
	disp(eq);
end
if opt.DNRA2
	disp(' ');disp('--- DNRA2 (ANA) ---');
	eq = [num2str(yOM_dnra2), '(OM) + ',...
		  num2str(yNO2_dnra2),'(NO2) = ',...
		  '1(B) + ',...
		  num2str(eDIC_dnra2),'(DIC) + ',...
		  num2str(eNH4_dnra2),'(NH4)'];
	disp(eq);
	disp(' ');
end
