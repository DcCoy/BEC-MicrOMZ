% Script to build theoretical 'coammox' type
clear all
load('../output/MicrOMZ_obligate_baseline_Km.mat');

% Calculate 'Y' (excreted NO2 from AOA / required NO2 from NOB);
Y = params.aoa.e_no2 ./ params.nob.y_red;

% Get required reactants
nh4 = (params.aoa.y_red + params.nob.y_nh4*Y); 
dic = (params.aoa.y_dic + params.nob.y_dic*Y);
o2  = (params.aoa.y_oxy + params.nob.y_oxy*Y);

% Get produced products
BC  = (1 + Y);
no3 = (params.aoa.e_no2);

% Display total equation
disp(' ');
disp(' --- ');
disp('TOTAL');
disp([ ...
	num2str(round(nh4,4)),'(NH4) + ',...
	num2str(round(dic,4)),'(DIC) + ',...
	num2str(round(o2,4)),'(O2) = ',...
	num2str(round(BC,4)),'(B_c) + ',...
	num2str(round(no3,4)),'(NO3)']);
disp(' ');

% Display normalized
nh4 = nh4/BC;
o2  = o2/BC;
dic = dic/BC;
no3 = no3/BC;
BC  = BC/BC;
disp(' ');
disp(' ------------------ ');
disp('Normalized to 1(B_c)');
disp([ ...
	num2str(round(nh4,4)),'(NH4) + ',...
	num2str(round(dic,4)),'(DIC) + ',...
	num2str(round(o2,4)),'(O2) = ',...
	num2str(round(BC,4)),'(B_c) + ',...
	num2str(round(no3,4)),'(NO3)']);
disp(' ');
return
