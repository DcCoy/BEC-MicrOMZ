function restoring = restoring(opt,inputs,grid,tr);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Specifies restoring source (linear relaxation) for selected dissolved tracers
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% tauZvar
if opt.TAUZVAR 
    tauh = ppval(inputs.restore.tauh_pre,grid.z_r);
else
    tauh = phy.tauh;
end

% PO4
if ~opt.rest_po4
    restoring.po4 = zeros(length(tr.po4),1);
else
    restoring.po4 = (inputs.restore.po4_cout-tr.po4)./tauh;
end

% NO3
if ~opt.rest_no3
    restoring.no3 = zeros(length(tr.no3),1);
else
    restoring.no3 = (inputs.restore.no3_cout-tr.no3)./tauh; 
end

% O2
if ~opt.rest_o2
    restoring.o2 = zeros(length(tr.o2),1);
else
%    if bgc.forceanoxic == 1
%        cout(find(grid.z_r==bgc.forceanoxic_bounds(2)):find(grid.z_r==bgc.forceanoxic_bounds(1)))=0;
%    end
    restoring.o2 = (inputs.restore.o2_cout-tr.o2)./tauh;
end

% N2O
if ~opt.rest_n2o
    restoring.n2o = zeros(length(tr.n2o),1);
else
    restoring.n2o = (inputs.restore.n2o_cout-tr.n2o)./tauh;
end

% NO2
if ~opt.rest_no2 
    restoring.no2 = zeros(length(tr.no2),1);
else
    restoring.no2 = (inputs.restore.no2_cout-tr.no2)./tauh;
end

% Fe
if ~opt.rest_fe
    restoring.fe = zeros(length(tr.fe),1);
else
    restoring.fe = (inputs.restore.fe_cout-tr.fe)./tauh;
end

% SiO3
if ~opt.rest_sio3
    restoring.sio3 = zeros(length(tr.sio3),1);
else
    restoring.sio3 = (inputs.restore.sio3_cout-tr.sio3)./tauh;
end
