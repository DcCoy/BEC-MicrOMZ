function [flux] = airsea_flux(tr,phy,frc,params,grid,opt);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Air-sea flux of dissolved gases 
% Adopted from Wanninkhof et al., 1992 (eq 3)
% Parameters via Sarmiento & Gruber, 2006
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get wind-stress from forcing
ws_used = wind_stress(frc,params);    

% Get piston velocity (Kw) conversion factor (from (cm/hr) to (m/s))
xkw = params.bec.a * ws_used * ws_used; 

% Get schmidt #s for dissolved gasses using SST
schmidt_o2  = cschmidt_o2(phy.temp(1));
schmidt_n2o = cschmidt_n2o(phy.temp(1));
schmidt_n2  = cschmidt_n2(phy.temp(1));

% Get atmospheric saturation concentrations at surface    
sat_o2_1atm  = o2satu(phy.temp(1),phy.salt(1));
sat_n2o_1atm = n2osatu(phy.temp(1),phy.salt(1),params);
if opt.N2_FULL
    sat_n2_1atm = n2satu(phy.temp(1),phy.salt(1),params);
else
    % Only track excess N2
    sat_n2_1atm = 0;
end

% Get piston velocities
pv_n2o = xkw * sqrt(660 / schmidt_n2o);  
pv_n2  = xkw * sqrt(660 / schmidt_n2);
pv_o2  = xkw * sqrt(660 / schmidt_o2);

% Multiply saturation by molar ratios of air
sat_n2o_1atm = sat_n2o_1atm .* params.bec.xn2o;
sat_n2_1atm  = sat_n2_1atm .* params.bec.xn2;

% Collect surface saturation differences
n2o_diff = sat_n2o_1atm - tr.n2o(1);
n2_diff  = sat_n2_1atm  - tr.n2(1);
o2_diff  = sat_o2_1atm  - tr.o2(1);

% Get fluxes in mmol/m2s (loss or gain term at top cell)
flux.n2o = (pv_n2o * n2o_diff); 
flux.n2  = (pv_n2  * n2_diff);  
flux.o2  = (pv_o2  * o2_diff); 

% Reduce flux if it leads to overshoots
% Applied due to large timesteps
if abs(flux.n2o/grid.Hz(1)).*opt.dt > abs(n2o_diff)
    flux.n2o = abs(n2o_diff./opt.dt).*grid.Hz(1).*sign(flux.n2o);    
end
if abs(flux.n2/grid.Hz(1)).*opt.dt > abs(n2_diff)
    flux.n2 = abs(n2_diff./opt.dt).*grid.Hz(1).*sign(flux.n2);    
end
if abs(flux.o2/grid.Hz(1)).*opt.dt > abs(o2_diff)
    flux.o2 = abs(o2_diff./opt.dt).*grid.Hz(1).*sign(flux.o2);    
end




