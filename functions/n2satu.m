function n2sat = n2satu(sst,sss,params);
%---------------------------------------------------------------------------
%
%   Computes dinitrogen saturation concentration at 1 atm total pressure
%   in mmol/m^3 given the temperature (t, in deg C) and the salinity (s  in permil) 
%
%   FROM Sarmiento and Gruber 2006
%   --> Solubility : table 3.2.2
%   --> apply correction for moist air following Panel 3.2.1
%
%   N2SAT IS DEFINED BETWEEN T(freezing) <= T <= 40(deg C) AND
%   Simon Yang, Feb 2018
%---------------------------------------------------------------------------

a_1 = -59.6274;
a_2 = 85.7761;
a_3 = 24.3696;
b_1 = -0.051580;
b_2 = 0.026329;
b_3 = -0.0037252;
vbarn2o = 22.4136;  % molar volume of N2 at standard temperature and pressure

% Avoid bad SSS
S_LOC = max(5,sss);

% Convert SST to Kelvin
TS    = sst + 273.16;

% Calc
bunsen = exp(a_1 + a_2 * (100.0/TS) + a_3 * log(TS/100.0) + ...
             S_LOC * (b_1 + b_2 * (TS/100.0) + b_3 * ...
             (TS/100.0).^2));
PH2O_to_P = exp(24.4543-67.4509*(100.0/TS)-4.8489*log(TS/100.0) ...
              -0.000544*S_LOC);
n2sat = ( bunsen / vbarn2o) * (1 - PH2O_to_P);
n2sat = n2sat .* 1e6;


