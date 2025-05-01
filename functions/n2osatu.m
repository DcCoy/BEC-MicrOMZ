function n2osat = n2osatu(temp,salt,params);
%---------------------------------------------------------------------------
%
%   Computes n2o saturation concentration at 1 atm n2o pressure
%   in mmol/m^3 given the temperature (t, in deg C) and the salinity (s, in permil) .
%
%   FROM WEISS AND PRICE (1980), Marine chemistry.
%
%---------------------------------------------------------------------------

% Convert temperature to kelvin, avoid bad salinity
T = temp + 273.15;
S = max(salt,5);

% Coefficients for the fit of solubility
% Page 78, chapter 3, S&G
a1 = -165.8806;
a2 = 222.8743;
a3 = 92.0792;
a4 = -1.48425;
b1 = -0.056235;
b2 = 0.031619;
b3 = -0.0048472;

% Equation in ROMS 
F = exp(a1 + a2.*(100./T) + a3.*log(T./100) + a4.*(T./100).^2 + S.*(b1 + b2.*(T./100) + b3.*(T./100).^2));

% Get pH2O/P
pH2O = exp(24.4543 - 67.4509.*(100./T) - 4.8489.*log(T./100) - 0.000544.*S);

% Get pAmoist (excluding xn2o)
pAmoist = (1-(pH2O./1));

% Get Solubility parameter
Sa = (F./(1-pH2O))*1000*1000;

% Get N2O
n2osat = Sa.*pAmoist;
n2osat(n2osat==-inf) = [];

