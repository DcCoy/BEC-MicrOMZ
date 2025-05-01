function o2sat = o2satu(sst,sss);
 %---------------------------------------------------------------------------
 %
 %   Computes oxygen saturation concentration at 1 atm total pressure
 %   in mmol/m^3 given the temperature (t, in deg C) and the salinity (s,
 %   in permil)
 %
 %   FROM GARCIA AND GORDON (1992), LIMNOLOGY and OCEANOGRAPHY.
 %   THE FORMULA USED IS FROM PAGE 1310, EQUATION (8).
 %
 %   *** NOTE: THE "A_3*TS^2" TERM (IN THE PAPER) IS INCORRECT. ***
 %   *** IT SHOULDN'T BE THERE.                                ***
 %
 %   O2SAT IS DEFINED BETWEEN T(freezing) <= T <= 40(deg C) AND
 %   0 permil <= S <= 42 permil
 %   CHECK VALUE:  T = 10.0 deg C, S = 35.0 permil,
 %   O2SAT = 282.015 mmol/m^3
 %
 %---------------------------------------------------------------------------

a_0 = 2.00907;
a_1 = 3.22014;
a_2 = 4.05010;
a_3 = 4.94457;
a_4 = -2.56847e-1;
a_5 = 3.88767;
b_0 = -6.24523e-3;
b_1 = -7.37614e-3;
b_2 = -1.03410e-2;
b_3 = -8.17083e-3;
c_0 = -4.88682e-7;

% Avoid bad SSS
S_LOC = max(1e-4,sss);

% Convert SST to Kelvin
TS = log(((273.16+25.0)-sst) / (273.16+sst));

% Calc
o2sat = exp(a_0+TS*(a_1+TS*(a_2+TS*(a_3+TS*(a_4+TS*a_5)))) + ...
         S_LOC*((b_0+TS*(b_1+TS*(b_2+TS*b_3))) + S_LOC*c_0 ));

% Convert from ml/l to mmol/m3
o2sat = o2sat .* 44.6596; 



