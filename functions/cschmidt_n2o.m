function schmidt_n2o = cschmidt_n2o(sst);
%---------------------------------------------------------------------------
%   Compute Schmidt number of N2O in seawater as function of SST
%
%   ref : Sarmiento and Gruber 2006 book (OBCD, table 3.3.1 page 85
%---------------------------------------------------------------------------

a = 2301.1;
b = 151.1;
c = 4.7364;
d = 0.059431;

schmidt_n2o = a + sst.*(-b + sst.*(c + sst.*(-d)));
