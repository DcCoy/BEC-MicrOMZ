function schmidt_n2 = cschmidt_n2o(sst);
%---------------------------------------------------------------------------
%   Compute Schmidt number of N2 in seawater as function of SST
%
%   ref : Sarmiento and Gruber 2006 (Table 3.3.1)
%---------------------------------------------------------------------------

a = 2206.1;
b = 144.86;
c = 4.5413;
d = 0.056988;

schmidt_n2 = a + sst.*(-b + sst.*(c + sst.*(-d)));
