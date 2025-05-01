function schmidt_o2 = cschmidt_o2(sst);
%---------------------------------------------------------------------------
%   Compute Schmidt number of O2 in seawater as function of SST
%
%   ref : Sarmiento and Gruber 2006 book (OBCD, table 3.3.1 page 85
%---------------------------------------------------------------------------

a = 1638.0;
b = 81.83;
c = 1.483;
d = 0.008004;

schmidt_o2 = a + sst.*(-b + sst.*(c + sst.*(-d)));
