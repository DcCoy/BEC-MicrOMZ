function f = calc_f(dGa, dGd, dGs, ep)
% f = calc_y(dGo, dGr, dGs, ep)
% where:
% OUTPUTS: 
% f = fraction of electrons routed to cell synthesis
%
% INPUTS:
% dGa  = 'delta' Gibbs of oxidation step (electron acceptor)
% dGd  = 'delta' Gibbs of reduction step (electron donor) 
% dGs  = 'delta' Gibbs of cell synthesis (penalty included!!!)
% ep   = efficiency of electron transfer (0 -> 1) 

f = 1/(1-dGs/ep/(dGa-dGd));

% f = 1 / (A + 1)
% where A = (-dGs)/(ep*dGr)
% where dGr = (dGa - dGd)
