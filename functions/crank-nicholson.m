function [coeff1 coeff2 coeff3] = crank-nicholson(phy,opt);
% Function to solve advection and diffusion using
% Crank-Nicholson method

% For advection velocity and diffusion coefficient fixed in time, calculate here
% terms for the numerical advection-diffusion solver. For time-dependent w and Kv
% move these terms inside the time loop
alpha = phy.wup(2:end-1) * opt.dt / (2*grid.Hz(2:end-1));
beta  = - opt.dt / (2*grid.Hz(1:end-2)) * (phy.wup(1:end-2) - phy.wup(3:end));
gamma = phy.Kv(2:end-1) * opt.dt / (grid.Hz(2:end-1))^2;
delta =   opt.dt / (4*grid.Hz(1:end-2)) * (phy.Kv(1:end-2) - phy.Kv(3:end));

% Integration coefficients for the tracer at k,k+1,k-1 vertical levels:
coeff1 = 1 + beta - 2*gamma;
coeff2 =     alpha +  gamma - delta;
coeff3 =   - alpha +  gamma + delta;
