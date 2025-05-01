function adv = upwind(c, wd, dz)
% Upwind advection function

% Pad concentration array with zero at bottom (no flux boundary)
c = [0;c;0];

% Get positive advection
positive_wd = (wd + abs(wd))./2.;

% Get negative advection
negative_wd = (wd - abs(wd))./2.;

% Get advective transport from upper cell
Fu = positive_wd(1:end-1,:).*c(1:end-2,:) + negative_wd(1:end-1,:).*c(2:end-1);

% Get advection from lower cell (if negative, flux into lower cell)
Fd = positive_wd(2:end).*c(2:end-1,:) + negative_wd(2:end).*c(1:end-2,:);

% Get divergence of fluxes
adv = (Fd - Fu)./dz;
