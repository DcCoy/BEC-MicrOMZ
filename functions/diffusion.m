function diff = diffusion(c, kappa_z, dz)
% Diffusion function

% Pad concentration array with zeros at surface/bottom
c = [transpose(zeros(size(c,2))); c; transpose(zeros(size(c,2)))];

% Get length of Kv profile
kappa_z_len = size(kappa_z,1);

% Get diffusion into upper cell
f_up = kappa_z(1:kappa_z_len-1) .* (c(2:kappa_z_len,:) - c(1:kappa_z_len-1,:)) ./ dz;

% Get diffusion from lower cell 
f_down = kappa_z(2:kappa_z_len) .* (c(3:kappa_z_len+1,:) - c(2:kappa_z_len,:)) ./ dz;

% Get divergence 
diff = (f_down - f_up) ./ dz;

