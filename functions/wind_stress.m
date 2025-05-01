function ws_used = wind_stress(frc,params);
% Script to calculate wind-stress based on forcing fields

% Get ambient air-density (kg/m3)
rho_air = 1.2;

% Get local variables
coef_1  = 0.0027;     % (m/s)
coef_2  = 0.000142;   % (non-dimensional)
coef_3  = 0.0000764;  % (s/m)
c_d     = 1.7e-3;     % drag coefficient for initial u

% Solve for WS in (m/s)
ustar_squared = (sqrt(frc.sustr.^2 + frc.svstr.^2) .* params.bec.rho0) ./ rho_air; 
ws_used = sqrt(ustar_squared ./ c_d);

% Perform 3 Newton iterations
ws_used = ws_used - ...
    (ws_used*(coef_1 + ws_used*(coef_2 + ws_used*coef_3)) - ustar_squared) / ...
    (coef_1 + ws_used*(2*coef_2 + ws_used*3*coef_3));

ws_used = ws_used - ...
    (ws_used*(coef_1 + ws_used*(coef_2 + ws_used*coef_3)) - ustar_squared) / ...
    (coef_1 + ws_used*(2*coef_2 + ws_used*3*coef_3));

ws_used = ws_used - ...
    (ws_used*(coef_1 + ws_used*(coef_2 + ws_used*coef_3)) - ustar_squared) / ...
    (coef_1 + ws_used*(2*coef_2 + ws_used*3*coef_3));
