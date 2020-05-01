% rSat = earthRotComp(rSat, T)
%
% Satellite position compensation for the earth rotation
%
% Parameters:
% rSat...... Satellite position (ECEF) [m] (3x1)
% T......... Time to compensate the rotation for [s]
%
% Returns:
% rSat...... Compensated satellite position (ECEF) [m] (3x1)
%
function rSat = earthRotComp(rSat, T)

% WGS 84 earth rotation rate (as copied from GPS ICD)
earthRotVel = 7.2921151467e-5;
a = earthRotVel * T;

rSat = [cos(a) sin(a) 0;
    -sin(a)  cos(a) 0;
    0        0      1] * rSat;
