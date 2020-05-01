% [az, el] = getAzEl(x, dx [, lat, lon ])
%
% Transformation of vector dx into topocentric coordinate system with
% origin at X.
%
% Parameters:
% x............... vector origin coordinates (3x1, ECEF) [m]
% dx.............. vector to be transformed (3x1, ECEF) [m]
% lat, lon........ (optional) latitude and longitude of x [rad]
%
% Returns:
% az.............. azimuth angle [rad]
% el.............. elevation angle [rad]
%
function [az, el] = getAzEl(x, dx, lat, lon)

if ~exist('lat', 'var') || ~exist('lon', 'var')
    [lat, lon] = xyz2llh(x);
end

cl  = cos(lon);
sl  = sin(lon);
cb  = cos(lat);
sb  = sin(lat);

F   = [-sl -sb*cl cb*cl;
    cl -sb*sl cb*sl;
    0    cb   sb];

dxLocal = F' * dx;
E   = dxLocal(1);
N   = dxLocal(2);
U   = dxLocal(3);

hor_dis = sqrt(E^2 + N^2);

az = atan2(E, N);
el = atan2(U, hor_dis);

if az < 0
    az = az + 2*pi;
end
