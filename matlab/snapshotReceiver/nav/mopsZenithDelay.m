% [tauw,taud] = mopsZenithDelay(lat, alt, d)
%
% Returns the estimated tropospheric zenith wet and dry delay after MOPS, 
% RTCA DO-229D, in meters.
%
% Parameters:
% lat....... receiver's latitude [rad]
% alt....... receiver's altitude [m]
% d......... day of year
% 
% Returns:
% tauw...... wet zenith delay [m]
% taud...... dry zenith delay [m]
% 
function [tauw,taud] = mopsZenithDelay(lat, alt, d)

latgrid = [ 15 30 45 60 75 ];
P0 = [ 1013.25 1017.25 1015.75 1011.75 1013.00 ]; % mbar
T0 = [  299.65  294.15  283.15  272.15  263.65 ]; % K
e0 = [   26.31   21.79   11.66    6.78    4.11 ]; % mbar
b0 = [ 6.30e-3 6.05e-3 5.58e-3 5.39e-3 4.53e-3 ]; % K/m
l0 = [    2.77    3.15    2.57    1.81    1.55 ]; % dimensionless

dP = [ 0.00 -3.75 -2.25 -1.75 -0.50 ]; % mbar
dT = [ 0.00  7.00 11.00 15.00 14.50 ]; % K
de = [ 0.00  8.85  7.24  5.36  3.39 ]; % mbar
db = [ 0.00 0.25e-3 0.32e-3 0.81e-3 0.62e-3 ]; % K/m
dl = [ 0.00  0.33  0.46  0.74  0.30 ]; % dimensionless

% find closest two grid points
lat = lat*180/pi;
if lat > 0
    Dmin = 28;
else
    Dmin = 211;
    lat  = -lat;
end

common = cos(2*pi*(d-Dmin)/365.25);

if lat <= 15
    P = P0(1) - dP(1)*common;
    T = T0(1) - dT(1)*common;
    e = e0(1) - de(1)*common;
    b = b0(1) - db(1)*common;
    l = l0(1) - dl(1)*common;
elseif lat >= 75
    P = P0(5) - dP(5)*common;
    T = T0(5) - dT(5)*common;
    e = e0(5) - de(5)*common;
    b = b0(5) - db(5)*common;
    l = l0(5) - dl(5)*common;
else
    i1 = find(lat-latgrid > 0, 1, 'last');
    i2 = find(latgrid-lat > 0, 1, 'first');
    lat1 = latgrid(i1);
    lat2 = latgrid(i2);
            
    P0_i = P0(i1) + (P0(i2)-P0(i1))*(lat-lat1)/(lat2-lat1);
    dP_i = dP(i1) + (dP(i2)-dP(i1))*(lat-lat1)/(lat2-lat1);
    P = P0_i - dP_i*common;

    T0_i = T0(i1) + (T0(i2)-T0(i1))*(lat-lat1)/(lat2-lat1);
    dT_i = dT(i1) + (dT(i2)-dT(i1))*(lat-lat1)/(lat2-lat1);
    T = T0_i - dT_i*common;

    e0_i = e0(i1) + (e0(i2)-e0(i1))*(lat-lat1)/(lat2-lat1);
    de_i = de(i1) + (de(i2)-de(i1))*(lat-lat1)/(lat2-lat1);
    e = e0_i - de_i*common;

    b0_i = b0(i1) + (b0(i2)-b0(i1))*(lat-lat1)/(lat2-lat1);
    db_i = db(i1) + (db(i2)-db(i1))*(lat-lat1)/(lat2-lat1);
    b = b0_i - db_i*common;

    l0_i = l0(i1) + (l0(i2)-l0(i1))*(lat-lat1)/(lat2-lat1);
    dl_i = dl(i1) + (dl(i2)-dl(i1))*(lat-lat1)/(lat2-lat1);
    l = l0_i - dl_i*common;
end

k1 = 77.604; % K/mbar
k2 = 382000; % K^2/mbar
Rd = 287.054; % J/kg/K
gm = 9.784; % m/s^2
g  = 9.80665; % m/s^2
Zwet = 1e-6 * k2 * Rd * e / ((gm*(l+1)-b*Rd) * T);
Zdry = 1e-6 * k1 * Rd * P / gm;

tauw = (1-b*alt/T)^(-1+(l+1)*g/(Rd*b))*Zwet;
taud = (1-b*alt/T)^(g/(Rd*b))*Zdry;
