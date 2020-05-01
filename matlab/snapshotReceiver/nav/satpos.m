% [rSat, clkBias, vSat, clkDrift] = satpos( t, eph )
%
% Computation of the satellite position, velocity, clock bias and clock
% drift using the broadcast ephemerides.
%
% Parameters:
% t.............. time for which the satellite position should be
%                 determined
% eph............ satellite's ephemeris parameters
%
% Returns:
% rSat........... Satellite position in ECEF coordinates (3x1) [m]
% clkBias........ Clock bias of the satellite at time t [s]
% vSat........... Velocity of the satellite in ECEF coordinates (3x1) [m/s]
% clkDrift....... Satellite's clock drift [s/s]
%
function [rSat, clkBias, vSat, clkDrift] = satpos(t, eph)

% INITIALIZATION OF CONSTANTS =============================================
Omegae_dot     = 7.2921151467e-5;  % Earth rotation rate, [rad/s]
GM             = 3.986004418e14;   % Earth's universal gravitational parameter, [m^3/s^2]
F              = -4.442807633e-10; % Constant, [sec/(meter)^(1/2)]

% POSITION DETERMINATION ==================================================
a   = eph.rootA * eph.rootA;
tk  = t - eph.t_oe;

n0  = sqrt(GM / a^3);  % Initial mean motion
n   = n0 + eph.deltan; % Mean motion

M   = rem( eph.M_0 + n*tk + 2*pi, 2*pi ); % Mean anomaly

% Recursive computation of eccentric anomaly (no for loop to be more eff.)
E = M + eph.e*sin(M);
E = M + eph.e*sin(E);
E = M + eph.e*sin(E);
E = M + eph.e*sin(E);
E = M + eph.e*sin(E);
E = M + eph.e*sin(E);

dtr = F * eph.e * eph.rootA * sin(E); % Relativistic correction
nu = atan2(sqrt(1 - eph.e^2) * sin(E), cos(E)-eph.e); % True anomaly

phi = rem( nu + eph.omega, 2*pi );

% Argument of latitude
u = phi + eph.C_uc * cos(2*phi) + eph.C_us * sin(2*phi);

r = a * (1 - eph.e*cos(E)) + eph.C_rc * cos(2*phi) + eph.C_rs * sin(2*phi);

i = eph.i_0 + eph.iDot * tk + eph.C_ic * cos(2*phi) + eph.C_is * sin(2*phi);

Omega = eph.Omega_0 + (eph.omegaDot - Omegae_dot)*tk - Omegae_dot * eph.t_oe;


% Satellite coordinates
rSat = [ cos(u)*r * cos(Omega) - sin(u)*r * cos(i)*sin(Omega);
    cos(u)*r * sin(Omega) + sin(u)*r * cos(i)*cos(Omega);
    sin(u)*r * sin(i) ];

% CLOCK OFFSET COMPUTATION ================================================
if nargout >= 2
    dt = t - eph.t_oc;
    
    clkBias = (eph.a_f2 * dt + eph.a_f1) * dt + eph.a_f0 - eph.T_GD + dtr;
end


% VELOCITY COMPUTATION ====================================================
if nargout >= 3
    Mdot = n;
    Edot = Mdot / ( 1-eph.e*cos(E) );
    
    nudot = sin(E)*Edot*(1+eph.e*cos(nu)) / ( (1-cos(E)*eph.e)*sin(nu) );
    phidot = nudot;
    
    udot = phidot + 2 * ( eph.C_us*cos(2*phi)-eph.C_uc*sin(2*phi) ) * phidot;
    
    rdot = a*eph.e*sin(E)*Edot + 2*( eph.C_rs*cos(2*phi)-eph.C_rc*sin(2*phi) )*phidot;
    idot = eph.iDot + 2*( eph.C_is*cos(2*phi)-eph.C_ic*sin(2*phi) )*phidot;
    
    Omegadot = eph.omegaDot - Omegae_dot;
    
    xidot = rdot * cos(u) - r * sin(u) * udot;
    yidot = rdot * sin(u) + r * cos(u) * udot;
    yi = r * cos(u);
    
    vSat = [ xidot*cos(Omega) - yidot*cos(i)*sin(Omega) + yidot*sin(i)*sin(Omega)*idot - rSat(2)*Omegadot;
        xidot*sin(Omega) + yidot*cos(i)*cos(Omega) - yidot*sin(i)*idot*cos(Omega) + rSat(1)*Omegadot;
        yidot*sin(i) + yi*cos(i)*idot ];
end


% CLOCK DRIFT COMPUTATION =================================================
if nargout >= 4
    clkDrift = eph.a_f2*dt + eph.a_f1;
end
