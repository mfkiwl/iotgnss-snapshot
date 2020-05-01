% [t,r,v] = snapshotRx(s, fs, fi, t0, r0, rinexNav)
%
% GNSS snapshot receiver a set of input samples.
%
% Firt an acquisition is performed. Second a coarse positioning
% algorithm is initialized. This algorithm solves the millisecond integer
% ambiguity by the algorithm given by Van Diggelen in his book A-GPS. This
% requires initial time and position (t0 and r0) to be known with an
% accuracy of few seconds and few kilometers, respectively.
%
% Parameters:
% s................ Input samples
% fs............... Sampling frequency of input samples [Hz]
% fi............... Intermediate frequency of input samples [Hz]
% t0............... Coarse time of the first sample in the received file
%                   (TOW) [s] (+/-100s)
% r0............... Coarse receiver position at the time of the first
%                   sample (ECEF) [m] (+/-80km)
% rinexNav......... Path to a RINEX navigation file containing the orbits
%
% Returns:
% t................ Time of the first sample (TOW) [s]
% r................ Receiver position (ECEF, 3x1) [m]
% v................ Velocity of the receiver (ECEF, 3x1) [m/s]
%
function [t,r,v] = snapshotRx(s, fs, fi, t0, r0, rinexNav)

addpath('nav');

% --------- ACQUISITION ---------------------------------------------------
[prn, codePhase, Tambg, doppler, cno] = acquisition(s, fs, fi);

% --------- POSITIONING WITH MILLISECOND AMBIGUITY RESOLUTION -------------
% Perform first a rough positioning to resolve the ambiguities
[eph,doy,klobAlpha,klobBeta] = parseRinexNavFile(rinexNav, prn, t0);

% selection of satellites with orbits
ix = find(arrayfun(@(x) ~isempty(x.PRN), eph));

% position estimation
[t,r,~,v] = coarsePVT(t0, r0, eph(ix), codePhase(ix), Tambg(ix), ...
    doppler(ix), cno(ix), klobAlpha, klobBeta, doy);

if nargout <= 1
    [lat,lon,alt] = ecef2lla(r(1), r(2), r(3));
    
    fprintf('%.8f, %.8f, %.1f\n', lat * 180/pi, lon * 180/pi, alt);
    fprintf('https://www.google.com/maps/@%.8f,%.8f,100m/data=!3m1!1e3\n', lat*180/pi, lon*180/pi);
end
