% [t0,rEst,clkbEst,vEst,clkdEst] = coarsePVT(t0, r0, eph, codePhase, Tambg,
%                                            fD, cno, alpha, beta, doy)
%
% Coarse snapshot positioning after Van Diggelen (A-GPS book). This
% includes estimation of a coarse time-offset, a fine clock bias and the
% position itself. The algorithm also resolves the ambiguities arising from
% the duration of one spreading code sequence. In contrast to what's
% covered in literature, this approach does not make use of the Doppler.
% Since the receiver may be moving (quite fast...) this information is not
% reliable. Therefore the Doppler-part steming only from the satellite
% movements (which is actually the one we're interested in) is estimated
% on-the-fly.
%
% Although only a coarse position solution is seeked in this function, the
% atmospheric corrections are applied - if available.
%
% Parameters:
% t0............ Coarse time of the receiver (TOW) [s]
% r0............ Coarse receiver position (ECEF) [m]
% eph........... Ephemerides for the satellites (vector of structs)
% codePhase..... Code phase measurements [s]
% Tambg......... Corresponding code phase ambiguities [s]
% fD............ Doppler frequencies [Hz]
% cno........... C/N0 [dB-Hz]
% alpha, beta... Klobuchar iono model parameters
% doy........... Day of year
%
% Returns:
% t0............ Time of the measurements (TOW) [s]
% rEst.......... Receiver position (ECEF) [m] (3x1)
% clkbEst....... Receiver clock bias [m]
% vEst.......... Receiver velocity (ECEF) [m/s] (3x1)
% clkdEst....... Receiver clock drift [m/s]
%
function [t0,rEst,clkbEst,vEst,clkdEst] = coarsePVT(t0, r0, eph, ...
    codePhase, Tambg, fD, cno, alpha, beta, doy)

% INITIALIZATION ----------------------------------------------------------
N_it = 5;
K    = length(codePhase);

c = 299792458; % speed of light [m/s]

rEst      = r0(:);
clkbEst = 0;

% pseudorange measurements and their ambiguities
prMeas = (Tambg(:) - codePhase(:)) * c;
ambg = Tambg(:) * c;

% measurement variance
prVar = 10.^(-cno/10);

H       = ones(K, 5); % the geometry matrix (fifth column for the coarse timing)
velSat  = zeros(3, K);
clkdSat = zeros(K, 1);
prEst   = zeros(K, 1);
r_d     = zeros(K, 1); % geometric range plus the satellite's clock bias
Nms     = NaN(K, 1);   % millisecond ambiguities

% FIRST ITERATION ---------------------------------------------------------
% First estimation of the pseudorange. If done outside the loop, the later
% loop will be more easy implementable.
for k=1:K
    [rSat, clkbSat] = satpos(t0-.065, eph(k));
    prEst(k) = norm(rEst-rSat) - c*clkbSat;
    
    % arbitrary ambiguity fixing for the first (=reference) satellite
    if k==1
        Nms(k) = round((prEst(k)-prMeas(k)) / ambg(k));
    end
end

% NEWTON'S METHOD ---------------------------------------------------------
for n=2:N_it
    % in the last two iterations the atmospheric corrections are
    % applied (if they are switched on)
    if n >= N_it-2
        compensateAtmosphere = 1;
        [lat,lon,alt] = xyz2llh(rEst);
    else
        compensateAtmosphere = 0;
    end
    
    for k=1:K
        % satellite positioning (position, clock and velocity)
        [rSat,clkbSat,velSat(:,k),clkdSat(k)] = ...
            satpos(t0-prEst(k)/c, eph(k));
        rSat = earthRotComp(rSat, norm(rEst-rSat)/c);
        
        % update of the geometry matrix
        d = rEst - rSat;
        H(k,1:3) = d' / norm(d);
        H(k,5) = -H(k,1:3) * velSat(:,k);
        r_d(k) = norm(d) - c*clkbSat;
        
        corr = 0;
        if compensateAtmosphere
            [az,el] = getAzEl(rEst, rSat-rEst);
            prVar(k) = 10^(-cno(k)/10) ./ sin(el)^2;
            
            % ionosphere
            corr = corr + ionoKlob(t0, alpha, beta, az, el, lat, lon);
            
            % troposphere
            m = mopsMappingFunc(el);
            [tauw,taud] = mopsZenithDelay(lat, alt, doy);
            corr = corr + m*(tauw + taud);
        end
        
        % estimation of the pseudorange
        prEst(k) = norm(d) + clkbEst - c*clkbSat + corr;
        
        % ambiguity resolution w.r.t. to the ambiguity fixed for the first
        % (i.e. reference) satellite
        if k > 1
            Nms(k) = round((prMeas(1)-prMeas(k) + Nms(1)*ambg(1) ...
                + r_d(k)-r_d(1))/ambg(k));
        end
    end
    
    % Update of the receiver position
    W = diag(1./prVar);
    dx = (H'*W*H) \ H'*W*(prMeas+Nms.*ambg - prEst);
    if any(~isreal(dx))
        break;
    end
    
    rEst      = rEst + dx(1:3);
    clkbEst = clkbEst + dx(4);
    t0        = t0 + dx(5);
end


% VELOCITY ESTIMATION -----------------------------------------------------
% Receiver velocity estimation
fc = 1.57542e9;
rrMeas = -(c/fc)*fD(:) + sum(velSat'.*H(:,1:3), 2) + c*clkdSat;
H4 = H(:,1:4);
xi_hat_dot = (H4'*W*H4) \ H4'*W*rrMeas;

vEst = xi_hat_dot(1:3);
clkdEst = xi_hat_dot(4);
