% ionoDelay = ionoKlob(TOW, alpha, beta, az, el, lat, lon )
% 
% computes the ionospheric delay (in meters) using the standard GPS
% Klobuchar model. 
% 
% Parameters:
% TOW......... Time of week (GPS time)
% alpha....... Alpha parameters ordered as [alpha0 alpha1 ... alpha3]
% beta........ Same but beta
% az.......... Azimuths of the satellites [rad]
% el.......... Elevations of the satellites [rad]
% lat......... Latitude of the receiver [rad]
% lon......... Longitutde of the receiver [rad]
% 
% Returns: 
% ionoDelay... Ionospheric slant delay  [m]
% 
function ionoDelay = ionoKlob(TOW, alpha, beta, az, el, lat, lon)

if length(alpha) < 4 || length(beta) < 4
    ionoDelay = 0;
    return;
end

c = 299792458;

F = ( 1.0 + 16.0*(0.53 - el/pi).^3 ) * c;

psi = 0.0137./(el/pi + 0.11) - 0.022;

phii = lat/pi + psi.*cos(az);
phii( (phii > 0.416) )  = 0.416;
phii( (phii < -0.416) ) = -0.416;

lambdai = lon/pi + psi.*sin(az)./cos(phii*pi);
phim = phii + 0.064*cos((lambdai-1.617)*pi);

t = mod( 4.32*1e4*lambdai + TOW, 24*3600 );

if size(beta,1) > size(beta,2)
    beta = beta';
end
if size(alpha,1) > size(alpha,2)
    alpha = alpha';
end

PER = beta(1) + beta(2)*phim + beta(3)*phim.^2 + beta(4)*phim.^3;
PER( PER < 72000 ) = 72000;
x = 2*pi*(t-50400)./PER;

AMP = alpha(1) + alpha(2)*phim + alpha(3)*phim.^2 + alpha(4)*phim.^3;
AMP( AMP < 0 ) = 0;

ionoDelay = F.*5.0e-9;

if abs(x) < 1.57
    ionoDelay = ionoDelay + F.*AMP*(1-(x.^2/2)+(x.^4/24));
end
