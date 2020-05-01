% m = mopsMappingFunc(el)
% 
% Returns the MOPS mapping functions as defined in RTCA DO-229D. 
% 
% Parameters:
% el......... elevations to the satellites [rad]
% 
% Returns:
% m.......... mapping function
% 
function m = mopsMappingFunc(el)

m = 1.001./sqrt( .002001+sin(el).^2 );

if el*180/pi < 4
    m = m.*( 1 + .015*max([0, 4*pi/180-el]).^2 );
end
