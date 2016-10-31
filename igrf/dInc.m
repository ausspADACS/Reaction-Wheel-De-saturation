%Inclination change orbit gravity
% Richard Rieber
% October 22, 2006
% rrieber@gmail.com
%
% function dV = dInc(V,dI,fpa)
%
% Purpose:  This function calculates the change of velocity needed to
%           change the inclination of an orbit given the current
%           velocity in km/s and the change in inclination in radians.
%           The flight path angle may also be provided if the orbit is not
%           circular.
%
% Inputs:  o V   - velocity of the orbit at which the inclination will be
%                  changed in km/s
%          o dI  - the chnage of inclination needed in radians
%          o fpa - Flight path angle of the orbit relative to the local
%                  horizon in radians [OPTIONAL]
%
% Outputs: o dV - The necessary change in velocity needed to complete
%                 the maneuver in km/s

function dV = dInc(V,dI,fpa)

if nargin < 2
    error('Too few inputs.  See help dInc')
elseif nargin == 2
    fpa = 0;
elseif nargin > 3
    error('Too many inputs.  See help dInc')
end

dV = 2*V*cos(fpa)*sin(dI/2);