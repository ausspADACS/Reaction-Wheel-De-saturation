% Orbit Hohmann transfer
% Richard Rieber
% October 17, 2006
% rrieber@gmail.com
%
% Revision 8/21/07: Supressed outputs from various calculations.
%                   Added H1 line for lookfor functionality
%
% Purpose:  This function calculates the two changes of velocity, time, 
%           and semi-major axis of a Hohmann transfer between two circular 
%           orbits.
%
% Inputs:  o R_init - Radius of initial circular orbit in km
%          o R_fin  - Radius of final circular orbit in km
%          o U      - Gravitational constant of body being orbited (km^3/s^2).
%                     Default is Earth at 398600.4415 km^3/s^2.  OPTIONAL
%
% Outputs: o dV   - A 1x2 vector of the changes of velocity needed at the initial burn
%                   and final burn in km/s
%          o T    - The time needed to complete the transfer
%          o a_tx - The semi-major axis of the tranfer orbit.
%

function [dV,T,a_tx] = Hohmann(R_init,R_fin,U)

if nargin < 2
    error('Too few inputs.  See help Hohmann')
elseif nargin > 3
    error('Too many inputs.  See help Hohmann')
elseif nargin == 2
    U = 398600.4415; %km^3/s^2
end

%Initial circular velocity
v_init = (U/R_init)^.5; %km/s

%Final circular velocity
v_fin = (U/R_fin)^.5; %km/s

a_tx = (R_init + R_fin)/2;  %Semi-major axis of transfer orbit (km)

% Required velocities at periapse and apoapse of the transfer orbit
V_trans_a = (2*U/R_init - U/a_tx)^.5; %Initial transfer vel. needed (km/s)
V_trans_b = (2*U/R_fin - U/a_tx)^.5;  %Final transfer vel. (km/s)

%Change in velocities needed
dVa = V_trans_a - v_init; %km/s
dVb = v_fin - V_trans_b;  %km/s

dV = [dVa, dVb]; %km/s

%Transfer time
T = pi*(a_tx^3/U)^.5; %sec