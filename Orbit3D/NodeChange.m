%Node change right ascension of the ascending node RAAN raan orbit
% Richard Rieber
% rrieber@gmail.com
% December 18, 2006
% 
% Revision 8/21/07: Added H1 line for lookfor functionality.
%
% function [dV,u_init,u_fin] = NodeChange(dO,inc,Vinit)
%
% Purpose: This function calculates the necessary dV and the initial
%          and final argument of longitude for a change in the right ascension
%          of the ascending node.
%
% Inputs: o dO    - The needed change in the right ascension of the ascending
%                   node (Omega) in radians.
%         o inc   - The inclination (i) of the orbit in radians.
%         o Vinit - The magnitude of the current velocity in km/s (or whatever
%                   units you like...perhaps furlongs per fortnight ^_^)
%
% Outputs: o dV     - The change in velocity needed to complete the manuever in km/s.
%          o u_init - The initial argument of longitude at which the burn will take place
%                     in radians.
%          o u_fin  - The argument of longitude of where the burn occurs in the final
%                     orbit in radians.


function [dV,u_init,u_fin] = NodeChange(dO,inc,Vinit)

if nargin ~= 3
    error('Incorrect number of inputs.  See help NodeChange.m')
end

% The angle between the initial and final velocity vectors in radians
Vang = acos(cos(inc)^2 + sin(inc)^2*cos(dO));

% The necessary change in velocity needed in km/s
% (or whatever units you provide as an input)
dV = 2*Vinit*sin(Vang/2);

% The initial and final argument of longitudes in radians
u_init = acos(tan(inc)*(cos(dO) - cos(Vang))/sin(Vang));
u_fin  = acos(cos(inc)*sin(inc)*(1-cos(dO))/sin(Vang));