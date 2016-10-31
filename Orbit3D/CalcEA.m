% Orbit eccentric anomaly, Kepler's equation keplers equation
% Richard Rieber
% 1/23/2005
% rrieber@gmail.com
% 
% Revision 8/21/07: Fixed typo in line 38 if statement
%                   Added H1 line for lookfor functionality
%
% function E = CalcEA(M,ecc,tol)
% 
% Purpose: Solves for eccentric anomaly, E, from a given mean anomaly, M,
%          and eccentricty, ecc.  Performs a simple Newton-Raphson iteration
%
% Inputs: o M   - Mean anomaly in radians.
%         o ecc - Eccentricity of the orbit.
%         o tol - a tolerance at which to terminate iterations; Default
%                 is 10^-8 radians. [OPTIONAL]
%
% Outputs: o E  - The eccentric anomaly in radians.
%
% E = CalcEA(M,ecc) uses default tolerances
%
% E = CalcEA(M,ecc,tol) will use a user specified tolerance, tol
% 

function E = CalcEA(M,ecc,tol)

%Checking for user inputed tolerance
if nargin == 2
    %using default value
    tol = 10^-8;
elseif nargin > 3
    error('Too many inputs.  See help CalcEA')
elseif nargin < 2
    error('Too few inputs.  See help CalcEA')
end

if (M > -pi && M < 0) || M > pi
    E = M - ecc;
else
    E = M + ecc;
end

Etemp  = E + (M - E + ecc*sin(E))/(1-ecc*cos(E));
Etemp2 = E;

while abs(Etemp - Etemp2) > tol
    Etemp = Etemp2;
    Etemp2 = Etemp + (M - Etemp + ecc*sin(Etemp))/(1-ecc*cos(Etemp));
end

E = Etemp2;