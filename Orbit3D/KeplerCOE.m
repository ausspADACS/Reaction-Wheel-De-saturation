%Orbit Kepler position velocity
% Richard Rieber
% November 9, 2006
% rrieber@gmail.com
% 
% Revision 08/21/07: Added H1 line for lookfor functionality
%          11/11/08: Fixed a really stupid bug in the call to elorb.m
%          09/27/09: Added an input for mu and tolerance
%                    Extensive modifications to make this work for
%                    equatorial circular orbits, equatorial elliptical
%                    orbits, and inclined circular orbits
%
% function [R,V] = KeplerCOE(Ro,Vo,dT,U,tol)
%
% Purpose: This function calculates position and velocity
%          at a given time based on an initial position
%          and velocity of an orbiting object.
%
% Inputs:  Ro  - Initial position of length 3 (km)
%          Vo  - Initial velocity of length 3 (km/s)
%          dT  - A time step at which to calculate the new R and V vectors (sec)
%          U   - Gravitational constant of body being orbited (km^3/s^2).  Default is Earth
%                at 398600.4415 km^3/s^2.  [OPTIONAL]
%          tol - Tolerance for CalcEA.m, defaults to 10^-8 [OPTIONAL]
%
% Outputs: R - Position at time dT of length 3 (km)
%          V - Velocity at time dT of length 3 (km/s)
%
% NOTE:  This function uses the subfunction CalcEA.m, randv.m, and elorb.m

function [R,V] = KeplerCOE(Ro,Vo,dT,U,tol)

if nargin < 3 || nargin > 5
    error('Incorrect number of inputs.  See help KeplerCOE')
elseif nargin == 3
    U = 398600.4415;  %km^3/s^2 Gravitational Constant of Earth
    tol = 10^-8;
elseif nargin == 4
    tol = 10^-8;
end

if length(Ro) ~= 3
    error('Position vector must be of length 3.  See help KeplerCOE')
elseif length(Vo) ~= 3
    error('Velocity vector must be of length 3.  See help KeplerCOE')
end

% Calculating kepler orbital elements at given position
[a,ecc,inc,O,w,nu,w_true,u_true,lambda_true] = elorb(Ro,Vo,U);
% note: * = unitless

% Mean motion of orbit
n = (U/a^3)^.5; %rad/s

%%%%% Equatorial elliptical circular orbit
if (abs(inc) < tol) && (ecc > tol)
    inc = 0;
    O   = 0;
    w   = w_true;
end

% Checking for circular orbits
if ecc < tol
    ecc = 0;
    
    %%%% Equatorial circular orbit
    if inc == 0
        % Assuming that the eccentric anomaly, mean anomaly, true anomaly, 
        % and true longitude of a circular, equatorial orbit are all equal,
        % we can say that the new true anomaly is simply the true longitude
        % plus the mean motion times dT
        nu2 = lambda_true + n*dT;
        w   = 0;
        O   = 0;
        
    %%%% Inclined circular orbit
    else
        % Assuming that the eccentric anomaly, mean anomaly, true anomaly,
        % and argument of latitude of an inclined circular orbit are all
        % equal, we can say that the new true anomaly is simply the
        % argument of latitude plus the mean motion times dT
        nu2 = u_true + n*dT;
        w   = 0;
    end
else
%   Initial eccentric anomaly
    Eo    = atan2(sin(nu)*(1-ecc^2), ecc+cos(nu)); %rad
%   Initial mean anomaly
    Mo    = Eo-ecc*sin(Eo); %rad
%   Final mean anomaly
    M     = Mo + n*dT; %rad
%   Final eccentric anomaly    
    E     = CalcEA(M,ecc,tol); %rad
%   Final true anomaly    
    nu2   = atan2(sin(E)*(1-ecc^2)^.5,cos(E)-ecc); %rad
end

[R,V] = randv(a,ecc,inc,O,w,nu2,U); %[km, km/s]
