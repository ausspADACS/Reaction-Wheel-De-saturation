%Kepler Orbit Anomaly true time periapse perigee
% Richard Rieber
% September 27, 2006
% rrieber@gmail.com
%
% Revision 9/25/07 - Fixed a grusome error in the default tolerance.
%                    Changed from 10^8 radians to 10^-8 radians.  Whoops.
% Revision 10/1/09 - Added H1 help line and fixed some references to
%                    nuFromM, which were copy&paste errors
%
% function nu = nuFromM(Tp,ecc,n,tol)
% 
% Purpose:  This function calculates the true anomaly (nu) of a position in an orbit given
%           the time since periapse passage in seconds, the eccentricity (ecc) of the orbit,
%           the mean motion (n) in radians/sec.
%           This uses another function, CalcEA.
%           
% Input: o Tp  - Time since periapse passage (seconds)
%        o ecc - Eccentricity of orbit
%        o n   - Mean motion of orbiting body (radians/sec)
%        o tol - A tolerance for calculating the eccentric anomaly (see help CalcEA.m)
%                Default is 10^-8 radians [OPTIONAL]
% 
% Output: o nu - True anomaly of position in radians

function nu = nuFromTp(Tp,ecc,n,tol)

if nargin < 3 || nargin > 4
    error('Incorrect number of inputs, see help nuFromTp.m')
elseif nargin == 3
    tol = 10^-8;
end

% Calculating the mean anomaly from the mean motion and time [see pg 53-54 of
% David A. Vallado's Fundamentals of Astrodynamics and Applications]
M = n*Tp;

E = CalcEA(M,ecc,tol);  %Determining eccentric anomaly from mean anomaly

% Since tan(x) = sin(x)/cos(x), we can use atan2 to ensure that the angle for nu
% is in the correct quadrant since we know both sin(nu) and cos(nu).  [see help atan2]
nu = atan2((sin(E)*(1-ecc^2)^.5),(cos(E)-ecc));