%Kepler Orbit Anomaly true mean
% Richard Rieber
% September 27, 2006
% rrieber@gmail.com
% 
% Revision 9/25/07 - Fixed a grusome error in the default tolerance.
%                    Changed from 10^8 radians to 10^-8 radians.  Whoops.
%
% function nu = nuFromM(M,ecc,tol)
% 
% Purpose:  This function calculates the true anomaly (nu) of a position in an orbit given
%           the mean anomaly of the position (M) and the eccentricity (ecc) of the orbit.
%           This uses another function, calcEA.
%           
% Inputs:  M   - mean anomaly of position in radians
%          ecc - eccentricity of orbit
%          tol - A tolerance for calculating the eccentric anomaly (see help calcEA.m)
%                Default is 10^-8 radians [OPTIONAL]
% 
% Output:  nu  - true anomaly of position in radians

function nu = nuFromM(M,ecc,tol)

if nargin < 2 || nargin > 3
    error('Incorrect number of inputs, see help nuFromM.m')
elseif nargin == 2
    tol = 10^-8;
end

E = CalcEA(M,ecc,tol);  %Determining eccentric anomaly from mean anomaly

% Since tan(x) = sin(x)/cos(x), we can use atan2 to ensure that the angle for nu
% is in the correct quadrant since we know both sin(nu) and cos(nu).  [see help atan2]
nu = atan2((sin(E)*(1-ecc^2)^.5),(cos(E)-ecc));