%Kepler orbital elements ECI Position orbit conversion
% Richard Rieber
% September 26, 2006
% rrieber@gmail.com
%
% function [R,V] = randv(a,ecc,inc,Omega,w,nu,U)
% 
% Updates:  9/27/09 - Preallocated memory before the for-loop.
%                     No longer save R_pqw and V_pqw in the for-loop.
%
% Purpose:  This function converts Kepler orbital elements (p,e,i,O,w,nu)
%           to ECI cartesian coordinates in km
%           This function is derived from algorithm 10 on pg. 125 of
%           David A. Vallado's Fundamentals of Astrodynamics and 
%           Applications (2nd Edition)
% 
% WARNING:  o If the orbit is near equatorial and near circular, set w = 0 and Omega = 0
%             and nu to the true longitude.
%           o If the orbit is inclined and near circular, set w = 0 and nu to the argument
%             of latitude.
%           o If the orbit is near equatorial and elliptical, set Omega = 0 and w to the true
%             longitude of periapse.
% 
% Inputs:  a:     Semi-major axis in km of length n
%          ecc:   Eccentricity of length n
%          inc:   Inclination of orbit in radians of length n
%          Omega: Right ascension of ascending node in radians of length n
%          w:     Argument of perigee in radians of length n
%          nu:    True anomaly in radians of length n
%          U:     Gravitational constant of body being orbited (km^3/s^2).  Default is Earth
%                 at 398600.4415 km^3/s^2.  OPTIONAL
%
% Outputs: R:  3 x n array of the radius in km
%          V:  3 x n array of the velocity in km/s

function [R,V] = randv(a,ecc,inc,Omega,w,nu,U)

if nargin < 6 || nargin > 7
    error('Number of inputs incorrect.  See help randv')
elseif nargin == 6
    U = 398600.4415; %Gravitational Constant for Earth (km^3/s^2)
end

if length(a) ~= length(ecc) && length(a) ~= length(inc) && length(a) ~= length(Omega)...
        && length(a) ~= length(w) && length(a) ~= length(nu)
    error('Input vectors are not the same size.  Check inputs')
end

%Pre-allocating memory
R = zeros(3,length(a));
V = zeros(3,length(a));

for j = 1:length(a)
    p = a(j)*(1-ecc(j)^2);
    
    % CREATING THE R VECTOR IN THE pqw COORDINATE FRAME
    R_pqw(1,1) = p*cos(nu(j))/(1 + ecc(j)*cos(nu(j)));
    R_pqw(2,1) = p*sin(nu(j))/(1 + ecc(j)*cos(nu(j)));
    R_pqw(3,1) = 0;
    
    % CREATING THE V VECTOR IN THE pqw COORDINATE FRAME    
    V_pqw(1,1) = -(U/p)^.5*sin(nu(j));
    V_pqw(2,1) =  (U/p)^.5*(ecc(j) + cos(nu(j)));
    V_pqw(3,1) =   0;
    
    % ROTATING THE pqw VECOTRS INTO THE ECI FRAME (ijk)
    R(:,j) = R3(-Omega(j))*R1(-inc(j))*R3(-w(j))*R_pqw;
    V(:,j) = R3(-Omega(j))*R1(-inc(j))*R3(-w(j))*V_pqw;
end