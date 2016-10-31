%Kepler orbital elements ECI Position orbit conversion
% Richard Rieber
% September 26, 2006
% rrieber@gmail.com
%
% Revision: 11/11/08 - Fixed input checks to take either a 1x3 or a 3x1
%           vector for R and V function
%
% [a,ecc,inc,Omega,w,nu,w_true,u_true,lambda_true] = elorb(R,V)
%
% function [a,ecc,inc,Omega,w,nu] = elorb(R,V,U)
%
% function [...] = elorb(R,V,U,tol)
%
% Purpose:  This function converts ECI Cartesian coordinates (R & V)
%           to Kepler orbital elements (a, e, i, O, w, nu).
%           This function is derived from algorithm 9 on pg. 120 of
%           David A. Vallado's Fundamentals of Astrodynamics and 
%           Applications (2nd Edition)
% 
% Inputs:  R:   n x 3 radius vector from the center of the Earth to the
%               spacecraft (ECI coordinates) in km.  n rows should be used
%               for multiple inputs
%          V:   A n x 3 velocity vector of the space craft in km/s.
%               n rows should be used for multiple inputs
%          U:   Gravitational constant of body being orbited (km^3/s^2).  Default is Earth
%               at 398600.4415 km^3/s^2.  [OPTIONAL]
%          tol: A tolerance value for checking equality.  Default is 10^-8.
%               [OPTIONAL]
%
% Outputs: The user has two options for outputs.  If 6 outputs are asked for, the
%          normal 6 Kepler elements will be returned (a,ecc,inc,Omega,w,nu).
%          If 9 outputs are asked for, w_true, u_true, and lambda_true will also 
%          be returned.  These parameters are important for near equatorial elliptical
%          orbits, inclined near circular orbits, and near equatorial, near circular
%          orbits
%
%          a:           Semi-major axis of orbit in km
%          ecc:         Eccentricity of orbit
%          inc:         Inclination of orbit in radians
%          Omega:       Right ascension of the ascending node in radians
%          w:           Argument of perigee in radians
%          nu:          True anomaly in radians
%          w_true:      True longitude of periapse in radians if equatorial
%                       elliptical orbit [OPTIONAL]
%          u_true:      Argument of Latitude in radians if orbit is circular
%                       inclined [OPTIONAL]
%          lambda_true: True longitude in radians if orbit is equatorial
%                       circular [OPTIONAL]
%

function [a,ecc,inc,Omega,w,nu,w_true,u_true,lambda_true] = elorb(R,V,U,tol)

if nargin < 2 || nargin > 4
    error('Incorrect number of inputs:  see help elorb')
elseif nargin == 2
    U   = 398600.4415; %km^3/s^2 Gravitational Constant of Earth
    tol = 10^-8;     %A tolerance for checking equality
elseif nargin == 3
    tol = 10^-8;  %A tolerance for checking equality
end

[x,y] = size(R);
[x2,y2] = size(V);

if y == 1
    R = R';
    [x,y] = size(R);
end

if y2 == 1
    V = V';
end
clear x2
clear y2

if y ~= 3
    error('Radius vector is not the right size:  see help elorb')
end

if size(R) ~= size(V)
    error('Velocity and Radius vectors are not the same size.  Check inputs')
end

ecc         = zeros(1,x);
a           = zeros(1,x);
w           = zeros(1,x);
nu          = zeros(1,x);
Omega       = zeros(1,x);
inc         = zeros(1,x);
w_true      = zeros(1,x);
u_true      = zeros(1,x);
lambda_true = zeros(1,x);

for j = 1:x
	h = cross(R(j,:),V(j,:));  %Specfic angular momentum vector
 	N = cross([0,0,1],h);
	
	%Eccenticity vector
	e_vec = ((norm(V(j,:))^2  - U/norm(R(j,:)))*R(j,:) - dot(R(j,:),V(j,:))*V(j,:))/U;
	ecc(j) = norm(e_vec);  %Magnitude of eccentricty vector
	
	zeta = norm(V(j,:))^2/2 - U/norm(R(j,:));  %Specific mechanical energy of orbit
	
	if (1 - abs(ecc(j))) > tol  %Checking to see if orbit is parabolic
        a(j) = -U/(2*zeta);     %Semi-major axis
    else
        a(j) = inf;
	end
	
	inc(j) = acos(h(3)/norm(h)); %inclination of orbit in radians
    
    Omega(j) = acos(N(1)/norm(N)); %Right ascension of ascending node in radians
    if N(2) < 0  %Checking for quadrant
        Omega(j) = 2*pi - Omega;
    end
    
    w(j) = acos(dot(N,e_vec)/(norm(N)*ecc(j)));  %Argument of perigee in radians
    if e_vec(3) < 0  %Checking for quadrant
        w(j) = 2*pi - w(j);
    end
    
    nu(j) = acos(dot(e_vec,R(j,:))/(ecc(j)*norm(R(j,:))));  %True anomaly in radians
    if dot(R(j,:),V(j,:)) < 0  %Checking for quadrant
        nu(j) = 2*pi - nu(j);
    end

    %%%%%%%%%%%%  Special Cases  %%%%%%%%%%%%%%%
    
    % Elliptical Equatorial
    w_true(j) = acos(e_vec(1)/ecc(j));  %True longitude of periapse
    if  e_vec(2) < 0  %Checking for quadrant
        w_true(j) = 2*pi - w_true(j);
    end
    
    % Circular Inclined
    u_true(j) = acos(dot(N,R(j,:))/(norm(N)*norm(R(j,:))));  %Argument of Latitude
    if R(j,3) < 0  %Checking for quadrant
        u_true(j) = 2*pi - u_true(j);
    end
    
    % Circular Equatorial
    lambda_true(j) = acos(R(j,1)/norm(R(j,:)));  %True Longitude
    if R(j,2) < 0  %Checking for quadrant
        lambda_true(j) = 2*pi - lambda_true(j);
    end
end