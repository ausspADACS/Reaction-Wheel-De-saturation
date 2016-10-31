%Orbit range elevation azimuth position ground station site latitude longitude
% Richard Rieber
% October 1, 2006
% rrieber@gmail.com
% 
% function topo = topo(ECEF, lat, long, h, Rp)
% 
% Revision: 10/1/09 - Added an H1 look up line. Also preallocated memory for 
%                     r, az, and el before for-loop to speed things up.
% Purpose:  This function computes range, elevation, and azimuth of an orbiting object
%           given the ECEF coordinates of the object and the latitude, longitude, and
%           altitude of the tracking station.
% 
% NOTE:  This function requires the use of the subfunction 'R2.m' which creates a
%        rotation matrix about the 2-axis (Y-axis) and 'R3.m' which creates a
%        rotation matrix about the 3-axis (Z-axis).
% 
% Inputs: o ECEF - A 3 x n matrix of position vectors in the ECEF frame of the orbiting
%                  object in km
%         o lat  - A vector of length n providing the latitude of the tracking station
%                  in radians
%         o long - A vector of length n providing the longitude of the tracking station
%                  in radians
%         o h    - A vector of length n providing the height of the tracking station
%                  in km
%         o Rp   - Radius of planet in km.  Default is Earth at 6378.1363 km [OPTIONAL]
% 
% Outputs: o topo - A 3 x n matrix consisting of range in km, elevation in radians,
%                   and azimuth in radians

function topo = topo(ECEF, lat, long, h, Rp)

%Checking number of inputs for errors
if nargin < 4 || nargin > 5
    error('Incorrect number of inputs.  See help eci2ecef.')
elseif nargin == 4
    Rp = 6378.1363; %Radius of Earth (km)
end

[b,n] = size(ECEF);

%Checking to see if the ECEF vector has 3 elements
if b ~= 3
    error('ECEF vector must have 3 elements to the vector (X,Y,Z) coordinates')
end

r  = zeros(1,n);
az = zeros(1,n);
el = zeros(1,n);

for j = 1:n  %Iterating thru the number of positions provided by user
    r_site = h(j) + Rp;
    R_site = [r_site*cos(lat(j))*cos(long(j));
              r_site*cos(lat(j))*sin(long(j));
              r_site*sin(lat(j))];
    
    R_sat = ECEF(:,j) - R_site;
    
    R_SEZ = R2(pi/2 - lat(j))*R3(long(j))*R_sat;
    
    r(j)   = norm(R_SEZ);
    sin_el = R_SEZ(3)/r(j);
    cos_el = norm(R_SEZ(1:2))/r(j);
    
    sin_az = R_SEZ(2)/norm(R_SEZ(1:2));
    cos_az = -R_SEZ(1)/norm(R_SEZ(1:2));
    
    el(j) = atan2(sin_el,cos_el);
    az(j) = atan2(sin_az,cos_az);
end

topo = [r, el, az];