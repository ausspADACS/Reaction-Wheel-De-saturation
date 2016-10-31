%Orbit gravity plot orbit spherical
% Richard Rieber
% November 16, 2006
% rrieber@gmail.com
% 
% Revision 8/21/07: Added functionality for equatorial and circular orbits
%                   Added H1 line for lookfor functionality
%
% Revision 9/25/07: Added functionality for orbiting bodies other than the
%                   earth.
%
% Revision 9/27/09: Added a string input to customize the plot
%                   Now have the function return h, the handle of the
%                   figure
%
% function h = plotorb(ECEF, V_ECEF, mu, Rbody, fig, s)
% 
% Purpose:  This function provides a visualization of an orbit with respect
%           by plotting the orbit around Earth.  
%
% Inputs:  o ECEF   - The ECEF position vector [1x3] of satellite (km)
%          o V_ECEF - The ECEF velocity vector [1x3] of satellite (km/s)
%          o mu     - Gravitational constant of the body being orbited (km^3/s^2)
%                     Default is Earth at 398600.4415 km^3/s^2 [OPTIONAL]
%          O Rbody  - Radius of the body being orbited (km)
%                   - Default is Earth at 6378.1363 km [OPTIONAL]
%          o fig    - Figure number where this will be plotted. [OPTIONAL]
%          o s      - String for customizing the plot (example: '--b').
%                     See, help plot, for more information [OPTIONAL]
%
% Outputs: o h      - Figure handle
%

function [h] = plotorb(ECEF, V_ECEF, mu, Rbody, fig, s)

if nargin < 2 || nargin > 6
    error('Incorrect number of inputs.  See help plotorb.m')
elseif nargin == 2
    h = figure;
    Rbody = 6378.1363; %Radius of Earth (km)
    mu = 398600.4415;  %Gravitational constant of Earth (km^3/s^2)
    s = 'b';
elseif nargin == 3
    warning('If providing mu, you may want to provide Rbody as well')
    warning('Earth''s radius (6378.1363 km) is being used')
    h = figure;
    Rbody = 6378.1363; %Radius of Earth (km)
    s = 'b';
elseif nargin == 4
    h = figure;
    s = 'b';
elseif nargin == 5
    h = figure(fig);
end

%% Creating sphere for Earth
[X,Y,Z] = sphere(20);
X = X.*Rbody;
Y = Y.*Rbody;
Z = Z.*Rbody;

surf(X,Y,Z)
colormap('white')
hold on

clear X
clear Y
clear Z

if Rbody == 6378.1363
    % Earth coastline - Loading, converting from Lat/Long to 3-D cartesian
    % coords, and plotting

    x = load('Coastline.dat');
    x = x*pi/180;
    R = Rbody*ones(length(x),1);
    [X,Y,Z] = sph2cart(x(:,1),x(:,2),R);

    hold on
    plot3(X,Y,Z,'g')
end

%% Plotting current location of the satellite
plot3(ECEF(1),ECEF(2),ECEF(3),'*r')
hold on

[a,ecc,inc,O,w,nutemp] = elorb(ECEF,V_ECEF,mu);

nu2 = linspace(0,2*pi,100);

% Checks for circular and equatorial orbits
if isnan(O)
    O = 0;
end
if isnan(w)
    w = 0;
end

R = randv(ones(1,length(nu2))*a,ones(1,length(nu2))*ecc,...
          ones(1,length(nu2))*inc,ones(1,length(nu2))*O,...
          ones(1,length(nu2))*w,nu2,mu);

plot3(R(1,:),R(2,:),R(3,:),s)
axis equal
hold off