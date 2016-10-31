%Inclination change velocity change orbit hohmann transfer
% Richard Rieber
% October 17, 2006
% rrieber@gmail.com
% 
% function [dV, dI] = dVdI(R_init,R_fin,Inc,U,Tol)
% 
% Purpose:  This function calculates the change of velocity and inclination needed
%           to perform a Hohmann transfer between two circular orbits of different
%           inclinations around earth. The equation to perform this is not
%           algebraically solveable, thus, this is solved with fzero.
%
% Inputs:  o R_init - Radius of initial circular orbit in km
%          o R_fin  - Radius of final circular orbit in km
%          o Inc    - Total change of inclination in radians
%          o U      - Gravitational constant of body being orbited (km^3/s^2).  Default is Earth
%                     at 398600.4415 km^3/s^2.  OPTIONAL
%          o Tol    - A tolerance factor to set iteration limits.  Default is 10^-8
%                     radians - OPTIONAL
%
% Outputs: o dV - A 1x2 vector of the changes of velocity needed at the initial burn
%                 and final burn in km/s
%            dI - A 1x2 vector of the changes of inclination needed at the initial
%                 burn and the final burn in radians
%
% NOTE: This function uses the sub-functions Hohmann.m

function [dV, dI] = dVdI(R_init,R_fin,Inc,U,Tol)

if nargin < 3
    error('Too few inputs.  See help dVdI')
elseif nargin == 3
    U = 398600.4415;  %km^3/s^2 Gravitational Constant of Earth
    Tol = 10^-8;
elseif nargin == 4
    Tol = 10^-8;
elseif nargin > 5
    error('Too many inputs.  See help dVdI')
end

% Initial and Final velocities of circular orbits
v_init = (U/R_init)^.5; %km/s
v_fin  = (U/R_fin)^.5;  %km/s

% Simple Hohmann transfer from initial orbit to final orbit
[dVh,T,a_tx] = Hohmann(R_init,R_fin);
clear T
clear a_tx

% Delta-V's to for the Hohmann transfer
v_tx_a = v_init + dVh(1);
v_tx_b = v_fin  - dVh(2);

opts = optimset('TolX',Tol);

s = fzero(@findS,0,opts,v_init,v_fin,v_tx_a,v_tx_b,Inc);
dI = [s*Inc, (1-s)*Inc];

dV(1) = (v_tx_a^2 + v_init^2 - 2*v_tx_a*v_init*cos(dI(1)))^.5;
dV(2) = (v_tx_b^2 + v_fin^2  - 2*v_tx_b*v_fin *cos(dI(2)))^.5;

end

function x = findS(s,v_init,v_fin,v_tx_a,v_tx_b,Inc)
	% This function finds the best change in inclination with each burn
	
	Va = (v_tx_a^2 + v_init^2 - 2*v_tx_a*v_init*cos(   s *Inc))^.5;
	Vb = (v_tx_b^2 + v_fin^2  - 2*v_tx_b*v_fin *cos((1-s)*Inc))^.5;
	
	x = Va*v_fin*v_tx_b*sin((1-s)*Inc)/(Vb*v_init*v_tx_a) - sin(s*Inc);
end