%Two body Orbit gravity
% Richard Rieber
% December 18, 2006
% rrieber@gmail.com
%
% Revision 8/21/07: Improved help and comments
%                   Added H1 line for lookfor functionality
%
% function [dX] = TwoBody(t,X,U)
%
% Purpose: This function is designed to work best with ODE45.
%          Given a time, t, and the initial position and velocity.
%          it returns the derivatives of the position and velocity
%          due to gravitational acceleration in a two-body system.
%
% Inputs: o t - time in seconds
%         o X - A 1x6 vector containing [R,V] where R is the cartesian
%               position in km of length 3 and V is the cartesian
%               velocity in km/s of length 3.
%         o U - Gravitational constant of body being orbited (km^3/s^2).
%               Default is Earth at 398600.4415 km^3/s^2.  [OPTIONAL]
%
% Outputs: o dX - A 1x6 vector containing [V,A] where V is the cartesian
%                 velocity in km/s of length 3 and A is the cartesian
%                 acceleration in km/s^2 of length 3.
%
% EXAMPLE: [t,X] = ode45(@TwoBody,tspan,Xo,[],U);
%          where tspan is the time, Xo is the initial conditions, [] are 
%          various ODE parameters which are detailed in "help ode45", and
%          U is the gravitational constant.  See "help ode45" for details
%          on the use this function.
%

function [dX] = TwoBody(t,X,U)

if nargin < 2
    error('Not enough inputs.  See help TwoBody.m')
elseif nargin == 2
    U = 398600.4415; %Earth in km^3/s^2
elseif nargin > 3
    error('Too many inputs.  See help TwoBody.m')
end

dX = zeros(1,6);

Accel_sc_E = -U*X(1:3)/norm(X(1:3))^3;

dX = [X(4:6); Accel_sc_E];