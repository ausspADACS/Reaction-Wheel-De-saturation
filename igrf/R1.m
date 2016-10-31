%Rotation matrix direction cosine matrix
% Richard Rieber
% September 21, 2006
% rrieber@gmail.com
%
% Revision 8/21/07: Deleted unneeded ; in function name.
%                   Added example.
%                   Added H1 line for lookfor functionality
%
% function A = R1(x)
% 
% This function creates a rotation matrix about the 1-axis (or the X-axis)
%
% A = [1      0       0;
%      0      cos(x)  sin(x);
%      0      -sin(x) cos(x)];
%
% Inputs:  x - rotation angle in radians
% Outputs: A - the rotation matrix about the X-axis
%
% EXAMPLE:
%
% R1(pi/4) = 
%     1.0000         0         0
%          0    0.7071    0.7071
%          0   -0.7071    0.7071

function A = R1(x)

if nargin > 1
    error('Too many inputs.  See help file')
end

A = [1      0       0;
     0      cos(x)  sin(x);
     0      -sin(x) cos(x)];