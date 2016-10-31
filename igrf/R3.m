%Rotation matrix direction cosine matrix
% Richard Rieber
% September 21, 2006
% rrieber@gmail.com
%
% Revision 8/21/07: Deleted unneeded ; in function name.
%                   Added example.
%                   Added H1 line for lookfor functionality
%
% function A = R3(x)
%
% This function creates a rotation matrix about the 3-axis (or the Z-axis)
%
% A = [cos(x)  sin(x)     0;
%      -sin(x) cos(x)     0;
%      0       0          1];
%
% Inputs:  x - rotation angle in radians
% Outputs: A - the rotation matrix about the Z-axis
%
% EXAMPLE:
%
% R3(pi/4) = 
%     0.7071    0.7071         0
%    -0.7071    0.7071         0
%          0         0    1.0000

function A = R3(x)

if nargin > 1
    error('Too many inputs.  See help file')
end

A = [cos(x)  sin(x)     0;
     -sin(x) cos(x)     0;
     0       0          1];