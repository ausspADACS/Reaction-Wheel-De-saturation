%Rotation matrix direction cosine matrix
% Richard Rieber
% September 21, 2006
% rrieber@gmail.com
%
% Revision 8/21/07: Deleted unneeded ; in function name.
%                   Added example.
%                   Added H1 line for lookfor functionality
%
% function A = R2(x)
%
% This function creates a rotation matrix about the 2-axis (or the Y-axis)
%
% A = [cos(x)  0      -sin(x);
%      0       1      0;
%      sin(x)  0      cos(x)];
%
% Inputs:  x - rotation angle in radians
% Outputs: A - the rotation matrix about the Y-axis
%
% EXAMPLE:
%
% R2(pi/4) = 
%    0.7071         0   -0.7071
%         0    1.0000         0
%    0.7071         0    0.7071

function A = R2(x)

if nargin > 1
    error('Too many inputs.  See help file')
end

A = [cos(x)  0      -sin(x);
     0       1      0;
     sin(x)  0      cos(x)];