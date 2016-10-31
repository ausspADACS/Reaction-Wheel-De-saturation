% Julian Date
% Richard Rieber
% November 9, 2006
% rrieber@gmail.com
%
% Revision 8/21/07: Added H1 line for lookfor functionality.
%
% function JD = JD(yr,day)
% 
% Purpose:  This function calculates the julian date given
%           the year and the day.
%
% Inputs:  yr  - The year of the given date
%          day - The day of the given date
%
% Outputs: JD  - Julian date

function JD = JD(yr,day)
if nargin ~= 2
    error('Incorrect number of inputs.  See help JD.m')
end

JD = 367*yr + day + 1721013.5;