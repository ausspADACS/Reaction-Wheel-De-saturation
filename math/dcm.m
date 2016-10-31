function [ matrix ] = dcm( angle_vec )
%Creates a rotation/attitude matrix from 3 rotation angles about the x, y
%and z axis in degrees; requires a vector containing the angles to be
%passed in row form

% Compute individual rotation matrix for each axis
xRot = [1,0,0; 0, cosd(angle_vec(1,1)), -sind(angle_vec(1,1)); 0, sind(angle_vec(1,1)), cosd(angle_vec(1,1))];
yRot = [cosd(angle_vec(1,2)), 0, sind(angle_vec(1,2)); 0, 1, 0; -sind(angle_vec(1,2)), 0, cosd(angle_vec(1,2))];
zRot = [cosd(angle_vec(1,3)), -sind(angle_vec(1,3)), 0; sind(angle_vec(1,3)), cosd(angle_vec(1,3)), 0; 0, 0, 1];
% Rotate around x, y then z axis
matrix = zRot * yRot * xRot;

end

