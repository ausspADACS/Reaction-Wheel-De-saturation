function [ quat ] = mat2quat( matrix )
%Converts a 3x3 matrix to a quaternion in row vector form
quat(1) = 0.5*matrix(1,1)+matrix(2,2)+matrix(3,3);
quat(2) = (matrix(3,2) - matrix(2,3))/(4*quat(1));
quat(3) = (matrix(1,3) - matrix(3,1))/(4*quat(1));
quat(4) = (matrix(2,1) - matrix(1,2))/(4*quat(1));

end

