clc, clear all

quat_init = [1 0 0 0];
angle = [0 0 360];
change_quat = mat2quat(dcm(angle));
new_quat = quat_mult2(change_quat/norm(change_quat), quat_init);
disp(new_quat)