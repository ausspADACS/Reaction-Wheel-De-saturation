%{
    Reaction Wheel De-Saturation Simulation
%}
tic
clc, clear all
addpath('Orbit3D')
addpath('math')
addpath('igrf')
format bank
disp(['Simulation started at ' datestr(now)])
% Inputs
sim_time = 1440; %minutes in simulation
    %95.63 = 1 orbit
time_step = 10; %seconds between data sample--use numbers below 1
    %DO NOT SET BELOW 0.01 LEST YOU WISH TO MEET THE KRAKEN
k = 500000; %De-saturation gain value
%attitude = [0 0 0 0];
B_field_body = zeros(sim_time/time_step,3);
B_field_orb = zeros(sim_time/time_step,3);
wheel_ang_vel = zeros(sim_time/time_step,3);
wheel_ang_vel(1,:) = [837.76 -837.76 837.76];
attitude = [0.71 0 0 0.71];


% Constants
g_parameter = 3.986E14; %graviational parameter of earth, m^3 / s^2
earth_radius = 6378; %radius of earth, km
conv_to_rad = 0.01745329251;
conv_to_deg = 57.29577951;
wheelsMOI = [0.0612 0 0; 0 0.1346 0; 0 0 0.1758];

% Orbit Specs
orb_alt = 550; %mean altitude of orbit, km
inclin = 0; %orbital inclination, degrees
    %0 < inclin < 89 means the satellite rotates with earth
    %91 < inclin < 179 means the satellite rotates against earth
eccent = 0; %eccentricity of orbit, 0 is circular
RA_ascend = 0; %right ascension of ascending node
arg_per   = 0; %argument of perigee
true_anom = 0;  %true anomaly of departure
semi_major = orb_alt + earth_radius; %semi_major axis height, km
orb_vel = sqrt(g_parameter/(semi_major * 1000)); %orbital velocity, km/s
period = (2*pi*(semi_major * 1000)) / orb_vel;
num_orb = (sim_time * 60) / period;
fprintf('Simulating %2.2f Orbits \n', num_orb)
time = datenum(datetime);

% Get latitudes/longitudes over entire orbit
[lat, long] = Orbit3D(RA_ascend, arg_per, true_anom, inclin, semi_major,...
    eccent, time_step, num_orb);

% Call IGRF to create a matrix of magnetic field values in the NED frame
[Bx_NED, By_NED, Bz_NED] = igrf(time, lat, long, orb_alt);
B_field_NED = [Bx_NED, By_NED, Bz_NED];

% Convert magnetic field vectors from NED to orbital frame
c = 1; %Counter for next for loop for matrix indices
for i = 0:(time_step/60):sim_time
    B_field_orb(c,:) = ned_to_orb(inclin, period, i*60)*B_field_NED(c,:)';
    c = c+1;
end

% Main loop
for j = 1:size(B_field_orb, 1)
    % Rotate mag vec from orbital to body, convert from nT to T
    B_field_body(j,:) = (v_rot_q(B_field_orb(j,:)', (attitude/norm(attitude))'))'*10^-9;
    wheel_ang_momentum = (wheelsMOI * wheel_ang_vel(j,:)')';
    mag_moment = k*cross(wheel_ang_momentum, B_field_body(j,:));
    ang_accel = (inv(wheelsMOI)*cross(mag_moment, B_field_body(j,:))'*time_step)';
    % Update wheel angular velocity
    if j < size(B_field_orb, 1)
        wheel_ang_vel(j+1,:) = wheel_ang_vel(j,:) + ang_accel;
    end
end

% Graph Velocity
total_time = (0:(time_step / 60):sim_time)';
figure()
plot(total_time, wheel_ang_vel)
title('Reaction Wheel Angular Velocity')
legend('\omegax',  '\omegay', '\omegaz')
xlabel('Time, minutes')
ylabel('Angular Velocity, deg/s')
grid on

runtime = toc; %end runtime counter
fprintf('Total Simulation Time = %4.2f \n', runtime)