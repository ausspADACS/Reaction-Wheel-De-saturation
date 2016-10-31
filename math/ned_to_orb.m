function dcm = ned_to_orb(orb_inc, orb_period, time_in_orbit)
    omega = 360 / orb_period;
    theta = 90 - orb_inc * cosd(omega * time_in_orbit);
    dcm = [cosd(theta), sind(theta), 0;
           -sind(theta), sind(theta), 0;
           0, 0, 1];
end

