% Omar Ebrahim 110076575
clear; clc;

% --------------------------------------------
% Calculate the velocity of point A
% --------------------------------------------

% Define symbolic variables
syms L R theta omega alpha

% Define position vector
r_atheta = R * [-sin(theta), cos(theta)];

% Calculate velocity vector from r*omega
v_a = r_atheta .* omega;
v_ax = v_a(1);
v_ay = v_a(2);

% --------------------------------------------
% Calculate the acceleration of point A
% --------------------------------------------

% Define position vector
r_ao = R * [cos(theta), sin(theta)];

% Calculate acceleration vector from -omega^2 * r
a_a = -omega^2 * r_ao;
a_ax = a_a(1);
a_ay = a_a(2);

% --------------------------------------------
% Calculate the velocities of the crank AC
% --------------------------------------------

% Define position vector
r_ca = L * [-sin(theta), cos(theta)];

% Calculate velocity vector from v_a + r*omega
v_c = v_a + omega * r_ca;
v_cx = v_c(1);
v_cy = v_c(2);

% Calculate the angular velocity of the crank from omega = v/r
omega_c = -v_ay / r_ca(2);

% --------------------------------------------
% Calculate the accelerations of the crank AC
% --------------------------------------------

% Calculate acceleration vector from a_a + alpha * r - omega^2 * r
a_c = a_a + alpha * r_ca - omega^2 * r_ca;
a_cx = a_c(1);
a_cy = a_c(2);

% Calculate the acceleration of the mass center of the crank AC
a_g_cx = a_ax + alpha * r_ca(1)/2 - omega^2 * r_ca(1);
a_g_cy = a_ay + alpha * r_ca(2)/2 - omega^2 * r_ca(2);

% --------------------------------------------
% Calculate the velocities of the connecting rod AB
% --------------------------------------------
% Calculate the angle phi of the connecting rod AB
phi = asin(R/L * sin(theta));

% Define position vector
r_ba = L * [-sin(phi), -cos(phi)];

% Calculate velocity vector from v_a + r*omega
v_b = v_a + omega * r_ba;
v_bx = v_b(1);
v_by = v_b(2);

% Calculate the angular velocity of the connecting rod from omega = v/r
omega_b = -v_ay / L*cos(phi);
