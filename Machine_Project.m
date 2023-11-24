% Omar Ebrahim 110076575
clear; clc;

% --------------------------------------------
% Calculate the velocity of point A
% --------------------------------------------
% Define symbolic variables
syms R theta omega t

% Define position vector
r_atheta = R * [-sin(theta), cos(theta)];
% Calculate velocity vector from r*omega
v_a = r_atheta .* omega;

% Extract components
v_ax = v_a(1)
v_ay = v_a(2)

% --------------------------------------------
% Calculate the acceleration of point A
% --------------------------------------------

% Define position vector
r_ao = R * [cos(theta), sin(theta)];

% Calculate acceleration vector from -omega^2 * r
a_a = -omega^2 * r_ao;

% Extract components
a_ax = a_a(1)
a_ay = a_a(2)
