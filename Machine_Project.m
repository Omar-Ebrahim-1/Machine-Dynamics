% Omar Ebrahim 110076575
clear; clc;

% --------------------------------------------
% Calculate the velocity and acceleration of A
% --------------------------------------------
% Define symbolic variables
syms R theta omega t

% Define position vector
r_a = R * [-sin(theta), cos(theta)];
% Calculate velocity vector from r*omega
v_a = r_a .* omega;

% Extract components
v_ax = v_a(1)
v_ay = v_a(2)
