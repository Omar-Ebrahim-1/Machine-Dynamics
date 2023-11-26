% Omar Ebrahim 110076575
clear; clc;

% --------------------------------------------
% --------------------------------------------
% Symbolic solution
% --------------------------------------------
% --------------------------------------------

% --------------------------------------------
% Calculate the velocity of point A
% --------------------------------------------

% Define symbolic variables
syms L R theta omega alpha m_c r_c phi
syms v_a_y a_a_y

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

% Calculate the angular velocity and acceleration of the crank
omega_c = -v_a_y / L * cos(phi);
alpha_c = -(a_a_y + omega^2 * L * sin(phi)) / L * cos(phi);


% --------------------------------------------
% Calculate the accelerations of the crank AC
% --------------------------------------------

% Calculate acceleration vector from a_a + alpha * r - omega^2 * r
a_c = a_a + alpha * r_ca - omega^2 * r_ca;
a_cx = a_c(1);
a_cy = a_c(2);

% Calculate the acceleration of the mass center of the crank AC
a_g_cx = -(omega^2 * R * cos(theta))/2
a_g_cy = -(omega^2 * R * sin(theta))/2

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
omega_b = -v_a_y / L*cos(phi);

% Calculate the angular acceleration of the connecting rod from alpha = a/r
alpha_b = -(a_a_y + omega^2 * L * sin(phi)) / L*cos(phi);

% --------------------------------------------
% Calculate the accelerations of the connecting rod AB
% --------------------------------------------

% Calculate acceleration vector from a_a + alpha * r - omega^2 * r
a_b = a_a + alpha * r_ba - omega^2 * r_ba;
a_bx = a_b(1);
a_by = a_b(2);

% Calculate the acceleration of the mass center of the connecting rod AB
a_g_bx = a_ax + alpha * r_ba(1)/2 - omega^2 * r_ba(1);
a_g_by = a_ay + alpha * r_ba(2)/2 - omega^2 * r_ba(2);

% --------------------------------------------
% Calculate the crank forces
% --------------------------------------------

% Calculate the crank forces
F_cx = m_c * r_c/2 * omega^2 * cos(theta);
F_cy = m_c * r_c/2 * omega^2 * sin(theta);

n = r_c/L;
F_I = 0;
F_II = (2 * m_c * r_c * omega^2 * cos(2 * theta))/n;


% --------------------------------------------
% --------------------------------------------
% Numerical solution
% --------------------------------------------
% --------------------------------------------
% Define constants
omega = 5500*2*pi/60; % rad/s
L = 0.4419; % m
R = 0.1016; % m
m_c = 0.98; % kg
m_r = 0.62; % kg
m_p = 0.78; % kg
x = -0.01:0.0001:0.01; % m
theta1 = omega*x; % rad

% --------------------------------------------------------
% The calculations use the symbolic solutions found above
% --------------------------------------------------------

% Calculate the angle phi
phi = matlabFunction(phi); % rad
phi = phi(L, R, theta1); % rad

% Calculate the velocity of point A
v_ax = matlabFunction(v_ax); % m/s
v_ax = v_ax(R, omega, theta1); % m/s

v_ay = matlabFunction(v_ay); % m/s
v_ay = v_ay(R, omega, theta1); % m/s

% Calculate the acceleration of point A
a_ax = matlabFunction(a_ax); % m/s^2
a_ax = a_ax(R, omega, theta1); % m/s^2

a_ay = matlabFunction(a_ay); % m/s^2
a_ay = a_ay(R, omega, theta1); % m/s^2

% Calculate the angular velocity and acceleration of the crank
omega_c = matlabFunction(omega_c); % rad/s
omega_c = omega_c(L, phi, v_ay); % rad/s

alpha_c = matlabFunction(alpha_c); % rad/s^2
alpha_c = alpha_c(L, a_ay, omega_c, phi); % rad/s^2

% Calculate the angular velocity and acceleration of AB
omega_b = matlabFunction(omega_b) % rad/s
omega_b = omega_b(L, phi, v_ay); % rad/s

alpha_b = matlabFunction(alpha_b) % rad/s^2
alpha_b = alpha_b(L, a_ay, omega_b, phi); % rad/s^2

% Calculate the acceleration of the mass center of the crank
a_g_cx = matlabFunction(a_g_cx) % m/s^2
a_g_cx = a_g_cx(omega, R, theta1); % m/s^2

a_g_cy = matlabFunction(a_g_cy) % m/s^2
a_g_cy = a_g_cy(omega, R, theta1); % m/s^2

% Calculate the acceleration of C
a_cx = @(a_ax, alpha_c, omega_c, L, phi) ...
  a_ax + alpha_c * L * sin(phi) - omega_c^2 * L * cos(phi); % m/s^2
a_cx = a_cx(a_ax, alpha_c, omega_c, L, phi); % m/s^2

% Calculate the acceleration of the mass center of the connecting rod
a_G_cx = @(a_ax, alpha_c, omega_c, L, phi) ...
  a_ax + alpha_c * L/2 * sin(phi) - omega_c^2 * L * cos(phi)
a_G_cx = a_G_cx(a_ax, alpha_c, omega_c, L, phi); % m/s^2

a_G_cy = @(a_ay, alpha_c, omega_c, L, phi) ...
  a_ay + alpha_c * L/2 * cos(phi) - omega_c^2 * L * sin(phi)
a_G_cy = a_G_cy(a_ay, alpha_c, omega_c, L, phi); % m/s^2

% Calculate the acceleration of B
a_bx = @(a_ax, alpha_b, omega_b, L, phi) ...
  a_ax + alpha_b * L * sin(phi) - omega_b^2 * L * cos(phi); % m/s^2
a_bx = a_bx(a_ax, alpha_b, omega_b, L, phi); % m/s^2

% Calculate the acceleration of the mass center of the connecting rod
a_G_bx = @(a_ax, alpha_b, omega_b, L, phi) ...
  a_ax + alpha_b * L/2 * sin(phi) - omega_b^2 * L * cos(phi)
a_G_bx = a_G_bx(a_ax, alpha_b, omega_b, L, phi); % m/s^2

a_G_by = @(a_ay, alpha_b, omega_b, L, phi) ...
  a_ay + alpha_b * L/2 * cos(phi) - omega_b^2 * L * sin(phi)
a_G_by = a_G_by(a_ay, alpha_b, omega_b, L, phi); % m/s^2
