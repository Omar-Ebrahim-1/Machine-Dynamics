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
syms v_a_y a_a_x a_a_y

% Define position vector
r_atheta = R * [-sin(theta), cos(theta)];

% Calculate velocity vector from r*omega
v_a = r_atheta .* omega;

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

% Calculate the acceleration of the mass center of the crank AC
a_g_cx = -(omega^2 * R * cos(theta))/2;
a_g_cy = -(omega^2 * R * sin(theta))/2;

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
phi = asin(R/L .* sin(theta1)); % rad

% Calculate the velocity of point A
v_ax = @(R, omega, theta1) ...
  R .* omega .* cos(theta1); % m/s
v_ax = v_ax(R, omega, theta1); % m/s

v_ay = @(R, omega, theta1) ...
  R .* omega .* sin(theta1); % m/s
v_ay = v_ay(R, omega, theta1); % m/s

% Calculate the acceleration of point A
a_ax = @(R, omega, theta1) ...
  -R .* omega.^2 .* cos(theta1); % m/s.^2
a_ax = a_ax(R, omega, theta1); % m/s.^2

a_ay = @(R, omega, theta1) ...
  -R .* omega.^2 .* sin(theta1); % m/s.^2
a_ay = a_ay(R, omega, theta1); % m/s.^2

% Calculate the angular velocity and acceleration of the crank
omega_c = matlabFunction(omega_c); % rad/s
omega_c = omega_c(L, phi, v_ay); % rad/s

alpha_c = matlabFunction(alpha_c); % rad/s^2
alpha_c = alpha_c(L, a_ay, omega_c, phi); % rad/s^2

% Calculate the angular velocity and acceleration of AB
omega_b = matlabFunction(omega_b); % rad/s
omega_b = omega_b(L, R, phi, v_ay); % rad/s

alpha_b = matlabFunction(alpha_b); % rad/s^2
alpha_b = alpha_b(L, R, a_ay, omega_b, phi); % rad/s^2

% Calculate the acceleration of point B
a_bx = @(a_ax, alpha_b, omega_b, L, phi) ...
  a_ax + alpha_b .* L .* sin(phi) - omega_b.^2 .* L .* cos(phi); % m/s.^2
a_bx = a_bx(a_ax, alpha_b, omega_b, L, phi); % m/s.^2

a_by = @(a_ay, alpha_b, omega_b, L, phi) ...
  a_ay + alpha_b .* L .* cos(phi) - omega_b.^2 .* L .* sin(phi); % m/s.^2
a_by = a_by(a_ay, alpha_b, omega_b, L, phi); % m/s.^2

% Calculate the acceleration of the crank
a_cx = @(a_ax, alpha_c, omega_c, L, phi) ...
  a_ax + alpha_c .* L .* sin(phi) - omega_c.^2 .* L .* cos(phi); % m/s.^2
a_cx = a_cx(a_ax, alpha_c, omega_c, L, phi); % m/s.^2

a_cy = @(a_ay, alpha_c, omega_c, L, phi) ...
  a_ay + alpha_c .* L .* cos(phi) - omega_c.^2 .* L .* sin(phi); % m/s.^2
a_cy = a_cy(a_ay, alpha_c, omega_c, L, phi); % m/s.^2

% Calculate the acceleration of the mass center of the crank
a_g_cx = matlabFunction(a_g_cx); % m/s.^2
a_g_cx = a_g_cx(omega, R, theta1); % m/s.^2

a_g_cy = matlabFunction(a_g_cy); % m/s.^2
a_g_cy = a_g_cy(omega, R, theta1); % m/s.^2

% Calculate the acceleration of C
a_cx = @(a_ax, alpha_c, omega_c, L, phi) ...
  a_ax + alpha_c .* L .* sin(phi) - omega_c.^2 .* L .* cos(phi); % m/s.^2
a_cx = a_cx(a_ax, alpha_c, omega_c, L, phi); % m/s.^2

% Calculate the acceleration of the mass center of the connecting rod
a_G_cx = @(a_ax, alpha_c, omega_c, L, phi) ...
  a_ax + alpha_c .* L/2 .* sin(phi) - omega_c.^2 .* L .* cos(phi);
a_G_cx = a_G_cx(a_ax, alpha_c, omega_c, L, phi); % m/s.^2

a_G_cy = @(a_ay, alpha_c, omega_c, L, phi) ...
  a_ay + alpha_c .* L/2 .* cos(phi) - omega_c.^2 .* L .* sin(phi);
a_G_cy = a_G_cy(a_ay, alpha_c, omega_c, L, phi); % m/s.^2

% Calculate the acceleration of B
a_bx = @(a_ax, alpha_b, omega_b, L, phi) ...
  a_ax + alpha_b .* L .* sin(phi) - omega_b.^2 .* L .* cos(phi); % m/s.^2
a_bx = a_bx(a_ax, alpha_b, omega_b, L, phi); % m/s.^2

% Calculate the acceleration of the mass center of the connecting rod
a_G_bx = @(a_ax, alpha_b, omega_b, L, phi) ...
  a_ax + alpha_b .* L/2 .* sin(phi) - omega_b.^2 .* L .* cos(phi);
a_G_bx = a_G_bx(a_ax, alpha_b, omega_b, L, phi); % m/s.^2

a_G_by = @(a_ay, alpha_b, omega_b, L, phi) ...
  a_ay + alpha_b .* L/2 .* cos(phi) - omega_b.^2 .* L .* sin(phi);
a_G_by = a_G_by(a_ay, alpha_b, omega_b, L, phi); % m/s.^2

% Calculate the crank forces
C = m_p .* a_cx;
F_cx = m_r .* a_g_cx - C;
F_cy = m_r .* a_g_cy;

% Calculate the connecting rod forces
B = m_p .* a_bx;
F_bx = m_r .* a_G_bx - B;
F_by = m_r .* a_G_by;

% Calculate the Torque
T = R.*cos(theta1) .* (F_by + F_cy) - R.*sin(theta1) .* (F_bx + F_cx);

% --------------------------------------------
% --------------------------------------------
% Plotting
% --------------------------------------------
% --------------------------------------------
% Plot the velocity of point A
figure(1)
plot(x, v_ax)
hold on
plot(x, v_ay)
hold off
title('Velocity of point A')
xlabel('x (m)')
ylabel('Velocity (m/s)')
legend('v_{ax}', 'v_{ay}')
grid on

% Plot the acceleration of point A
figure(2)
plot(x, a_ax)
hold on
plot(x, a_ay)
hold off
title('Acceleration of point A')
xlabel('x (m)')
ylabel('Acceleration (m/s.^2)')
legend('a_{ax}', 'a_{ay}')
grid on

% Plot the acceleration of point B
figure(3)
plot(x, a_bx)
hold on
plot(x, a_by)
hold off
title('Acceleration of point B')
xlabel('x (m)')
ylabel('Acceleration (m/s.^2)')
legend('a_{bx}', 'a_{by}')
grid on


% Plot the angular velocity and acceleration of the crank
figure(4)
plot(x, omega_c)
hold on
plot(x, alpha_c)
hold off
title('Angular velocity and acceleration of the crank')
xlabel('x (m)')
ylabel('Angular velocity (rad/s)')
legend('omega_c', 'alpha_c')
grid on

% Plot the angular velocity and acceleration of the connecting rod
figure(5)
plot(x, omega_b)
hold on
plot(x, alpha_b)
hold off
title('Angular velocity and acceleration of the connecting rod')
xlabel('x (m)')
ylabel('Angular velocity (rad/s)')
legend('omega_b', 'alpha_b')
grid on

% Plot the acceleration of the crank
figure(6)
plot(x, a_cx)
hold on
plot(x, a_cy)
hold off
title('Acceleration of the crank')
xlabel('x (m)')
ylabel('Acceleration (m/s.^2)')
legend('a_{cx}', 'a_{cy}')
grid on

% Plot the acceleration of the mass center of the crank
figure(7)
plot(x, a_g_cx)
hold on
plot(x, a_g_cy)
hold off
title('Acceleration of the mass center of the crank')
xlabel('x (m)')
ylabel('Acceleration (m/s.^2)')
legend('a_{gx}', 'a_{gy}')
grid on

% Plot the acceleration of the mass center of the connecting rod
figure(8)
plot(x, a_G_cx)
hold on
plot(x, a_G_cy)
hold off
title('Acceleration of the mass center of the connecting rod')
xlabel('x (m)')
ylabel('Acceleration (m/s.^2)')
legend('a_{Gx}', 'a_{Gy}')
grid on

% Plot the crank forces
figure(9)
plot(x, F_cx)
hold on
plot(x, F_cy)
hold off
title('Crank forces')
xlabel('x (m)')
ylabel('Force (N)')
legend('F_{cx}', 'F_{cy}')
grid on

% Plot the connecting rod forces
figure(10)
plot(x, F_bx)
hold on
plot(x, F_by)
hold off
title('Connecting rod forces')
xlabel('x (m)')
ylabel('Force (N)')
legend('F_{bx}', 'F_{by}')
grid on

% Plot the torque
figure(11)
plot(x, T)
title('Torque')
xlabel('x (m)')
ylabel('Torque (N.*m)')
legend('T')
grid on
