clear all;
addpath('model');



% Physical constants.
g = 9.81;
m = 0.05;
L = 0.025;
k = 3e-9;
b = 1e-9;
I = diag([5e-6, 5e-6, 10e-6]);
kd = 0.025;

% Simulation times, in seconds.
tstart = 0;
tend = 50;
dt = 0.1;

ts = tstart:dt:tend;

% Number of points in the simulation.
N = numel(ts);

% Output values, recorded as the simulation runs.
xout = zeros(3, N);
xdotout = zeros(3, N);
thetaout = zeros(3, N);
thetadotout = zeros(3, N);
inputout = zeros(4, N);

% Struct given to the controller. Controller may store its persistent state in it.
controller_params = struct('dt', dt, 'I', I, 'k', k, 'L', L, 'b', b, 'm', m, 'g', g);

% Initial system state.
x = [0; 0; 10];
xdot = zeros(3, 1);
theta = zeros(3, 1);

Kp_ang = 10;
Ki_ang = 2;
Kd_ang = 3;

ang_ctrl = PIDController(Kp_ang, Ki_ang, Kd_ang, 1);
ang_ctrl.Initialize(0, 0);
ang_ctrl.SetSampleTime(dt*1000.0);
ang_ctrl.SetOutputLimits(-500, 500);

% Angular rate controller
rate_ctrl = controller('pid', 10, 3, 5);

% Step Values
roll_step = 1; %deg

thetadot_desired=[0;0;0];
thetadot_actual=[0;0;0];
ind = 0;
for t = ts
    ind = ind + 1;
    
    roll_ang_desired = 0;
    
    if t > 1.0
        roll_ang_desired = roll_step;
    end
    
    roll_rate_desired = ang_ctrl.compute(theta(1), roll_ang_desired);
    
    thetadot_desired(1) = deg2rad(roll_rate_desired);
    
    
    % Get input from controller.
    thetadot_error = thetadot_actual - thetadot_desired;
    [i, controller_params] = rate_ctrl(controller_params, thetadot_error);


    % Compute forces, torques, and accelerations.
    omega = thetadot2omega(thetadot_actual, theta);
    a = acceleration(i, theta, xdot, m, g, k, kd);
    omegadot = angular_acceleration(i, omega, I, L, b, k);

    % Advance system state.
    omega = omega + dt * omegadot;
    thetadot_actual = omega2thetadot(omega, theta); 
    theta = theta + dt * thetadot_actual;
    xdot = xdot + dt * a;
    x = x + dt * xdot;

    % Store simulation state for output.
    xout(:, ind) = x;
    xdotout(:, ind) = xdot;
    thetaout(:, ind) = theta;
    thetadotout(:, ind) = thetadot_actual;
    inputout(:, ind) = i;
end

% Put all simulation variables into an output struct.
result = struct('x', xout, 'theta', thetaout, 'vel', xdotout, ...
                'angvel', thetadotout, 't', ts, 'dt', dt, 'input', inputout);
            
disp('hello')

figure(1); clf(1); grid on; hold on;
plot(result.t, result.angvel(1,:))
plot(result.t, result.angvel(2,:))
plot(result.t, result.angvel(3,:))
legend('Roll', 'Pitch', 'Yaw');
title('Angular Velocities');
xlabel('Time (s)');
ylabel('Rotation Rate (rad/s)');

figure(2); clf(2); grid on; hold on;
plot(result.t, result.theta(1,:))
plot(result.t, result.theta(2,:))
plot(result.t, result.theta(3,:))
legend('Roll', 'Pitch', 'Yaw');
title('Angles');
xlabel('Time (s)');
ylabel('Degrees');
            
%% Local Functions ------------------------------------------          
% Arbitrary test input.
function in = input(t)
    in = zeros(4, 1);
    in(:) = 700;
    in(1) = in(1) + 150;
    in(3) = in(3) + 150;
    in = in .^ 2;
end

% Compute thrust given current inputs and thrust coefficient.
function T = thrust(inputs, k)
    T = [0; 0; k * sum(inputs)];
end

% Compute torques, given current inputs, length, drag coefficient, and thrust coefficient.
function tau = torques(inputs, L, b, k)
    tau = [
        L * k * (inputs(1) - inputs(3))
        L * k * (inputs(2) - inputs(4))
        b * (inputs(1) - inputs(2) + inputs(3) - inputs(4))
    ];
end

% Compute acceleration in inertial reference frame
% Parameters:
%   g: gravity acceleration
%   m: mass of quadcopter
%   k: thrust coefficient
%   kd: global drag coefficient
function a = acceleration(inputs, angles, vels, m, g, k, kd)
    gravity = [0; 0; -g];
    R = rotation(angles);
    T = R * thrust(inputs, k);
    Fd = -kd * vels;
    a = gravity + 1 / m * T + Fd;
end

% Compute angular acceleration in body frame
% Parameters:
%   I: inertia matrix
function omegad = angular_acceleration(inputs, omega, I, L, b, k)
    tau = torques(inputs, L, b, k);
    omegad = inv(I) * (tau - cross(omega, I * omega));
end

% Convert derivatives of roll, pitch, yaw to omega.
function omega = thetadot2omega(thetadot, angles)
    phi = angles(1);
    theta = angles(2);
    psi = angles(3);
    W = [
        1, 0, -sin(theta)
        0, cos(phi), cos(theta)*sin(phi)
        0, -sin(phi), cos(theta)*cos(phi)
    ];
    omega = W * thetadot;
end

% Convert omega to roll, pitch, yaw derivatives
function thetadot = omega2thetadot(omega, angles)
    phi = angles(1);
    theta = angles(2);
    psi = angles(3);
    W = [
        1, 0, -sind(theta)
        0, cosd(phi), cosd(theta)*sind(phi)
        0, -sind(phi), cosd(theta)*cosd(phi)
    ];
    thetadot = inv(W) * omega;
end
