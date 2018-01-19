% This file is intended to plot the results of a known system simulation
% against that solved by the AI software. In reality, this is only used for
% some intermediary testing and debugging steps, so don't pay much
% attention to this code. 

%% Manual Solve ODE Step Response
a = [-8.202, -2.029; -0.149 -3.25];
b = [1.14; -1.23];
c = [1 0];
d = 0;

x0 = [0; 0];
dt = 0.01;
t0 = 0;
tf = 1.5;
numSteps = tf/dt;
%[t,y] = ode45(@func, [0 10], y0, [], a,b,1);

setPoint = 1;
Yout = [];
Tout = [];

Kp = 12.05;
Ki = 68.38;
Kd = 0.0;

previousError = 0;
integral = 0;
u = 0;
for i=1:numSteps
    dxdt = RK4(@func,t0,x0,dt, a,b,u);
    
    error = setPoint-dxdt(1,:);
    integral = integral+error*dt;
    derivative = (error-previousError)/dt;
    
    u = Kp*error + Ki*integral + Kd*derivative;
    
    % Logging 
    Yout = [Yout x0];
    Tout = [Tout t0];
    
    % Update the variables
    x0 = dxdt;
    t0 = t0+dt;
    previousError = error;
end

figure(2); clf(2);
plot(Tout, Yout(1,:))
title('Step Response: PID Tuner');
ylabel('Amplitude');
xlabel('Time (seconds)');
grid on;

%% DAI Solved Step Response
path = 'C:\LocalTempProjects\ValkyrieAI\ValkyrieAI\x64\Debug\';

%----------OP1-----------
filename = 'rawOutput';
extension = '.csv';
full_path = strcat(path,filename,extension);
array = csvread(full_path);
y1 = array(1,1:end-2);
t1 = array(2,1:end-2);

figure(3); clf(3);
plot(t1, y1, 'Linewidth', 2)
title('Step Response: DAI OP1');
ylabel('Amplitude');
xlabel('Time (seconds)');
grid on;


