clear all;
addpath('C:\git\GitHub\ValkyrieAI\Matlab');
addpath('C:\git\GitHub\ValkyrieRNN\DroneData\csv');

noisyData = csvread('timeSeriesDataNoisy.csv',1,0);

limit = 60000;

pitch = noisyData(1:limit,2);
roll = noisyData(1:limit,3);
asp = noisyData(1:limit,18);
asr = noisyData(1:limit,19);
asy = noisyData(1:limit,20);

time = linspace(0, limit*0.002, limit)';

asp_sim = [time, asp];

quad_initialization(0);

quadModel.g = 9.81;
quadModel.d = 0.2223;
quadModel.mass = 0.150;
quadModel.ct = 1.4865e-7;
quadModel.cq = 2.9250e-9;
quadModel.Jx = 0.0045;
quadModel.Jy = 0.0045;
quadModel.Jz = 0.0086;
quadModel.Jm = 3.7882e-7;
quadModel.Jb = [quadModel.Jx,0,0; 0,quadModel.Jy,0; 0,0,quadModel.Jz];
quadModel.Jbinv = inv(quadModel.Jb);
quadModel.dctcq = [0,3.30374625e-08,0,-3.30374625e-08;
    -3.30374625e-08,0,3.30374625e-08,0;
    -2.925e-09,2.925e-09,-2.925e-09,2.925e-09];

quadModel.motor_m = 0.0730;
quadModel.motor_dm = 0.2223;
quadModel.motor_h = 0.0318;
quadModel.motor_r = 0.0140;
quadModel.ESC_m = 0.0300;
quadModel.ESC_a = 0.0254;
quadModel.ESC_b = 0.0572;
quadModel.ESC_ds = 0.0826;
quadModel.HUB_m = 0.4310;
quadModel.HUB_r = 0.0325;
quadModel.HUB_H = 0.0429;
quadModel.arms_m = 0.0450;
quadModel.arms_r = 0.0325;
quadModel.arms_L = 0.115;
quadModel.arms_da = 0.0508;
quadModel.T = 0.0560;
quadModel.minThr = 5;
quadModel.cr = 80.5840;
quadModel.b = 976.20;
quadModel.plusConfig = 1;

    
input = struct();
input.startTime = 0.0;
input.endTime = 10.0;
input.axis = 'pitch';
input.stepMagnitude = 10;
input.stepEnable = 0.1;

sim_start_time = 0;
sim_end_time = 15;

% Controller parameters
ANGLE_KP_PHI = 2.2;
ANGLE_KI_PHI = 4.5;
ANGLE_KD_PHI = 1.3;

RATE_KP_PHI = 1.5;
RATE_KI_PHI = 4.0; 
RATE_KD_PHI = 0.05;

ANGLE_CTRL_OUT_LIMIT = 100;     % Max rotation rate desired cmd
RATE_CTRL_OUT_LIMIT = 500;      % Max motor cmd delta
 
%{
figure(1);clf(1);hold on; grid on;
plot(time, pitch);
plot(time, asp);
%}