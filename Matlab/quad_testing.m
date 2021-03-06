clear all;
USE_RAW = 0;

if USE_RAW
    addpath('C:\git\GitHub\ValkyrieAI\Matlab\Quadcopter Dynamic Modeling and Simulation\Simulation Files\Simulink Models');
    
    %% Initial Conditions
    IC = struct();

    IC.P = 0;       % Roll Rate   (rad/s)
    IC.Q = 0;       % Pitch Rate  (rad/s)
    IC.R = 0;       % Yaw Rate    (rad/s)
    IC.Phi = 0;     % Roll Angle  (deg)
    IC.The = 0;     % Pitch Angle (deg)
    IC.Psi = 0;     % Yaw Angle   (deg)
    IC.U = 0;       % X Velocity  (m/s)
    IC.V = 0;       % Y Velocity  (m/s)
    IC.W = 0;       % Z Velocity  (m/s)
    IC.X = 0;       % X Position  (m)
    IC.Y = 0;       % Y Position  (m)
    IC.Z = 0;       % Z Position  (m)
    IC.w1 = 4000;   % Motor 1 RPM
    IC.w2 = 4000;   % Motor 2 RPM
    IC.w3 = 4000;   % Motor 3 RPM
    IC.w4 = 4000;   % Motor 4 RPM

    %% Quadrotor Model Settings
    quadModel = struct();

    quadModel.g = 9.81;
    quadModel.d = 0.2223;
    quadModel.mass = 0.450;
    quadModel.ct = 1.4865e-7;
    quadModel.cq = 2.9250e-9;
    quadModel.Jx = 0.0095;
    quadModel.Jy = 0.0095;
    quadModel.Jz = 0.0186;
    quadModel.Jm = 3.7882e-6;
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
    quadModel.T = 0.0760;
    quadModel.minThr = 5;
    quadModel.cr = 80.5840;
    quadModel.b = 976.20;
    quadModel.plusConfig = 1;

    %% PID Settings
    Kp_phi = 2;
    Ki_phi = 1.1;
    Kd_phi = 1.2;

    Kp_theta = 2;
    Ki_theta = 1.1;
    Kd_theta = 1.2;

    Kp_psi = 4;
    Ki_psi = 0.5;
    Kd_psi = 3.5;

    %% Simulation Parameters
    sim_start_time = 0;
    sim_end_time = 10;
    control_frequency_hz = 250;

    phi_step_time = 0.5;                % sec
    phi_init_val = deg2rad(0);          % rad
    phi_final_val = deg2rad(45);        % rad

    theta_step_time = 0.6;
    theta_init_val = deg2rad(0);
    theta_final_val = deg2rad(0);

    psi_step_time = 0.9;
    psi_init_val = deg2rad(0);
    psi_final_val = deg2rad(0);

    sim('AC_Quadcopter_Simulation.slx');

    %% Performance Analysis
    PHI = 4;
    THE = 5;
    PSI = 6;

    phi_step = rad2deg(yout(:,PHI));
    the_step = yout(:,THE);
    psi_step = yout(:,PSI);

    S = stepinfo(phi_step, tout, 45);
end

if ~USE_RAW
    quad_initialization(0);
    
    input = struct();
    input.startTime = 0.0;
    input.endTime = 10.0;
    input.axis = 'pitch';
    input.stepMagnitude = 10;
    input.stepEnable = 0.1;
    
    input.Kp = 81.19;
    input.Ki = 1.15;
    input.Kd = 34.92;
    
    sim_out = quad_simulation(input);
    
    figure(1); clf(1); grid on;
    plot(sim_out.rawData(2,:), sim_out.rawData(1,:));
    
    
end