function result = quad_initialization(config_param)
    % Note: Any use of 'assignin' is to push data to the global workspace
    % so that the simulink model can "see" the variables. Apparently
    % setting the vars manually inside a function does not allow simulink
    % to see it.

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
    assignin('base', 'IC', IC);
    
    % Initialize the control input blocks
    assignin('base', 'phi_step_time',   0);        
    assignin('base', 'phi_init_val',    deg2rad(0));   
    assignin('base', 'phi_final_val',   deg2rad(0));        
    assignin('base', 'theta_step_time', 0);
    assignin('base', 'theta_init_val',  deg2rad(0));
    assignin('base', 'theta_final_val', deg2rad(0));
    assignin('base', 'psi_step_time',   0);
    assignin('base', 'psi_init_val',    deg2rad(0));
    assignin('base', 'psi_final_val',   deg2rad(0));

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
    quadModel.plusConfig = 0;
    assignin('base', 'quadModel', quadModel);
    
    
    %% Simulation Parameters
    assignin('base', 'sim_start_time', 0);
    assignin('base', 'sim_end_time', 0);
   
    % Indicate initialization was a success.
    result = 1;
end