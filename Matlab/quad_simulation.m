function result = quad_simulation(x)
    addpath('C:\git\GitHub\ValkyrieAI\Matlab\Quadcopter Dynamic Modeling and Simulation\Simulation Files\Simulink Models');

    %% PID Settings
    assignin('base', 'Kp_phi', 2);
    assignin('base', 'Ki_phi', 1.1);
    assignin('base', 'Kd_phi', 1.2);
    assignin('base', 'Kp_theta', 2);
    assignin('base', 'Ki_theta', 1);
    assignin('base', 'Kd_theta', 1.2);
    assignin('base', 'Kp_psi', 4);
    assignin('base', 'Ki_psi', 0.5);
    assignin('base', 'Kd_psi', 3.5);


    %% Simulation Parameters
    assignin('base', 'x', x);

    evalin('base', 'sim_start_time = x.startTime;');
    evalin('base', 'sim_end_time = x.endTime;');

    if strcmp(x.axis,'roll')
        evalin('base','phi_step_time = x.stepEnable;');
        evalin('base','phi_init_val = deg2rad(0);');
        evalin('base','phi_final_val = deg2rad(x.stepMagnitude);');

        evalin('base','Kp_phi = x.Kp;');
        evalin('base','Ki_phi = x.Ki;');
        evalin('base','Kd_phi = x.Kd;');

    elseif strcmp(x.axis,'pitch')
        evalin('base','theta_step_time = 0;');
        evalin('base','theta_init_val = deg2rad(0);');
        evalin('base','theta_final_val = deg2rad(x.stepMagnitude);');

        evalin('base','Kp_theta = x.Kp;');
        evalin('base','Ki_theta = x.Ki;');
        evalin('base','Kd_theta = x.Kd;');

    elseif strcmp(x.axis,'yaw')
        evalin('base','psi_step_time = 0;');
        evalin('base','psi_init_val = deg2rad(0);');
        evalin('base','psi_final_val = deg2rad(x.stepMagnitude);');

        evalin('base','Kp_psi = x.Kp;');
        evalin('base','Ki_psi = x.Ki;');
        evalin('base','Kd_psi = x.Kd;');
    end

    % This function is evaluated in the base workspace, but outputs data to
    % the function workspace.
    sim('AC_Quadcopter_Simulation.slx');


    %% Performance Analysis
    result = struct();
    step_data = 0;
    sim_data = 0;
    PHI_SIM_DATA = 4;
    THE_SIM_DATA = 5;
    PSI_SIM_DATA = 6;

    if strcmp(x.axis,'roll')
        phi_step = rad2deg(yout(:,PHI_SIM_DATA));
        step_data = stepinfo(phi_step, tout, x.stepMagnitude);
        sim_data = phi_step';

    elseif strcmp(x.axis,'pitch')
        the_step = rad2deg(yout(:,THE_SIM_DATA));
        step_data = stepinfo(the_step, tout, x.stepMagnitude);
        sim_data = the_step';

    elseif strcmp(x.axis,'yaw')
        psi_step = rad2deg(yout(:,PSI_SIM_DATA));
        step_data = stepinfo(psi_step, tout, x.stepMagnitude);
        sim_data = psi_step';
    end

    result.riseTime = step_data.RiseTime;
    result.settlingTime = step_data.SettlingTime;
    result.percentOvershoot = step_data.Overshoot;
    result.steadyStateError = x.stepMagnitude - sim_data(1,end);
    result.rawData = [sim_data; tout'];



end % quad_simulation()
