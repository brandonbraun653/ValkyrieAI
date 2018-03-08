classdef PIDController < handle
    properties (Access = public)
        dispKp = 0;
        dispKi = 0;
        dispKd = 0;
        
        kp=0;
        ki=0;
        kd=0;
        
        controllerDirection = 1;
        
        lastTime = 0;
        outputSum = 0;
        lastInput = 0;
        SampleTime_mS = 100;
        outMin = 0;
        outMax = 0;
    end
    
    methods (Access = public)
        % Constructor
        function obj = PIDController(KP, KI, KD, DIRECTION)
            obj.SetControllerDirection(DIRECTION)
            obj.SetTunings(KP, KI, KD)
        end
        
        function Initialize(obj, inState, outState)
            obj.outputSum = outState;
            obj.lastInput = inState;
            
            if obj.outputSum > obj.outMax
                obj.outputSum = obj.outMax;
            elseif obj.outputSum < obj.outMin    
                obj.outputSum = obj.outMin;
            end
        end
        
        % Calculates the next output
        function output = compute(obj, input, setpoint)
            % Calculate local working variables 
            error = setpoint - input;
            dInput = (input - obj.lastInput);
            
            obj.outputSum = obj.outputSum + obj.ki * error;
            
            % Force limit the integral error
            if obj.outputSum > obj.outMax
                obj.outputSum = obj.outMax;
            elseif obj.outputSum < obj.outMin    
                obj.outputSum = obj.outMin;
            end
            
            % Add proportional on error
            output = obj.kp * error;
            
            % Compute rest of output 
            output = output + (obj.outputSum - obj.kd * dInput);
            
            % Ensure the output can't exceed bounds
            if output > obj.outMax
                output = obj.outMax;
            elseif output < obj.outMin    
                output = obj.outMin;
            end
            
            obj.lastInput = input;
        end
        
        function SetOutputLimits(obj, min, max)
            if min < max
                obj.outMin = min;
                obj.outMax = max;
            end
        end
        
        function SetTunings(obj, kp, ki, kd)
            obj.dispKp = kp;
            obj.dispKi = ki;
            obj.dispKd = kd;
            
            SampleTimeInSec = obj.SampleTime_mS/1000;
            
            obj.kp = kp;
            obj.ki = ki * SampleTimeInSec;
            obj.kd = kd * SampleTimeInSec;
            
            if obj.controllerDirection == 0
                obj.kp = (0 - obj.kp);
                obj.ki = (0 - obj.ki);
                obj.kd = (0 - obj.kd);
            end
            
        end
        
        function SetControllerDirection(obj, direction)
            if obj.controllerDirection ~= direction
                obj.kp = (0 - obj.kp);
                obj.ki = (0 - obj.ki);
                obj.kd = (0 - obj.kd);
            end
            
            obj.controllerDirection = direction;
        end
        
        function SetSampleTime(obj, sample_time)
            if sample_time > 0
                ratio = sample_time / obj.SampleTime_mS;
               
                obj.ki = obj.ki*ratio;
                obj.kd = obj.kd/ratio;
                obj.SampleTime_mS = sample_time;
            end
        end
    end
    
end