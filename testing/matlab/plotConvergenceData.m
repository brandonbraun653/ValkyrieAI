% plotConvergenceData.m
%
% This script serves to plot various convergence parameter logs from the AI
% software such as fitness values, MSE, statistical data, etc.
clear all;
clc;

%% Add File Paths
logPath = 'C:\Users\Valkyrie\Desktop\TEMP\AILogs\';
addpath(logPath);

filenames = cell(4,1);
filenames{1} = 'dynamicOptimizer.txt';
filenames{2} = 'static1.txt';
filenames{3} = 'static2.txt';
filenames{4} = 'static3.txt';
filenames{5} = 'pitchVer1Log.txt';
filenames{6} = 'rollVer1Log.txt';
filenames{7} = 'yawVer1Log.txt';

extensions = cell(6,1);
extensions{1} = '_bestPerformers.csv';
extensions{2} = '_avgFit.csv';
extensions{3} = '_avgPOSFit.csv';
extensions{4} = '_avgSSERFit.csv';
extensions{5} = '_avgTSFit.csv';
extensions{6} = '_avgTRFit.csv';



%% Import Data
fileData = cell(size(filenames,1),1);

for i=1:size(fileData,1)
    fileData{i}.bestFit = csvread(strcat(filenames{i}, extensions{1}),0,0); 
    fileData{i}.avgFit = csvread(strcat(filenames{i}, extensions{2}),0,0); 
    fileData{i}.avgPOS = csvread(strcat(filenames{i}, extensions{3}),0,0); 
    fileData{i}.avgSSER = csvread(strcat(filenames{i}, extensions{4}),0,0); 
    fileData{i}.avgTS = csvread(strcat(filenames{i}, extensions{5}),0,0); 
    fileData{i}.avgTR = csvread(strcat(filenames{i}, extensions{6}),0,0); 
end

%% Plot
% FILE: DYNAMIC OPTIMIZER
figure(1); clf(1); hold on; grid on;
plot(fileData{1}.bestFit,'-','Linewidth',1.5);
plot(fileData{1}.avgFit,'--','Linewidth',1.0);
title('Avg. & Top Performance: Adaptive');
xlabel('Generation');
ylabel('Score');
ylim([0 1.0]);
legend('Top','Avg','Location','northwest');
print(strcat(logPath,'DynamicOptimizerPerformance'), '-dpdf');

figure(2); clf(2); hold on; grid on;
plot(fileData{1}.avgPOS,'-','Linewidth',1.0);
plot(fileData{1}.avgSSER,'--','Linewidth',1.0);
plot(fileData{1}.avgTS,':','Linewidth',1.0);
plot(fileData{1}.avgTR,'-.','Linewidth',1.0);
title('Avg. Performance Metric Score: Adaptive');
xlabel('Generation');
ylabel('Score');
ylim([0 1.0]);
legend('Percent Overshoot','Steady State Error', 'Settling Time', 'Rise Time',...
    'Location','northwest');
print(strcat(logPath,'DynamicOptimizerMetrics'), '-dpdf');



% FILE: STATIC 1
figure(3); clf(3); hold on; grid on;
plot(fileData{2}.bestFit,'-','Linewidth',1.5);
plot(fileData{2}.avgFit,'--','Linewidth',1.0);
title('Avg. & Top Performance: Static 1');
xlabel('Generation');
ylabel('Score');
ylim([0 1.0]);
legend('Top','Avg','Location','northwest');
print(strcat(logPath,'Static1Performance'), '-dpdf');

figure(4); clf(4); hold on; grid on;
plot(fileData{2}.avgPOS,'-','Linewidth',1.0);
plot(fileData{2}.avgSSER,'--','Linewidth',1.0);
plot(fileData{2}.avgTS,':','Linewidth',1.0);
plot(fileData{2}.avgTR,'-.','Linewidth',1.0);
title('Avg. Performance Metric Score: Static 1');
xlabel('Generation');
ylabel('Score');
ylim([0 1.0]);
legend('Percent Overshoot','Steady State Error', 'Settling Time', 'Rise Time',...
    'Location','northwest');
print(strcat(logPath,'Static1Metrics'), '-dpdf');


% FILE: STATIC 2
figure(5); clf(5); hold on; grid on;
plot(fileData{3}.bestFit,'-','Linewidth',1.5);
plot(fileData{3}.avgFit,'--','Linewidth',1.0);
title('Avg. & Top Performance: Static 2');
xlabel('Generation');
ylabel('Score');
ylim([0 1.0]);
legend('Top','Avg','Location','northwest');
print(strcat(logPath,'Static2Performance'), '-dpdf');

figure(6); clf(6); hold on; grid on;
plot(fileData{3}.avgPOS,'-','Linewidth',1.0);
plot(fileData{3}.avgSSER,'--','Linewidth',1.0);
plot(fileData{3}.avgTS,':','Linewidth',1.0);
plot(fileData{3}.avgTR,'-.','Linewidth',1.0);
title('Avg. Performance Metric Score: Static 2');
xlabel('Generation');
ylabel('Score');
ylim([0 1.0]);
legend('Percent Overshoot','Steady State Error', 'Settling Time', 'Rise Time',...
    'Location','northwest');
print(strcat(logPath,'Static2Metrics'), '-dpdf');




% FILE: STATIC 3
figure(7); clf(7); hold on; grid on;
plot(fileData{4}.bestFit,'-','Linewidth',1.5);
plot(fileData{4}.avgFit,'--','Linewidth',1.0);
title('Avg. & Top Performance: Static 3');
xlabel('Generation');
ylabel('Score');
ylim([0 1.0]);
legend('Top','Avg','Location','northwest');
print(strcat(logPath,'Static3Performance'), '-dpdf');

figure(8); clf(8); hold on; grid on;
plot(fileData{4}.avgPOS,'-','Linewidth',1.0);
plot(fileData{4}.avgSSER,'--','Linewidth',1.0);
plot(fileData{4}.avgTS,':','Linewidth',1.0);
plot(fileData{4}.avgTR,'-.','Linewidth',1.0);
title('Avg. Performance Metric Score: Static 3');
xlabel('Generation');
ylabel('Score');
ylim([0 1.0]);
legend('Percent Overshoot','Steady State Error', 'Settling Time', 'Rise Time',...
    'Location','northwest');
print(strcat(logPath,'Static3Metrics'), '-dpdf');


% File Pitch Version 1
figure(9); clf(9); hold on; grid on;
plot(fileData{5}.bestFit,'-','Linewidth',1.5);
plot(fileData{5}.avgFit,'--','Linewidth',1.0);
title('Avg. & Top Performance: Pitch');
xlabel('Generation');
ylabel('Score');
ylim([0 1.0]);
legend('Top','Avg','Location','northwest');
print(strcat(logPath,'PitchPerformance'), '-dpdf');

figure(10); clf(10); hold on; grid on;
plot(fileData{5}.avgPOS,'-','Linewidth',1.0);
plot(fileData{5}.avgSSER,'--','Linewidth',1.0);
plot(fileData{5}.avgTS,':','Linewidth',1.0);
plot(fileData{5}.avgTR,'-.','Linewidth',1.0);
title('Avg. Performance Metric Score: Pitch');
xlabel('Generation');
ylabel('Score');
ylim([0 1.0]);
legend('Percent Overshoot','Steady State Error', 'Settling Time', 'Rise Time',...
    'Location','northwest');
print(strcat(logPath,'Static3Metrics'), '-dpdf');


% File Roll Version 1
figure(11); clf(11); hold on; grid on;
plot(fileData{6}.bestFit,'-','Linewidth',1.5);
plot(fileData{6}.avgFit,'--','Linewidth',1.0);
title('Avg. & Top Performance: Roll');
xlabel('Generation');
ylabel('Score');
ylim([0 1.0]);
legend('Top','Avg','Location','northwest');
print(strcat(logPath,'PitchPerformance'), '-dpdf');

figure(12); clf(12); hold on; grid on;
plot(fileData{6}.avgPOS,'-','Linewidth',1.0);
plot(fileData{6}.avgSSER,'--','Linewidth',1.0);
plot(fileData{6}.avgTS,':','Linewidth',1.0);
plot(fileData{6}.avgTR,'-.','Linewidth',1.0);
title('Avg. Performance Metric Score: Roll');
xlabel('Generation');
ylabel('Score');
ylim([0 1.0]);
legend('Percent Overshoot','Steady State Error', 'Settling Time', 'Rise Time',...
    'Location','northwest');
print(strcat(logPath,'Static3Metrics'), '-dpdf');


% File Yaw Version 1
figure(13); clf(13); hold on; grid on;
plot(fileData{7}.bestFit,'-','Linewidth',1.5);
plot(fileData{7}.avgFit,'--','Linewidth',1.0);
title('Avg. & Top Performance: Yaw');
xlabel('Generation');
ylabel('Score');
ylim([0 1.0]);
legend('Top','Avg','Location','northwest');
print(strcat(logPath,'PitchPerformance'), '-dpdf');

figure(14); clf(14); hold on; grid on;
plot(fileData{7}.avgPOS,'-','Linewidth',1.0);
plot(fileData{7}.avgSSER,'--','Linewidth',1.0);
plot(fileData{7}.avgTS,':','Linewidth',1.0);
plot(fileData{7}.avgTR,'-.','Linewidth',1.0);
title('Avg. Performance Metric Score: Yaw');
xlabel('Generation');
ylabel('Score');
ylim([0 1.0]);
legend('Percent Overshoot','Steady State Error', 'Settling Time', 'Rise Time',...
    'Location','northwest');
print(strcat(logPath,'Static3Metrics'), '-dpdf');

