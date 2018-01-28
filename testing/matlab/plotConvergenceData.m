% plotConvergenceData.m
%
% This script serves to plot various convergence parameter logs from the AI
% software such as fitness values, MSE, statistical data, etc.

%% Add File Paths
addpath('DataDumps');

close all;
clear;
clc;

%% Plotting Options



%% Import Data
averageFitness = csvread('avgFit.csv', 0, 0);
fitPOS = csvread('avgPOSFit.csv', 0, 0);
fitSSER = csvread('avgSSERFit.csv', 0, 0);
fitTS = csvread('avgTSFit.csv', 0, 0);
fitTR = csvread('avgTRFit.csv', 0, 0);


%% Plot
figure(1); hold on;
plot(averageFitness);
plot(fitPOS);
plot(fitSSER);
plot(fitTS);
plot(fitTR);