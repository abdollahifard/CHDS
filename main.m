% CHDS: Conflict Handling in Direct Sampling for Stochastic Simulation of Spatial Variables
% This code can be used for both unconditional and conditional conflict handling direct sampling (CHDS)
% Authors: Hesam Soltan Mohammadi, Mohammad Javad Abdollahifard, Faramarz Doulati Ardejani 2019
%% Number of realizations
clear; clc; ReaNum = input('Please enter number of realizations: ');
%% Loading files
load ti_binary
load Is
%% TI
ti = ti_binary; %name of training image
%% Conditioning
%%conditional simulation with specific hard data%%
Is=Is; %the simulation grid (SG) which contains hard conditioning data

%%conditional simulatin with random hard data%%
%J=100; %number of conditioning data
% Is=nan(size(ti));
% ind=randperm(numel(Is));
% ind=ind(1:J);
% Is(ind)=ti(ind);

%%unconditional%%
% Is = nan(size(ti));
%% Other Input Parameters
params.search_radius = 30; %maximum extension of the data events
params.n = 30; %maximum number of points in the data event
params.alpha = 1.2; %value of parameter 'alpha'
params.beta = 1.5; %value of parameter 'beta'
params.m = 4; %devides the SG to m*m windows
params.disp= 0;
params.hr_type = 0;% 0 or 2
params.simul_type=1;% 1- Fast CHDS 2- CHDS 3- DS
%% Main Func.
tic; for i = 1:ReaNum
im = do_simulation_final(Is,ti,params);
Y(:,:,i)= im;
     end; toc