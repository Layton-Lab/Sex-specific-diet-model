% clear environment
clear;
clc;

% add helper functions to plot figures
addpath( './helper_functions' );

% run a simulation for a specific diet
t0=0; % start time of the first meal; I recommend to keep as is. 
      % Initial conditions represent fasted individuals after an overnight fast (8-12h)
WR=0; % this work rate, which is 0, as there is no exercise for this simulation
tstart=0; % unused, only for exercise. Start time of exercise
tend=0;   % unused, only for exercise. End time of exercise      
ndays=1;  % number of days of IF (intermittent fasting). Setting ndays=1 means running a single day simulation
lag = 24; % lag=24hrs means simulating a single meal at the start of a 24hr window, 
          % i.e., feeding starts at t=0 and is followed by a 24hrs fast. 
          % For example, lag=8 implies three meals/day (1 meal every 8 hours). 

% OPTION 1 - set amount of carbs and fat directly
carbs=96; fats=33; % amounts in grams per day 
sex = 0; % male (female, 1)
    [Tm1,Ym1,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm]=call_IFMod(t0,sex,WR,carbs,fats,ndays,lag,tstart,tend);
    Tm1=(Tm1-tstart)/60; % Subtracting tstart adjusts for cases where the action (feeding or exercise) 
                         % doesn't start at 0. This adjustment simplifies plotting by displaying negative
                         % times for the hours or minutes prior to the start of the activity, with the 
                         % activity consistently occurring at the marker t=0 on plots.
sex = 1; % female (male,0)
    [Tf1,Yf1,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm]=call_IFMod(t0,sex,WR,carbs,fats,ndays,lag,tstart,tend);
    Tf1=(Tf1-tstart)/60;

% OPTION 2 - call_diet.m for a specific diet modality. The file can be modified to your liking.
% Note that call_diet.m will launch call_IFmod.m directly
    sex = 0; % male (female, 1)
    diet = 'HiC'; % options are {'HiC','HiF','Bal'}
    [T,Y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend] = call_diet(sex,diet);
    % extract data
    my_field=''; % you can specify a field name if you want to name the ouput data accordingly
    driver_getdata; % will take care of extracting all data per organ
