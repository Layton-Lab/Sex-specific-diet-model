function [T,Y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend] = call_diet(sex,diet)
% set up simulation
t0=0; 
WR=0; % NO EXERCISE
tstart=0; % unused, only for exercise
tend=0; % unused, only for exercise
ndays=1; % number of days of to simulate
lag = 24; % single meal

% we are assuming a 400 Kcal meal
cal = 800;
switch diet
    case 'HiC'
        carbs = (cal*0.9)/4; fats = (cal*0.1)/9; % amounts in grams per day 
    case 'HiF' 
        carbs = (cal*0.5)/4; fats = (cal*0.5)/9; % "
    case 'Bal'
        carbs = (cal*0.7)/4; fats = (cal*0.3)/9; % "
end

[T,Y,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm]=call_IFMod(t0,sex,WR,carbs,fats,ndays,lag,tstart,tend);
end