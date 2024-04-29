clear; clc;


%% Male HiC
sex = 0; % male (female, 1)

    diet = 'HiC'; % options are {'HiC','HiF','Bal'}
    
    % control data
    sm = 0; % change in muscle mass; expressed as percentage
    sa = 0; % change in adipose mass; expressed as percentage
    [T0,Y0,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend] = call_diet(sex,diet,sm,sa);
    % extract data
    T=T0;
    Y=Y0;
    [~,~,~,~,~,~,F,~,~,~] = cellfun(@(t,y) IFMod(t,y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend,sm,sa),num2cell(T),num2cell(Y,2),'uni',0);   
    times = [3 9 24]; 
    [~, idx] = min( abs( T-times ) ); 
    %metabolic fluxes
    F  = cell2mat(F);
    F10 = F(1:6:end,:);
    F20 = F(2:6:end,:);
    F30 = F(3:6:end,:);
    F40 = F(4:6:end,:);
    F50 = F(5:6:end,:);
    F60 = F(6:6:end,:);
    F0 = {F10 F20 F30 F40 F50 F60};

    %% CASE 1: 5% increase in muscle mass
    sm = 0.05; % change in muscle mass; expressed as percentage
    sa = 0; % change in adipose mass; expressed as percentage
    [T,Y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend] = call_diet(sex,diet,sm,sa);    
    % extract data
    [~,~,~,~,~,~,F,~,~,~] = cellfun(@(t,y) IFMod(t,y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend,sm,sa),num2cell(T),num2cell(Y,2),'uni',0);   
    %metabolic fluxes
    F  = cell2mat(F);
    F1 = F(1:6:end,:);
    F2 = F(2:6:end,:);
    F3 = F(3:6:end,:);
    F4 = F(4:6:end,:);
    F5 = F(5:6:end,:);
    F6 = F(6:6:end,:);
    F = {F1 F2 F3 F4 F5 F6};

    % store output of interest
    cdata1 = cdata(sex,sm,sa,F0,F,Y0,Y);

    %% CASE 2: 5% decrease in muscle mass
    sm = -0.05; % change in muscle mass; expressed as percentage
    sa = 0; % change in adipose mass; expressed as percentage
    [T,Y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend] = call_diet(sex,diet,sm,sa);    
    % extract data
    [~,~,~,~,~,~,F,~,~,~] = cellfun(@(t,y) IFMod(t,y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend,sm,sa),num2cell(T),num2cell(Y,2),'uni',0);   
    %metabolic fluxes
    F  = cell2mat(F);
    F1 = F(1:6:end,:);
    F2 = F(2:6:end,:);
    F3 = F(3:6:end,:);
    F4 = F(4:6:end,:);
    F5 = F(5:6:end,:);
    F6 = F(6:6:end,:);
    F = {F1 F2 F3 F4 F5 F6};

    % store output of interest
    cdata2 = cdata(sex,sm,sa,F0,F,Y0,Y);

    %% CASE 3: 5% increase in fat mass
    sm = 0; % change in muscle mass; expressed as percentage
    sa = 0.05; % change in adipose mass; expressed as percentage
    [T,Y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend] = call_diet(sex,diet,sm,sa);    
    % extract data
    [~,~,~,~,~,~,F,~,~,~] = cellfun(@(t,y) IFMod(t,y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend,sm,sa),num2cell(T),num2cell(Y,2),'uni',0);   
    %metabolic fluxes
    F  = cell2mat(F);
    F1 = F(1:6:end,:);
    F2 = F(2:6:end,:);
    F3 = F(3:6:end,:);
    F4 = F(4:6:end,:);
    F5 = F(5:6:end,:);
    F6 = F(6:6:end,:);
    F = {F1 F2 F3 F4 F5 F6};

    % store output of interest
    cdata3 = cdata(sex,sm,sa,F0,F,Y0,Y);

    %% CASE 4: 5% decrease in fat mass
    sm = 0; % change in muscle mass; expressed as percentage
    sa = -0.05; % change in adipose mass; expressed as percentage
    [T,Y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend] = call_diet(sex,diet,sm,sa);    
    % extract data
    [~,~,~,~,~,~,F,~,~,~] = cellfun(@(t,y) IFMod(t,y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend,sm,sa),num2cell(T),num2cell(Y,2),'uni',0);   
    %metabolic fluxes
    F  = cell2mat(F);
    F1 = F(1:6:end,:);
    F2 = F(2:6:end,:);
    F3 = F(3:6:end,:);
    F4 = F(4:6:end,:);
    F5 = F(5:6:end,:);
    F6 = F(6:6:end,:);
    F = {F1 F2 F3 F4 F5 F6};

    % store output of interest
    cdata4 = cdata(sex,sm,sa,F0,F,Y0,Y);

    male_HiC_3 = [cdata1(idx(1),:); cdata2(idx(1),:); cdata3(idx(1),:); cdata4(idx(1),:)];
    male_HiC_9 = [cdata1(idx(2),:); cdata2(idx(2),:); cdata3(idx(2),:); cdata4(idx(2),:)];
    male_HiC_24 = [cdata1(idx(3),:); cdata2(idx(3),:); cdata3(idx(3),:); cdata4(idx(3),:)];

%% Male HiF
sex = 0; % male (female, 1)

    diet = 'HiF'; % options are {'HiC','HiF','Bal'}
    
    % control data
    sm = 0; % change in muscle mass; expressed as percentage
    sa = 0; % change in adipose mass; expressed as percentage
    [T0,Y0,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend] = call_diet(sex,diet,sm,sa);
    % extract data
    T=T0;
    Y=Y0;
    [~,~,~,~,~,~,F,~,~,~] = cellfun(@(t,y) IFMod(t,y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend,sm,sa),num2cell(T),num2cell(Y,2),'uni',0);   
    times = [3 9 24]; 
    [~, idx] = min( abs( T-times ) ); 
    %metabolic fluxes
    F  = cell2mat(F);
    F10 = F(1:6:end,:);
    F20 = F(2:6:end,:);
    F30 = F(3:6:end,:);
    F40 = F(4:6:end,:);
    F50 = F(5:6:end,:);
    F60 = F(6:6:end,:);
    F0 = {F10 F20 F30 F40 F50 F60};

    %% CASE 1: 5% increase in muscle mass
    sm = 0.05; % change in muscle mass; expressed as percentage
    sa = 0; % change in adipose mass; expressed as percentage
    [T,Y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend] = call_diet(sex,diet,sm,sa);    
    % extract data
    [~,~,~,~,~,~,F,~,~,~] = cellfun(@(t,y) IFMod(t,y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend,sm,sa),num2cell(T),num2cell(Y,2),'uni',0);   
    %metabolic fluxes
    F  = cell2mat(F);
    F1 = F(1:6:end,:);
    F2 = F(2:6:end,:);
    F3 = F(3:6:end,:);
    F4 = F(4:6:end,:);
    F5 = F(5:6:end,:);
    F6 = F(6:6:end,:);
    F = {F1 F2 F3 F4 F5 F6};

    % store output of interest
    cdata1 = cdata(sex,sm,sa,F0,F,Y0,Y);

    %% CASE 2: 5% decrease in muscle mass
    sm = -0.05; % change in muscle mass; expressed as percentage
    sa = 0; % change in adipose mass; expressed as percentage
    [T,Y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend] = call_diet(sex,diet,sm,sa);    
    % extract data
    [~,~,~,~,~,~,F,~,~,~] = cellfun(@(t,y) IFMod(t,y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend,sm,sa),num2cell(T),num2cell(Y,2),'uni',0);   
    %metabolic fluxes
    F  = cell2mat(F);
    F1 = F(1:6:end,:);
    F2 = F(2:6:end,:);
    F3 = F(3:6:end,:);
    F4 = F(4:6:end,:);
    F5 = F(5:6:end,:);
    F6 = F(6:6:end,:);
    F = {F1 F2 F3 F4 F5 F6};

    % store output of interest
    cdata2 = cdata(sex,sm,sa,F0,F,Y0,Y);

    %% CASE 3: 5% increase in fat mass
    sm = 0; % change in muscle mass; expressed as percentage
    sa = 0.05; % change in adipose mass; expressed as percentage
    [T,Y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend] = call_diet(sex,diet,sm,sa);    
    % extract data
    [~,~,~,~,~,~,F,~,~,~] = cellfun(@(t,y) IFMod(t,y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend,sm,sa),num2cell(T),num2cell(Y,2),'uni',0);   
    %metabolic fluxes
    F  = cell2mat(F);
    F1 = F(1:6:end,:);
    F2 = F(2:6:end,:);
    F3 = F(3:6:end,:);
    F4 = F(4:6:end,:);
    F5 = F(5:6:end,:);
    F6 = F(6:6:end,:);
    F = {F1 F2 F3 F4 F5 F6};

    % store output of interest
    cdata3 = cdata(sex,sm,sa,F0,F,Y0,Y);

    %% CASE 4: 5% decrease in fat mass
    sm = 0; % change in muscle mass; expressed as percentage
    sa = -0.05; % change in adipose mass; expressed as percentage
    [T,Y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend] = call_diet(sex,diet,sm,sa);    
    % extract data
    [~,~,~,~,~,~,F,~,~,~] = cellfun(@(t,y) IFMod(t,y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend,sm,sa),num2cell(T),num2cell(Y,2),'uni',0);   
    %metabolic fluxes
    F  = cell2mat(F);
    F1 = F(1:6:end,:);
    F2 = F(2:6:end,:);
    F3 = F(3:6:end,:);
    F4 = F(4:6:end,:);
    F5 = F(5:6:end,:);
    F6 = F(6:6:end,:);
    F = {F1 F2 F3 F4 F5 F6};

    % store output of interest
    cdata4 = cdata(sex,sm,sa,F0,F,Y0,Y);

    male_HiF_3 = [cdata1(idx(1),:); cdata2(idx(1),:); cdata3(idx(1),:); cdata4(idx(1),:)];
    male_HiF_9 = [cdata1(idx(2),:); cdata2(idx(2),:); cdata3(idx(2),:); cdata4(idx(2),:)];
    male_HiF_24 = [cdata1(idx(3),:); cdata2(idx(3),:); cdata3(idx(3),:); cdata4(idx(3),:)];

%% Female HiC
sex = 1; 

    diet = 'HiC'; % options are {'HiC','HiF','Bal'}
    
    % control data
    sm = 0; % change in muscle mass; expressed as percentage
    sa = 0; % change in adipose mass; expressed as percentage
    [T0,Y0,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend] = call_diet(sex,diet,sm,sa);
    % extract data
    T=T0;
    Y=Y0;
    [~,~,~,~,~,~,F,~,~,~] = cellfun(@(t,y) IFMod(t,y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend,sm,sa),num2cell(T),num2cell(Y,2),'uni',0);   
    times = [3 9 24]; 
    [~, idx] = min( abs( T-times ) ); 
    %metabolic fluxes
    F  = cell2mat(F);
    F10 = F(1:6:end,:);
    F20 = F(2:6:end,:);
    F30 = F(3:6:end,:);
    F40 = F(4:6:end,:);
    F50 = F(5:6:end,:);
    F60 = F(6:6:end,:);
    F0 = {F10 F20 F30 F40 F50 F60};

    %% CASE 1: 5% increase in muscle mass
    sm = 0.05; % change in muscle mass; expressed as percentage
    sa = 0; % change in adipose mass; expressed as percentage
    [T,Y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend] = call_diet(sex,diet,sm,sa);    
    % extract data
    [~,~,~,~,~,~,F,~,~,~] = cellfun(@(t,y) IFMod(t,y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend,sm,sa),num2cell(T),num2cell(Y,2),'uni',0);   
    %metabolic fluxes
    F  = cell2mat(F);
    F1 = F(1:6:end,:);
    F2 = F(2:6:end,:);
    F3 = F(3:6:end,:);
    F4 = F(4:6:end,:);
    F5 = F(5:6:end,:);
    F6 = F(6:6:end,:);
    F = {F1 F2 F3 F4 F5 F6};

    % store output of interest
    cdata1 = cdata(sex,sm,sa,F0,F,Y0,Y);

    %% CASE 2: 5% decrease in muscle mass
    sm = -0.05; % change in muscle mass; expressed as percentage
    sa = 0; % change in adipose mass; expressed as percentage
    [T,Y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend] = call_diet(sex,diet,sm,sa);    
    % extract data
    [~,~,~,~,~,~,F,~,~,~] = cellfun(@(t,y) IFMod(t,y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend,sm,sa),num2cell(T),num2cell(Y,2),'uni',0);   
    %metabolic fluxes
    F  = cell2mat(F);
    F1 = F(1:6:end,:);
    F2 = F(2:6:end,:);
    F3 = F(3:6:end,:);
    F4 = F(4:6:end,:);
    F5 = F(5:6:end,:);
    F6 = F(6:6:end,:);
    F = {F1 F2 F3 F4 F5 F6};

    % store output of interest
    cdata2 = cdata(sex,sm,sa,F0,F,Y0,Y);

    %% CASE 3: 5% increase in fat mass
    sm = 0; % change in muscle mass; expressed as percentage
    sa = 0.05; % change in adipose mass; expressed as percentage
    [T,Y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend] = call_diet(sex,diet,sm,sa);    
    % extract data
    [~,~,~,~,~,~,F,~,~,~] = cellfun(@(t,y) IFMod(t,y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend,sm,sa),num2cell(T),num2cell(Y,2),'uni',0);   
    %metabolic fluxes
    F  = cell2mat(F);
    F1 = F(1:6:end,:);
    F2 = F(2:6:end,:);
    F3 = F(3:6:end,:);
    F4 = F(4:6:end,:);
    F5 = F(5:6:end,:);
    F6 = F(6:6:end,:);
    F = {F1 F2 F3 F4 F5 F6};

    % store output of interest
    cdata3 = cdata(sex,sm,sa,F0,F,Y0,Y);

    %% CASE 4: 5% decrease in fat mass
    sm = 0; % change in muscle mass; expressed as percentage
    sa = -0.05; % change in adipose mass; expressed as percentage
    [T,Y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend] = call_diet(sex,diet,sm,sa);    
    % extract data
    [~,~,~,~,~,~,F,~,~,~] = cellfun(@(t,y) IFMod(t,y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend,sm,sa),num2cell(T),num2cell(Y,2),'uni',0);   
    %metabolic fluxes
    F  = cell2mat(F);
    F1 = F(1:6:end,:);
    F2 = F(2:6:end,:);
    F3 = F(3:6:end,:);
    F4 = F(4:6:end,:);
    F5 = F(5:6:end,:);
    F6 = F(6:6:end,:);
    F = {F1 F2 F3 F4 F5 F6};

    % store output of interest
    cdata4 = cdata(sex,sm,sa,F0,F,Y0,Y);

    female_HiC_3 = [cdata1(idx(1),:); cdata2(idx(1),:); cdata3(idx(1),:); cdata4(idx(1),:)];
    female_HiC_9 = [cdata1(idx(2),:); cdata2(idx(2),:); cdata3(idx(2),:); cdata4(idx(2),:)];
    female_HiC_24 = [cdata1(idx(3),:); cdata2(idx(3),:); cdata3(idx(3),:); cdata4(idx(3),:)];

%% Female HiF
sex = 1; 

    diet = 'HiF'; % options are {'HiC','HiF','Bal'}
    
    % control data
    sm = 0; % change in muscle mass; expressed as percentage
    sa = 0; % change in adipose mass; expressed as percentage
    [T0,Y0,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend] = call_diet(sex,diet,sm,sa);
    % extract data
    T=T0;
    Y=Y0;
    [~,~,~,~,~,~,F,~,~,~] = cellfun(@(t,y) IFMod(t,y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend,sm,sa),num2cell(T),num2cell(Y,2),'uni',0);   
    times = [3 9 24]; 
    [~, idx] = min( abs( T-times ) ); 
    %metabolic fluxes
    F  = cell2mat(F);
    F10 = F(1:6:end,:);
    F20 = F(2:6:end,:);
    F30 = F(3:6:end,:);
    F40 = F(4:6:end,:);
    F50 = F(5:6:end,:);
    F60 = F(6:6:end,:);
    F0 = {F10 F20 F30 F40 F50 F60};

    %% CASE 1: 5% increase in muscle mass
    sm = 0.05; % change in muscle mass; expressed as percentage
    sa = 0; % change in adipose mass; expressed as percentage
    [T,Y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend] = call_diet(sex,diet,sm,sa);    
    % extract data
    [~,~,~,~,~,~,F,~,~,~] = cellfun(@(t,y) IFMod(t,y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend,sm,sa),num2cell(T),num2cell(Y,2),'uni',0);   
    %metabolic fluxes
    F  = cell2mat(F);
    F1 = F(1:6:end,:);
    F2 = F(2:6:end,:);
    F3 = F(3:6:end,:);
    F4 = F(4:6:end,:);
    F5 = F(5:6:end,:);
    F6 = F(6:6:end,:);
    F = {F1 F2 F3 F4 F5 F6};

    % store output of interest
    cdata1 = cdata(sex,sm,sa,F0,F,Y0,Y);

    %% CASE 2: 5% decrease in muscle mass
    sm = -0.05; % change in muscle mass; expressed as percentage
    sa = 0; % change in adipose mass; expressed as percentage
    [T,Y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend] = call_diet(sex,diet,sm,sa);    
    % extract data
    [~,~,~,~,~,~,F,~,~,~] = cellfun(@(t,y) IFMod(t,y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend,sm,sa),num2cell(T),num2cell(Y,2),'uni',0);   
    %metabolic fluxes
    F  = cell2mat(F);
    F1 = F(1:6:end,:);
    F2 = F(2:6:end,:);
    F3 = F(3:6:end,:);
    F4 = F(4:6:end,:);
    F5 = F(5:6:end,:);
    F6 = F(6:6:end,:);
    F = {F1 F2 F3 F4 F5 F6};

    % store output of interest
    cdata2 = cdata(sex,sm,sa,F0,F,Y0,Y);

    %% CASE 3: 5% increase in fat mass
    sm = 0; % change in muscle mass; expressed as percentage
    sa = 0.05; % change in adipose mass; expressed as percentage
    [T,Y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend] = call_diet(sex,diet,sm,sa);    
    % extract data
    [~,~,~,~,~,~,F,~,~,~] = cellfun(@(t,y) IFMod(t,y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend,sm,sa),num2cell(T),num2cell(Y,2),'uni',0);   
    %metabolic fluxes
    F  = cell2mat(F);
    F1 = F(1:6:end,:);
    F2 = F(2:6:end,:);
    F3 = F(3:6:end,:);
    F4 = F(4:6:end,:);
    F5 = F(5:6:end,:);
    F6 = F(6:6:end,:);
    F = {F1 F2 F3 F4 F5 F6};

    % store output of interest
    cdata3 = cdata(sex,sm,sa,F0,F,Y0,Y);

    %% CASE 4: 5% decrease in fat mass
    sm = 0; % change in muscle mass; expressed as percentage
    sa = -0.05; % change in adipose mass; expressed as percentage
    [T,Y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend] = call_diet(sex,diet,sm,sa);    
    % extract data
    [~,~,~,~,~,~,F,~,~,~] = cellfun(@(t,y) IFMod(t,y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend,sm,sa),num2cell(T),num2cell(Y,2),'uni',0);   
    %metabolic fluxes
    F  = cell2mat(F);
    F1 = F(1:6:end,:);
    F2 = F(2:6:end,:);
    F3 = F(3:6:end,:);
    F4 = F(4:6:end,:);
    F5 = F(5:6:end,:);
    F6 = F(6:6:end,:);
    F = {F1 F2 F3 F4 F5 F6};

    % store output of interest
    cdata4 = cdata(sex,sm,sa,F0,F,Y0,Y);

    female_HiF_3 = [cdata1(idx(1),:); cdata2(idx(1),:); cdata3(idx(1),:); cdata4(idx(1),:)];
    female_HiF_9 = [cdata1(idx(2),:); cdata2(idx(2),:); cdata3(idx(2),:); cdata4(idx(2),:)];
    female_HiF_24 = [cdata1(idx(3),:); cdata2(idx(3),:); cdata3(idx(3),:); cdata4(idx(3),:)];

%% Functions
function val = cdata(sex,sm,sa,F0,F,Y0,Y)

    if sex == 0 && sm ~= 0
        P = 70*(0.4)*0.8;
        Pnew = 70*(0.4+sm)*0.8;
        deltaP = Pnew-P;
    elseif sex == 1 && sm ~= 0
        P = 58*(0.3)*0.8;
        Pnew = 58*(0.3+sm)*0.8;
        deltaP = Pnew-P;
    elseif sex == 0 && sa ~= 0
        P = 70*(0.16);
        Pnew = 70*(0.16+sa);
        deltaP = Pnew-P;
     elseif sex == 1 && sa ~= 0
        P = 58*(0.295);
        Pnew = 58*(0.295+sa);
        deltaP = Pnew-P;  
    end
    
    Fsens = cell(1,6);
    for i=1:6
        Yold = F0{i};
        Ynew = F{i};
        Fsens{i} = sensitivity(Yold, Ynew, deltaP, P);
    end

val = [Fsens{5}(:,[2,5,7,8]) Fsens{6}(:,19)];
end


function val = sensitivity(Yold, Ynew, deltaP, P) 
    val = ( P ./ Yold ) .* ( ( Ynew - Yold ) ./ abs( deltaP ) );
end

