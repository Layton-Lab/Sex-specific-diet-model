%% Whole-body sex-specific results
clear;
clc;
addpath( './helper_functions' );

%% Whole-body RQ 
% MALE
sex = 0; % male (female, 1)
    diet = 'HiC'; % options are {'HiC','HiF','Bal'}
    [T,Y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend] = call_diet(sex,diet);
    % extract data
    [~,~,~,~,~,uprel,F,~,~,~] = cellfun(@(t,y) IFMod(t,y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend),num2cell(T*60),num2cell(Y,2),'uni',0);
    %uptake-release rates
    uprel = cell2mat(uprel);
    UR1 = uprel(1:7:end,:);
    UR2 = uprel(2:7:end,:);
    UR3 = uprel(3:7:end,:);
    UR4 = uprel(4:7:end,:);
    UR5 = uprel(5:7:end,:);
    UR6 = uprel(6:7:end,:);
    UR7 = uprel(7:7:end,:); % other tissues
    UR = {UR1 UR2 UR3 UR4 UR5 UR6 UR7}; % override UR
    %metabolic fluxes
    F  = cell2mat(F);
    F1 = F(1:6:end,:);
    F2 = F(2:6:end,:);
    F3 = F(3:6:end,:);
    F4 = F(4:6:end,:);
    F5 = F(5:6:end,:);
    F6 = F(6:6:end,:);
    F = {F1 F2 F3 F4 F5 F6}; % override F
    % convert simulation time from minutes to hours 
    Tm_HiC=(T-tstart)/60;
    % compute RER values for each organ
    % RER = CO2 production / O2 consumption; 
    CO2_prod = UR1(:,9)+UR2(:,9)+UR3(:,9)+UR4(:,9)+UR5(:,9)+UR6(:,9)+UR7(:,9);
    O2_util = UR1(:,8)+UR2(:,8)+UR3(:,8)+UR4(:,8)+UR5(:,8)+UR6(:,8)+UR7(:,8);
    RQm_wb_HiC = -CO2_prod./O2_util; % instantaneous values
    % in mmol/hr
    CHOm_oxidrate_HiC = (RQm_wb_HiC - 0.7)/0.3;
    FATm_oxidrate_HiC = 1 - CHOm_oxidrate_HiC;
    % wholebody metabolic fluxes
    glycogenolysis = F1(:,8) + F2(:,8) + F3(:,8) + F4(:,8) + F5(:,8) + F6(:,8); 
    glycogenesis = F1(:,7) + F2(:,7) + F3(:,7) + F4(:,7) + F5(:,7) + F6(:,7);   
    
    gluconeogenesis_II = F5(:,5); % only the liver produces glucose
    glycolysis_II =  F5(:,2)+F1(:,2)+F2(:,2)+F3(:,2)+F4(:,2)+F6(:,2);
    
    lipolysis = F1(:,19) + F2(:,19) + F3(:,19) + F4(:,19) + F5(:,19) + F6(:,19);
    TGsynthesis = F1(:,20) + F2(:,20) + F3(:,20) + F4(:,20) + F5(:,20) + F6(:,20);

    net_gyg_HiC_m = glycogenolysis - glycogenesis;
    net_gng_HiC_m = gluconeogenesis_II/2 - glycolysis_II;
    net_lpl_HiC_m = lipolysis - (TGsynthesis/3); % same as net TG breakdown
    dnl_HiC_m = F1(:,12) + F2(:,12) + F3(:,12) + F4(:,12) + F5(:,12) + F6(:,12); % de novo lipogenesis

 
    diet = 'HiF'; % options are {'HiC','HiF','Bal'}
    [T,Y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend] = call_diet(sex,diet);
    % extract data
    [~,~,~,~,~,uprel,F,~,~,~] = cellfun(@(t,y) IFMod(t,y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend),num2cell(T*60),num2cell(Y,2),'uni',0);
    %uptake-release rates
    uprel = cell2mat(uprel);
    UR1 = uprel(1:7:end,:);
    UR2 = uprel(2:7:end,:);
    UR3 = uprel(3:7:end,:);
    UR4 = uprel(4:7:end,:);
    UR5 = uprel(5:7:end,:);
    UR6 = uprel(6:7:end,:);
    UR7 = uprel(7:7:end,:); % other tissues
    UR = {UR1 UR2 UR3 UR4 UR5 UR6 UR7}; % override UR
    %metabolic fluxes
    F  = cell2mat(F);
    F1 = F(1:6:end,:);
    F2 = F(2:6:end,:);
    F3 = F(3:6:end,:);
    F4 = F(4:6:end,:);
    F5 = F(5:6:end,:);
    F6 = F(6:6:end,:);
    F = {F1 F2 F3 F4 F5 F6}; % override F
    % convert simulation time from minutes to hours 
    Tm_HiF=(T-tstart)/60;
    % compute RER values for each organ
    % RER = CO2 production / O2 consumption; 
    CO2_prod = UR1(:,9)+UR2(:,9)+UR3(:,9)+UR4(:,9)+UR5(:,9)+UR6(:,9)+UR7(:,9);
    O2_util = UR1(:,8)+UR2(:,8)+UR3(:,8)+UR4(:,8)+UR5(:,8)+UR6(:,8)+UR7(:,8);
    RQm_wb_HiF = -CO2_prod./O2_util; % instantaneous values
    % in mmol/hr
    CHOm_oxidrate_HiF = (RQm_wb_HiF - 0.7)/0.3;
    FATm_oxidrate_HiF = 1 - CHOm_oxidrate_HiF;

    % wholebody metabolic fluxes
    glycogenolysis = F1(:,8) + F2(:,8) + F3(:,8) + F4(:,8) + F5(:,8) + F6(:,8); 
    glycogenesis = F1(:,7) + F2(:,7) + F3(:,7) + F4(:,7) + F5(:,7) + F6(:,7);   
    
    gluconeogenesis_II = F5(:,5); % only the liver produces glucose
    glycolysis_II =  F5(:,2)+F1(:,2)+F2(:,2)+F3(:,2)+F4(:,2)+F6(:,2);
    
    lipolysis = F1(:,19) + F2(:,19) + F3(:,19) + F4(:,19) + F5(:,19) + F6(:,19);
    TGsynthesis = F1(:,20) + F2(:,20) + F3(:,20) + F4(:,20) + F5(:,20) + F6(:,20);

    net_gyg_HiF_m = glycogenolysis - glycogenesis;
    net_gng_HiF_m = gluconeogenesis_II/2 - glycolysis_II;
    net_lpl_HiF_m = lipolysis - (TGsynthesis/3);
    dnl_HiF_m = F1(:,12) + F2(:,12) + F3(:,12) + F4(:,12) + F5(:,12) + F6(:,12); % de novo lipogenesis
   

% FEMALE
sex = 1; % female (male, 0)
    diet = 'HiC'; % options are {'HiC','HiF','Bal'}
    [T,Y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend] = call_diet(sex,diet);
    % extract data
    [~,~,~,~,~,uprel,F,~,~,~] = cellfun(@(t,y) IFMod(t,y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend),num2cell(T*60),num2cell(Y,2),'uni',0);
    %uptake-release rates
    uprel = cell2mat(uprel);
    UR1 = uprel(1:7:end,:);
    UR2 = uprel(2:7:end,:);
    UR3 = uprel(3:7:end,:);
    UR4 = uprel(4:7:end,:);
    UR5 = uprel(5:7:end,:);
    UR6 = uprel(6:7:end,:);
    UR7 = uprel(7:7:end,:); % other tissues
    UR = {UR1 UR2 UR3 UR4 UR5 UR6 UR7}; % override UR
    %metabolic fluxes
    F  = cell2mat(F);
    F1 = F(1:6:end,:);
    F2 = F(2:6:end,:);
    F3 = F(3:6:end,:);
    F4 = F(4:6:end,:);
    F5 = F(5:6:end,:);
    F6 = F(6:6:end,:);
    F = {F1 F2 F3 F4 F5 F6}; % override F
    % convert simulation time from minutes to hours 
    Tf_HiC=(T-tstart)/60;
    % compute RER values for each organ
    % RER = CO2 production / O2 consumption; 
    CO2_prod = UR1(:,9)+UR2(:,9)+UR3(:,9)+UR4(:,9)+UR5(:,9)+UR6(:,9)+UR7(:,9);
    O2_util = UR1(:,8)+UR2(:,8)+UR3(:,8)+UR4(:,8)+UR5(:,8)+UR6(:,8)+UR7(:,8);
    RQf_wb_HiC = -CO2_prod./O2_util; % instantaneous values
    % in mmol/hr
    CHOf_oxidrate_HiC = (RQf_wb_HiC - 0.7)/0.3;
    FATf_oxidrate_HiC = 1 - CHOf_oxidrate_HiC;

    % wholebody metabolic fluxes
    glycogenolysis = F1(:,8) + F2(:,8) + F3(:,8) + F4(:,8) + F5(:,8) + F6(:,8); 
    glycogenesis = F1(:,7) + F2(:,7) + F3(:,7) + F4(:,7) + F5(:,7) + F6(:,7);   
    
    gluconeogenesis_II = F5(:,5); % only the liver produces glucose
    glycolysis_II =  F5(:,2)+F1(:,2)+F2(:,2)+F3(:,2)+F4(:,2)+F6(:,2);
    
    lipolysis = F1(:,19) + F2(:,19) + F3(:,19) + F4(:,19) + F5(:,19) + F6(:,19);
    TGsynthesis = F1(:,20) + F2(:,20) + F3(:,20) + F4(:,20) + F5(:,20) + F6(:,20);

    net_gyg_HiC_f = glycogenolysis - glycogenesis;
    net_gng_HiC_f = gluconeogenesis_II/2 - glycolysis_II;
    net_lpl_HiC_f = lipolysis - (TGsynthesis/3);
    dnl_HiC_f = F1(:,12) + F2(:,12) + F3(:,12) + F4(:,12) + F5(:,12) + F6(:,12); % de novo lipogenesis


    diet = 'HiF'; % options are {'HiC','HiF','Bal'}
    [T,Y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend] = call_diet(sex,diet);
    % extract data
    [~,~,~,~,~,uprel,F,~,~,~] = cellfun(@(t,y) IFMod(t,y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend),num2cell(T*60),num2cell(Y,2),'uni',0);
    %uptake-release rates
    uprel = cell2mat(uprel);
    UR1 = uprel(1:7:end,:);
    UR2 = uprel(2:7:end,:);
    UR3 = uprel(3:7:end,:);
    UR4 = uprel(4:7:end,:);
    UR5 = uprel(5:7:end,:);
    UR6 = uprel(6:7:end,:);
    UR7 = uprel(7:7:end,:); % other tissues
    UR = {UR1 UR2 UR3 UR4 UR5 UR6 UR7}; % override UR
    %metabolic fluxes
    F  = cell2mat(F);
    F1 = F(1:6:end,:);
    F2 = F(2:6:end,:);
    F3 = F(3:6:end,:);
    F4 = F(4:6:end,:);
    F5 = F(5:6:end,:);
    F6 = F(6:6:end,:);
    F = {F1 F2 F3 F4 F5 F6}; % override F
    % convert simulation time from minutes to hours 
    Tf_HiF=(T-tstart)/60;
    % compute RER values for each organ
    % RER = CO2 production / O2 consumption; 
    CO2_prod = UR1(:,9)+UR2(:,9)+UR3(:,9)+UR4(:,9)+UR5(:,9)+UR6(:,9)+UR7(:,9);
    O2_util = UR1(:,8)+UR2(:,8)+UR3(:,8)+UR4(:,8)+UR5(:,8)+UR6(:,8)+UR7(:,8);
    RQf_wb_HiF = -CO2_prod./O2_util; % instantaneous values
    % in mmol/hr
    CHOf_oxidrate_HiF = (RQf_wb_HiF - 0.7)/0.3;
    FATf_oxidrate_HiF = 1 - CHOf_oxidrate_HiF;

    % wholebody metabolic fluxes
    glycogenolysis = F1(:,8) + F2(:,8) + F3(:,8) + F4(:,8) + F5(:,8) + F6(:,8); 
    glycogenesis = F1(:,7) + F2(:,7) + F3(:,7) + F4(:,7) + F5(:,7) + F6(:,7);   
    
    gluconeogenesis_II = F5(:,5); % only the liver produces glucose
    glycolysis_II =  F5(:,2)+F1(:,2)+F2(:,2)+F3(:,2)+F4(:,2)+F6(:,2);
    
    lipolysis = F1(:,19) + F2(:,19) + F3(:,19) + F4(:,19) + F5(:,19) + F6(:,19);
    TGsynthesis =  F1(:,20) + F2(:,20) + F3(:,20) + F4(:,20) + F5(:,20) + F6(:,20);

    net_gyg_HiF_f = glycogenolysis - glycogenesis;
    net_gng_HiF_f = gluconeogenesis_II/2 - glycolysis_II;
    net_lpl_HiF_f = lipolysis - (TGsynthesis/3);
    dnl_HiF_f = F1(:,12) + F2(:,12) + F3(:,12) + F4(:,12) + F5(:,12) + F6(:,12); % de novo lipogenesis

%% FIGURES
%--------------------------------------------------------------------------
% FIGURE 5
% Time profile of whole-body respiratory quotient (RQ) in response 
% to a single 800 kcal meal.
%--------------------------------------------------------------------------
figure(1)
N = 2;
C = brewermap(N,'-BrBG');
cols = C;
colororder(cols)
times = [0 0.5 1 1.5 2 3 4 5 6 8 10 12]; % in hours
% male
hold on
[~, idx] = min( abs( Tm_HiC-times ) ); % indices for time points
plot(Tm_HiC(idx),RQm_wb_HiC(idx),'-','LineWidth',2,'Color',cols(1,:))
plot(Tm_HiC(idx),RQm_wb_HiC(idx),'o','MarkerSize',6,'MarkerFaceColor',cols(1,:),'Color',cols(1,:))
[~, idx] = min( abs( Tm_HiF-times ) );
plot(Tm_HiF(idx),RQm_wb_HiF(idx),'-','LineWidth',2,'Color',cols(2,:))
plot(Tm_HiF(idx),RQm_wb_HiF(idx),'o','MarkerSize',6,'MarkerFaceColor',cols(2,:),'Color',cols(2,:))
% female
[~, idx] = min( abs( Tf_HiC-times ) );
plot(Tf_HiC(idx),RQf_wb_HiC(idx),':','LineWidth',2,'Color',cols(1,:))
plot(Tf_HiC(idx),RQf_wb_HiC(idx),'o','MarkerSize',6,'MarkerFaceColor',cols(1,:),'Color',cols(1,:))
[~, idx] = min( abs( Tf_HiF-times ) );
plot(Tf_HiF(idx),RQf_wb_HiF(idx),':','LineWidth',2,'Color',cols(2,:))
plot(Tf_HiF(idx),RQf_wb_HiF(idx),'o','MarkerSize',6,'MarkerFaceColor',cols(2,:),'Color',cols(2,:))
% dummy plots for legend
colororder(cols)
L1 =plot(nan,nan,'o','MarkerSize',6,'MarkerFaceColor',cols(1,:));
L2 = plot(nan,nan,'o','MarkerSize',6,'MarkerFaceColor',cols(2,:));
L3 = plot(nan,nan,'-k','LineWidth',2);
L4 = plot(nan,nan,':k','LineWidth',2);
plot(nan,nan,'o','LineWidth',2)
hold off
xlabel('Time, hr')
ylabel('RQ')
ylim([0.79,0.91])
xlim([0,12])
grid on
box on
legend([L1,L2,L3,L4],{'HiC','HiF','Male','Female'})
set(gca,'FontSize',18,'FontName','Times New Roman')
% exportgraphics(gcf,'wb_rq.png','resolution',600)

%--------------------------------------------------------------------------
% FIGURE 6
% Carbohydrates and fat oxidation fractions in response to 
% a single 800 kcal meal.
%--------------------------------------------------------------------------
fg = figure(2);
N = 2;
cols = brewermap(N,'-BrBG');
% Tile 1
subplot(1,2,1)
% male
hold on
colororder(cols)
plot(Tm_HiC,CHOm_oxidrate_HiC,'-','LineWidth',2)
plot(Tm_HiF,CHOm_oxidrate_HiF,'-','LineWidth',2)
% female
colororder(cols)
plot(Tf_HiC,CHOf_oxidrate_HiC,':','LineWidth',2)
plot(Tf_HiF,CHOf_oxidrate_HiF,':','LineWidth',2)
% dummy plots for legend
colororder(cols)
L1 =plot(nan,nan,'-','LineWidth',4);
L2 = plot(nan,nan,'-','LineWidth',4);
L3 = plot(nan,nan,'-k','LineWidth',2);
L4 = plot(nan,nan,':k','LineWidth',2);
plot(nan,nan,'o','LineWidth',2)
hold off
xlabel('Time, hr')
ylabel('CHO oxidation fraction')
xlim([0,12])
ylim([0.25, 0.75])
grid on
box on
legend([L1,L2,L3,L4],{'HiC','HiF','Male','Female'},'Location','east')
set(gca,'FontSize',18,'FontName','Times New Roman')
% Tile 2
subplot(1,2,2)
% male
hold on
colororder(cols)
plot(Tm_HiC,FATm_oxidrate_HiC,'-','LineWidth',2)
plot(Tm_HiF,FATm_oxidrate_HiF,'-','LineWidth',2)
% female
colororder(cols)
plot(Tf_HiC,FATf_oxidrate_HiC,':','LineWidth',2)
plot(Tf_HiF,FATf_oxidrate_HiF,':','LineWidth',2)
% dummy plots for legend
colororder(cols)
L1 =plot(nan,nan,'-','LineWidth',4);
L2 = plot(nan,nan,'-','LineWidth',4);
L3 = plot(nan,nan,'-k','LineWidth',2);
L4 = plot(nan,nan,':k','LineWidth',2);
plot(nan,nan,'o','LineWidth',2)
hold off
xlabel('Time, hr')
ylabel('Fat oxidation fraction')
xlim([0,12])
ylim([0.25, 0.75])
grid on
box on
legend([L1,L2,L3,L4],{'HiC','HiF','Male','Female'},'Location','east')
set(gca,'FontSize',18,'FontName','Times New Roman')
% add letters to captions
    AddLetters2Plots(fg, {'(a)', '(b)'},...
    'HShift', -0.08, 'VShift', -0.04, 'Direction', 'LeftRight', ...
    'FontSize', 20, 'FontWeight', 'normal')
% export to png
% % exportgraphics(gcf,'ocid-fraction.png','resolution',600)


%% TABLES

%--------------------------------------------------------------------------
% TABLE 3
% Whole body metabolic fluxes
%--------------------------------------------------------------------------
disp('Absorptive phase')
times = [0 6]; 
T             = Tm_HiC;
[~, idx]      = min( abs( T-times ) );
flux_HiC_m    = mean([net_gyg_HiC_m(idx(1):idx(2)), net_gng_HiC_m(idx(1):idx(2)), ...
                net_lpl_HiC_m(idx(1):idx(2)), dnl_HiC_m(idx(1):idx(2))], 1);

T             = Tm_HiF;
[~, idx]      = min( abs( T-times ) );
flux_HiF_m    = mean([net_gyg_HiF_m(idx(1):idx(2)), net_gng_HiF_m(idx(1):idx(2)), ...
                 net_lpl_HiF_m(idx(1):idx(2)), dnl_HiF_m(idx(1):idx(2))], 1);

T             = Tf_HiC;
[~, idx]      = min( abs( T-times ) );
flux_HiC_f    = mean([net_gyg_HiC_f(idx(1):idx(2)), net_gng_HiC_f(idx(1):idx(2)), ...
                net_lpl_HiC_f(idx(1):idx(2)), dnl_HiC_f(idx(1):idx(2))], 1);

T             = Tf_HiF;
[~, idx]      = min( abs( T-times ) );
flux_HiF_f    = mean([net_gyg_HiF_f(idx(1):idx(2)), net_gng_HiF_f(idx(1):idx(2)), ...
                net_lpl_HiF_f(idx(1):idx(2)), dnl_HiF_f(idx(1):idx(2))], 1);

fluxname = ["net glycogenolysis"; "net gluconeogenesis"; "net lipolysis"; "denovo lipogenesis"];

tab = table(fluxname,flux_HiC_m',flux_HiC_f', 100*((flux_HiC_f'./flux_HiC_m')-1));
% Only operate on the numeric variables
table_HiC = varfun(@(x) round(x, 3, 'significant'), tab, 'InputVariables', @isnumeric);
table_HiC.Properties.VariableNames = ["Male", "Female", "\Delta F/M"];
% Add back in the LastName variable
table_HiC = addvars(table_HiC, tab.fluxname, 'Before', 1, 'NewVariableNames', 'Metabolic flux')

tab = table(fluxname,flux_HiF_m',flux_HiF_f', 100*((flux_HiF_f'./flux_HiF_m')-1));
% Only operate on the numeric variables
table_HiF = varfun(@(x) round(x, 3, 'significant'), tab, 'InputVariables', @isnumeric);
table_HiF.Properties.VariableNames = ["Male", "Female", "\Delta F/M"];
% Add back in the LastName variable
table_HiF = addvars(table_HiF, tab.fluxname, 'Before', 1, 'NewVariableNames', 'Metabolic flux')

disp('Postabsorptive phase')
times = [6 12]; 
T             = Tm_HiC;
[~, idx]      = min( abs( T-times ) );
flux_HiC_m    = mean([net_gyg_HiC_m(idx(1):idx(2)), net_gng_HiC_m(idx(1):idx(2)), ...
                net_lpl_HiC_m(idx(1):idx(2)), dnl_HiC_m(idx(1):idx(2))], 1);

T             = Tm_HiF;
[~, idx]      = min( abs( T-times ) );
flux_HiF_m    = mean([net_gyg_HiF_m(idx(1):idx(2)), net_gng_HiF_m(idx(1):idx(2)), ...
                 net_lpl_HiF_m(idx(1):idx(2)), dnl_HiF_m(idx(1):idx(2))], 1);

T             = Tf_HiC;
[~, idx]      = min( abs( T-times ) );
flux_HiC_f    = mean([net_gyg_HiC_f(idx(1):idx(2)), net_gng_HiC_f(idx(1):idx(2)), ...
                net_lpl_HiC_f(idx(1):idx(2)), dnl_HiC_f(idx(1):idx(2))], 1);

T             = Tf_HiF;
[~, idx]      = min( abs( T-times ) );
flux_HiF_f    = mean([net_gyg_HiF_f(idx(1):idx(2)), net_gng_HiF_f(idx(1):idx(2)), ...
                net_lpl_HiF_f(idx(1):idx(2)), dnl_HiF_f(idx(1):idx(2))], 1);

fluxname = ["net glycogenolysis"; "net gluconeogenesis"; "net lipolysis"; "denovo lipogenesis"];

tab = table(fluxname,flux_HiC_m',flux_HiC_f', 100*((flux_HiC_f'./flux_HiC_m')-1));
% Only operate on the numeric variables
table_HiC = varfun(@(x) round(x, 3, 'significant'), tab, 'InputVariables', @isnumeric);
table_HiC.Properties.VariableNames = ["Male", "Female", "\Delta F/M"];
% Add back in the LastName variable
table_HiC = addvars(table_HiC, tab.fluxname, 'Before', 1, 'NewVariableNames', 'Metabolic flux')

tab = table(fluxname,flux_HiF_m',flux_HiF_f', 100*((flux_HiF_f'./flux_HiF_m')-1));
% Only operate on the numeric variables
table_HiF = varfun(@(x) round(x, 3, 'significant'), tab, 'InputVariables', @isnumeric);
table_HiF.Properties.VariableNames = ["Male", "Female", "\Delta F/M"];
% Add back in the LastName variable
table_HiF = addvars(table_HiF, tab.fluxname, 'Before', 1, 'NewVariableNames', 'Metabolic flux')
