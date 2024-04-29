%% Organ sex-specific results
clear;
clc;
addpath( './helper_functions' );

%% run simulations and store relevant outputs
% MALE
sex = 0; % male (female, 1)
    diet = 'HiC'; % options are {'HiC','HiF','Bal'}
    [T,Y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend] = call_diet(sex,diet);
    % extract data
    [~,~,~,~,~,uprel,F,~,~,~] = cellfun(@(t,y) IFMod(t,y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend),num2cell(T*60),num2cell(Y,2),'uni',0);   
    %uptake-release rates
    uprel = cell2mat(uprel);
    UR1m_HiC = uprel(1:7:end,:);
    UR2m_HiC = uprel(2:7:end,:);
    UR3m_HiC = uprel(3:7:end,:);
    UR4m_HiC = uprel(4:7:end,:);
    UR5m_HiC = uprel(5:7:end,:);
    UR6m_HiC = uprel(6:7:end,:);
    UR7m_HiC = uprel(7:7:end,:); % other tissues
    URm_HiC = {UR1m_HiC UR2m_HiC UR3m_HiC UR4m_HiC UR5m_HiC UR6m_HiC UR7m_HiC}; % override UR
    %metabolic fluxes
    F  = cell2mat(F);
    F1m_HiC = F(1:6:end,:);
    F2m_HiC = F(2:6:end,:);
    F3m_HiC = F(3:6:end,:);
    F4m_HiC = F(4:6:end,:);
    F5m_HiC = F(5:6:end,:);
    F6m_HiC = F(6:6:end,:);
    Fm_HiC = {F1m_HiC F2m_HiC F3m_HiC F4m_HiC F5m_HiC F6m_HiC}; % override F
    % convert simulation time from minutes to hours 
    Tm_HiC=(T-tstart)/60;
    Ym_HiC=Y;
 

    diet = 'HiF'; % options are {'HiC','HiF','Bal'}
    [T,Y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend] = call_diet(sex,diet);
    % extract data
    [~,~,~,~,~,uprel,F,~,~,~] = cellfun(@(t,y) IFMod(t,y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend),num2cell(T*60),num2cell(Y,2),'uni',0);
    %uptake-release rates
    uprel = cell2mat(uprel);
    UR1m_HiF = uprel(1:7:end,:);
    UR2m_HiF = uprel(2:7:end,:);
    UR3m_HiF = uprel(3:7:end,:);
    UR4m_HiF = uprel(4:7:end,:);
    UR5m_HiF = uprel(5:7:end,:);
    UR6m_HiF = uprel(6:7:end,:);
    UR7m_HiF = uprel(7:7:end,:); % other tissues
    URm_HiF = {UR1m_HiF UR2m_HiF UR3m_HiF UR4m_HiF UR5m_HiF UR6m_HiF UR7m_HiF}; % override UR
    %metabolic fluxes
    F  = cell2mat(F);
    F1m_HiF = F(1:6:end,:);
    F2m_HiF = F(2:6:end,:);
    F3m_HiF = F(3:6:end,:);
    F4m_HiF = F(4:6:end,:);
    F5m_HiF = F(5:6:end,:);
    F6m_HiF = F(6:6:end,:);
    Fm_HiF = {F1m_HiF F2m_HiF F3m_HiF F4m_HiF F5m_HiF F6m_HiF}; % override F
    % convert simulation time from minutes to hours 
    Tm_HiF=(T-tstart)/60;
    Ym_HiF=Y;

% FEMALE
sex = 1; % female (male, 0)
    diet = 'HiC'; % options are {'HiC','HiF','Bal'}
    [T,Y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend] = call_diet(sex,diet);
    % extract data
    [~,~,~,~,~,uprel,F,~,~,~] = cellfun(@(t,y) IFMod(t,y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend),num2cell(T*60),num2cell(Y,2),'uni',0);
    %uptake-release rates
    uprel = cell2mat(uprel);
    UR1f_HiC = uprel(1:7:end,:);
    UR2f_HiC = uprel(2:7:end,:);
    UR3f_HiC = uprel(3:7:end,:);
    UR4f_HiC = uprel(4:7:end,:);
    UR5f_HiC = uprel(5:7:end,:);
    UR6f_HiC = uprel(6:7:end,:);
    UR7f_HiC = uprel(7:7:end,:); % other tissues
    URf_HiC = {UR1f_HiC UR2f_HiC UR3f_HiC UR4f_HiC UR5f_HiC UR6f_HiC UR7f_HiC}; % override UR
    %metabolic fluxes
    F  = cell2mat(F);
    F1f_HiC = F(1:6:end,:);
    F2f_HiC = F(2:6:end,:);
    F3f_HiC = F(3:6:end,:);
    F4f_HiC = F(4:6:end,:);
    F5f_HiC = F(5:6:end,:);
    F6f_HiC = F(6:6:end,:);
    Ff_HiC = {F1f_HiC F2f_HiC F3f_HiC F4f_HiC F5f_HiC F6f_HiC}; % override F
    % convert simulation time from minutes to hours 
    Tf_HiC=(T-tstart)/60;
    Yf_HiC=Y;

    diet = 'HiF'; % options are {'HiC','HiF','Bal'}
    [T,Y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend] = call_diet(sex,diet);
    % extract data
    [~,~,~,~,~,uprel,F,~,~,~] = cellfun(@(t,y) IFMod(t,y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend),num2cell(T*60),num2cell(Y,2),'uni',0);
    %uptake-release rates
    uprel = cell2mat(uprel);
    UR1f_HiF = uprel(1:7:end,:);
    UR2f_HiF = uprel(2:7:end,:);
    UR3f_HiF = uprel(3:7:end,:);
    UR4f_HiF = uprel(4:7:end,:);
    UR5f_HiF = uprel(5:7:end,:);
    UR6f_HiF = uprel(6:7:end,:);
    UR7f_HiF = uprel(7:7:end,:); % other tissues
    URf_HiF = {UR1f_HiF UR2f_HiF UR3f_HiF UR4f_HiF UR5f_HiF UR6f_HiF UR7f_HiF}; % override UR
    %metabolic fluxes
    F  = cell2mat(F);
    F1f_HiF = F(1:6:end,:);
    F2f_HiF = F(2:6:end,:);
    F3f_HiF = F(3:6:end,:);
    F4f_HiF = F(4:6:end,:);
    F5f_HiF = F(5:6:end,:);
    F6f_HiF = F(6:6:end,:);
    Ff_HiF = {F1f_HiF F2f_HiF F3f_HiF F4f_HiF F5f_HiF F6f_HiF}; % override F
    % convert simulation time from minutes to hours 
    Tf_HiF=(T-tstart)/60;
    Yf_HiF=Y;

%% Bar graphs for glycogen in liver and muscle
times = [0 6]; % absorptive

T = Tm_HiC;
[~, idx] = min( abs( T-times ) ); 
GLYm_HiC = Ym_HiC(idx(2),[77 121])-Ym_HiC(idx(1),[77 121]);
TGm_HiC  = Ym_HiC(idx(2),[73 117])-Ym_HiC(idx(1),[73 117]);

T = Tm_HiF;
[~, idx] = min( abs( T-times ) ); 
GLYm_HiF = Ym_HiF(idx(2),[77 121])-Ym_HiF(idx(1),[77 121]);
TGm_HiF  = Ym_HiF(idx(2),[73 117])-Ym_HiF(idx(1),[73 117]);

T = Tf_HiC;
[~, idx] = min( abs( T-times ) ); 
GLYf_HiC = Yf_HiC(idx(2),[77 121])-Yf_HiC(idx(1),[77 121]);
TGf_HiC  = Yf_HiC(idx(2),[73 117])-Yf_HiC(idx(1),[73 117]);

T = Tf_HiF;
[~, idx] = min( abs( T-times ) ); 
GLYf_HiF = Yf_HiF(idx(2),[77 121])-Yf_HiF(idx(1),[77 121]);
TGf_HiF  = Yf_HiF(idx(2),[73 117])-Yf_HiF(idx(1),[73 117]);

bars1 = [GLYm_HiC(2) GLYm_HiF(2); GLYf_HiC(2) GLYf_HiF(2)]; % liver
bars2 = [GLYm_HiC(1) GLYm_HiF(1); GLYf_HiC(1) GLYf_HiF(1)]; % muscle
bars3 = [TGm_HiC(2) TGm_HiF(2); TGf_HiC(2) TGf_HiF(2)]; % liver
bars4 = [TGm_HiC(1) TGm_HiF(1); TGf_HiC(1) TGf_HiF(1)]; % muscle

%% FIGURES
%--------------------------------------------------------------------------
% FIGURE 7
% Carbohydrates (glycogen) and fat (TG) storage during 
% the absorptive phase (0--6h).
%--------------------------------------------------------------------------
fg = figure(1);
N = 8;
cols = brewermap(N,'-Spectral');

% Tile 1
subplot(2,2,1)
colororder(cols)
b = bar(bars1); 
b(1).BarWidth = 0.9;
legend('HiC', 'HiF','Location','Northwest','Orientation','horizontal','NumColumns',2)
ylim([0 250])
ylabel('\Delta Glycogen, mM','interpreter','tex')
set(gca,'XTickLabel', {'M' 'F'});
set(gca,'FontSize',18,'FontName','Times New Roman')

% Tile 2
subplot(2,2,2)
colororder(cols)
b = bar(bars2); 
b(1).BarWidth = 0.9;
legend('HiC', 'HiF','Location','Northwest','Orientation','horizontal','NumColumns',2)
ylim([0 17])
ylabel('\Delta Glycogen, mM','interpreter','tex')
set(gca,'XTickLabel', {'M' 'F'});
set(gca,'FontSize',18,'FontName','Times New Roman')

% Tile 3
subplot(2,2,3)
% colororder(cols)
b = bar(bars3); 
b(1).BarWidth = 0.9;
legend('HiC', 'HiF','Location','Northwest','Orientation','horizontal','NumColumns',2)
ylim([0 0.8])
ylabel('\Delta TG, mM','interpreter','tex')
set(gca,'XTickLabel', {'M' 'F'});
set(gca,'FontSize',18,'FontName','Times New Roman')

% Tile 4
subplot(2,2,4)
% colororder(cols)
b = bar(bars4); 
b(1).BarWidth = 0.9;
legend('HiC', 'HiF','Location','Northwest','Orientation','horizontal','NumColumns',2)
ylim([0 8])
ylabel('\Delta TG, mM','interpreter','tex')
set(gca,'XTickLabel', {'M' 'F'});
set(gca,'FontSize',18,'FontName','Times New Roman')
% add letters to captions
    AddLetters2Plots(fg, {'(a)', '(b)', '(c)', '(d)'},...
    'HShift', -0.08, 'VShift', -0.04, 'Direction', 'LeftRight', ...
    'FontSize', 20, 'FontWeight', 'normal')
% export to png
% % exportgraphics(gcf,'Fig7.png','resolution',600)

%--------------------------------------------------------------------------
% FIGURE 8
% Metabolism of carbohydrates during a short-term fast (24h) following 
% a single 800~kcal meal. HiC, high-carbohydrate meal; HiF, high-fat meal.
%--------------------------------------------------------------------------
times = [12 24]; 
j=6; % choose substrate; FFA-->6; 
for i=1:7
    T = Tm_HiF;
    [~, idx_m] = min( abs( T-times ) ); 
    T = Tf_HiF;
    [~, idx_f] = min( abs( T-times ) ); 
    UR_glc_HiF_m(i) = mean( URm_HiF{i}( idx_m(1):idx_m(2) ,j) );
    UR_glc_HiF_f(i) = mean( URf_HiF{i}( idx_f(1):idx_f(2) ,j) );

    T = Tm_HiC;
    [~, idx_m] = min( abs( T-times ) ); 
    T = Tf_HiC;
    [~, idx_f] = min( abs( T-times ) ); 
    UR_glc_HiC_m(i) = mean( URm_HiC{i}( idx_m(1):idx_m(2) ,j) ,1 );
    UR_glc_HiC_f(i) = mean( URf_HiC{i}( idx_f(1):idx_f(2) ,j) ,1 );
end

bars = [UR_glc_HiC_m(1) UR_glc_HiC_f(1) UR_glc_HiF_m(1) UR_glc_HiF_f(1); 
        UR_glc_HiC_m(2) UR_glc_HiC_f(2) UR_glc_HiF_m(2) UR_glc_HiF_f(2);
        UR_glc_HiC_m(3) UR_glc_HiC_f(3) UR_glc_HiF_m(3) UR_glc_HiF_f(3); 
        UR_glc_HiC_m(4) UR_glc_HiC_f(4) UR_glc_HiF_m(4) UR_glc_HiF_f(4);
        UR_glc_HiC_m(5) UR_glc_HiC_f(5) UR_glc_HiF_m(5) UR_glc_HiF_f(5); 
        UR_glc_HiC_m(6) UR_glc_HiC_f(6) UR_glc_HiF_m(6) UR_glc_HiF_f(6);
        UR_glc_HiC_m(7) UR_glc_HiC_f(7) UR_glc_HiF_m(7) UR_glc_HiF_f(7)];

figure(2)
ax1=subplot(2,2,[1,2]);
cols = [0.20,0.53,0.74; 0.49,0.75,0.90; 0.20,0.73,0.74; 0.62,1.00,0.93];
colororder(cols)
b = bar(bars); 
b(1).BarWidth = 1;
hatchfill2(b(2),'single','HatchAngle',45,'HatchDensity',80,'hatchcolor','k');
hatchfill2(b(4),'single','HatchAngle',45,'HatchDensity',80,'hatchcolor','k');
ylabel('Mean uptake/release rate')
% Draw the legend
legendData = {'HiC Male', 'HiC Female', 'HiF Male', 'HiF Female'};
[legend_h, object_h, plot_h, text_str] = legendflex(b, legendData, 'Padding', [2, 2, 10], 'FontSize', 16, 'Location', 'northwest');
% object_h(1) is the first bar's text
% object_h(2) is the second bar's text
% object_h(3) is the third bar's text
% object_h(4) is the fourth bar's text
% object_h(5) is the first bar's patch
% object_h(6) is the second bar's patch
% object_h(7) is the third bar's patch
% object_h(8) is the fourth bar's patch
%
% Set the two patches within the legend
hatchfill2(object_h(6), 'single', 'HatchAngle', 45, 'HatchDensity', 20, 'HatchColor', 'k');
hatchfill2(object_h(8), 'single', 'HatchAngle', 45, 'HatchDensity', 20, 'HatchColor', 'k');
set(gca,'XTickLabel', {'B' 'H' 'M' 'G' 'L' 'A' 'O'});
set(gca,'FontSize',18,'FontName','Times New Roman')

ax2=subplot(2,2,3);
j=1; % substrate number, GLC-->1
C = [0.18,0.52,0.73; 0.20,0.73,0.74]; 
hold on
times = 0:3:24;
T = Tm_HiC;
[~, idx] = min( abs( T-times ) ); 
plot(Tm_HiC(idx), -UR5m_HiC(idx,j),'o','Color',C(1,:),'MarkerSize',8,'markerfacecolor',C(1,:)) % "
plot(Tm_HiC(idx), -UR5m_HiC(idx,j),'-','Color',C(1,:),'LineWidth',3) % "
T = Tf_HiC;
[~, idx] = min( abs( T-times ) ); 
plot(Tf_HiC(idx), -UR5f_HiC(idx,j),'o','Color',C(2,:),'MarkerSize',8,'markerfacecolor',C(2,:)) % "
plot(Tf_HiC(idx), -UR5f_HiC(idx,j),'-','Color',C(2,:),'LineWidth',3) % "
T = Tm_HiF;
[~, idx] = min( abs( T-times ) ); 
plot(Tm_HiF(idx), -UR5m_HiF(idx,j),'d','Color',C(1,:),'MarkerSize',8,'markerfacecolor',C(1,:)) % "
plot(Tm_HiF(idx), -UR5m_HiF(idx,j),'--','Color',C(1,:),'LineWidth',3) % "
T = Tf_HiF;
[~, idx] = min( abs( T-times ) ); 
plot(Tf_HiF(idx), -UR5f_HiF(idx,j),'d','Color',C(2,:),'MarkerSize',8,'markerfacecolor',C(2,:)) % "
plot(Tf_HiF(idx), -UR5f_HiF(idx,j),'--','Color',C(2,:),'LineWidth',3) % "
% dummy plots for legend
L1 = plot(nan,nan,'-o','Color','k','MarkerSize',6,'markerfacecolor','k','LineWidth',2);
L2 =plot(nan,nan,'--d','Color','k','MarkerSize',6,'markerfacecolor','k','LineWidth',2);
% colororder(cols2)
L4 = plot(nan,nan,'-','Color',C(1,:),'LineWidth',3);
L5 = plot(nan,nan,'-','Color',C(2,:),'LineWidth',3);
hold off
xlim([6.1 24])
xticks(6:3:24)
xlabel('Time, hr')
ylabel('Hepatic glucose output, mmol/min')
legend([L1,L2,L4,L5],{'HiC','HiF','Male','Female'},'location','southwest')
set(gca,'FontSize',18,'FontName','Times New Roman')
grid on
box on

ax3=subplot(2,2,4);
j=121; % substrate number, hepatic GLY-->121
C = [0.18,0.52,0.73; 0.20,0.73,0.74]; 
plot(Tm_HiC, Ym_HiC(:,j),'-','Color',C(1,:),'LineWidth',3) % "
hold on
plot(Tf_HiC, Yf_HiC(:,j),'-','Color',C(2,:),'LineWidth',3) % "
plot(Tm_HiF, Ym_HiF(:,j),'--','Color',C(1,:),'LineWidth',3) % "
plot(Tf_HiF, Yf_HiF(:,j),'--','Color',C(2,:),'LineWidth',3) % "
% dummy plots for legend
L1 = plot(nan,nan,'-k','LineWidth',2);
L2 =plot(nan,nan,'--k','LineWidth',2);
L4 = plot(nan,nan,'-','Color',C(1,:),'LineWidth',3);
L5 = plot(nan,nan,'-','Color',C(2,:),'LineWidth',3);
hold off
xlim([0 24])
xticks(0:3:24)
xlabel('Time, hr')
ylabel('Hepatic glycogen, mM')
legend([L1,L2,L4,L5],{'HiC','HiF','Male','Female'},'location','southwest')
set(gca,'FontSize',18,'FontName','Times New Roman')
grid on
box on
% add letters to captions
    AddLetters2Plots({ax1, ax2, ax3}, {'(a)','(b)', '(c)'},...
    'HShift', -0.1, 'VShift', -0.06, 'Direction', 'LeftRight', ...
    'FontSize', 20, 'FontWeight', 'normal')
% export to png
% % exportgraphics(gcf,'Fig8.png','resolution',600)


%--------------------------------------------------------------------------
% FIGURE 11
% Fat metabolism during a short-term fast (24h) following a single 800 
% kcal meal. HiC, high-carbohydrate meal; HiF, high-fat meal.
%--------------------------------------------------------------------------
times = [12 24]; 

j=7; % choose substrate; TG-->7
for i=1:7
    T = Tm_HiF;
    [~, idx_m] = min( abs( T-times ) ); 
    T = Tf_HiF;
    [~, idx_f] = min( abs( T-times ) ); 
    UR_glc_HiF_m(i) = mean( URm_HiF{i}( idx_m(1):idx_m(2) ,j) );
    UR_glc_HiF_f(i) = mean( URf_HiF{i}( idx_f(1):idx_f(2) ,j) );

    T = Tm_HiC;
    [~, idx_m] = min( abs( T-times ) ); 
    T = Tf_HiC;
    [~, idx_f] = min( abs( T-times ) ); 
    UR_glc_HiC_m(i) = mean( URm_HiC{i}( idx_m(1):idx_m(2) ,j) ,1 );
    UR_glc_HiC_f(i) = mean( URf_HiC{i}( idx_f(1):idx_f(2) ,j) ,1 );
end

bars = [UR_glc_HiC_m(1) UR_glc_HiC_f(1) UR_glc_HiF_m(1) UR_glc_HiF_f(1); 
        UR_glc_HiC_m(2) UR_glc_HiC_f(2) UR_glc_HiF_m(2) UR_glc_HiF_f(2);
        UR_glc_HiC_m(3) UR_glc_HiC_f(3) UR_glc_HiF_m(3) UR_glc_HiF_f(3); 
        UR_glc_HiC_m(4) UR_glc_HiC_f(4) UR_glc_HiF_m(4) UR_glc_HiF_f(4);
        UR_glc_HiC_m(5) UR_glc_HiC_f(5) UR_glc_HiF_m(5) UR_glc_HiF_f(5); 
        UR_glc_HiC_m(6) UR_glc_HiC_f(6) UR_glc_HiF_m(6) UR_glc_HiF_f(6);
        UR_glc_HiC_m(7) UR_glc_HiC_f(7) UR_glc_HiF_m(7) UR_glc_HiF_f(7)];

figure(4)
ax3=subplot(2,2,[3,4]);
cols = [0.20,0.53,0.74; 0.49,0.75,0.90; 0.20,0.73,0.74; 0.62,1.00,0.93];
colororder(cols)
b = bar(bars); 
b(1).BarWidth = 1;
hatchfill2(b(2),'single','HatchAngle',45,'HatchDensity',80,'hatchcolor','k');
hatchfill2(b(4),'single','HatchAngle',45,'HatchDensity',80,'hatchcolor','k');
ylabel('Mean uptake/release rate')
% Draw the legend
legendData = {'HiC Male', 'HiC Female', 'HiF Male', 'HiF Female'};
[legend_h, object_h, plot_h, text_str] = legendflex(b, legendData, 'Padding', [2, 2, 10], 'FontSize', 16, 'Location', 'northeast');
% object_h(1) is the first bar's text
% object_h(2) is the second bar's text
% object_h(3) is the third bar's text
% object_h(4) is the fourth bar's text
% object_h(5) is the first bar's patch
% object_h(6) is the second bar's patch
% object_h(7) is the third bar's patch
% object_h(8) is the fourth bar's patch
%
% Set the two patches within the legend
hatchfill2(object_h(6), 'single', 'HatchAngle', 45, 'HatchDensity', 20, 'HatchColor', 'k');
hatchfill2(object_h(8), 'single', 'HatchAngle', 45, 'HatchDensity', 20, 'HatchColor', 'k');
set(gca,'XTickLabel', {'B' 'H' 'M' 'G' 'L' 'A' 'O'});
set(gca,'FontSize',18,'FontName','Times New Roman')
% exportgraphics(gcf,'UR-TG.png','resolution',600)

ax1=subplot(2,2,1);
j=139; % Adipose TG
C = [0.18,0.52,0.73; 0.20,0.73,0.74]; 
hold on
plot(Tm_HiC, Ym_HiC(:,j),'-','Color',C(1,:),'LineWidth',3) % "
plot(Tf_HiC, Yf_HiC(:,j),'-','Color',C(2,:),'LineWidth',3) % "
plot(Tm_HiF, Ym_HiF(:,j),'--','Color',C(1,:),'LineWidth',3) % "
plot(Tf_HiF, Yf_HiF(:,j),'--','Color',C(2,:),'LineWidth',3) % "
% dummy plots for legend
L1 = plot(nan,nan,'-k','LineWidth',2);
L2 =plot(nan,nan,'--k','LineWidth',2);
L4 = plot(nan,nan,'-','Color',C(1,:),'LineWidth',3);
L5 = plot(nan,nan,'-','Color',C(2,:),'LineWidth',3);
hold off
xlim([0 24])
xticks(0:3:24)
xlabel('Time, hr')
ylabel('Adipose TG, mM')
legend([L1,L2,L4,L5],{'HiC','HiF','Male','Female'},'location','southwest')
set(gca,'FontSize',18,'FontName','Times New Roman')
grid on
box on

ax2=subplot(2,2,2);
j=19; % Adipose lipolysis
k = 20; % Adipose TG synthesis
C = [0.18,0.52,0.73; 0.20,0.73,0.74]; 
plot(Tm_HiC, F6m_HiC(:,j)-F6m_HiC(:,k)/3,'-','Color',C(1,:),'LineWidth',3) % "
hold on
plot(Tf_HiC, F6f_HiC(:,j)-F6f_HiC(:,k)/3,'-','Color',C(2,:),'LineWidth',3) % "
plot(Tm_HiF, F6m_HiF(:,j)-F6m_HiF(:,k)/3,'--','Color',C(1,:),'LineWidth',3) % "
plot(Tf_HiF, F6f_HiF(:,j)-F6f_HiF(:,k)/3,'--','Color',C(2,:),'LineWidth',3) % "
% dummy plots for legend
L1 = plot(nan,nan,'-k','LineWidth',2);
L2 =plot(nan,nan,'--k','LineWidth',2);
L4 = plot(nan,nan,'-','Color',C(1,:),'LineWidth',3);
L5 = plot(nan,nan,'-','Color',C(2,:),'LineWidth',3);
hold off
xlim([0 24])
xticks(0:3:24)
xlabel('Time, hr')
ylabel('Adipose net lipolysis mmol/min')
legend([L1,L2,L4,L5],{'HiC','HiF','Male','Female'},'location','southeast')
set(gca,'FontSize',18,'FontName','Times New Roman')
grid on
box on
% add letters to captions
    AddLetters2Plots({ax1,ax2,ax3}, {'(a)', '(b)','(c)'},...
    'HShift', -0.08, 'VShift', -0.04, 'Direction', 'LeftRight', ...
    'FontSize', 20, 'FontWeight', 'normal')

%--------------------------------------------------------------------------
% FIGURE 9
% Change in hepatic energy metabolism with fasting. Values represent 
% averages over the last 12 hours of a 24-hour fast following a 
% single 800~kcal meal.
%--------------------------------------------------------------------------
times = [12 24]; % postabsorptive

% GLYCEROL 
% Note: these are the values reported as glycerol uptake rates in Fig9
j = 5; % glycerol is 5
    % HiC
    T = Tm_HiC;
    [~, idx] = min( abs( T-times ) ); 
    URm_GLR_HiC=mean(UR5m_HiC(idx(1):idx(2),j))
    T = Tf_HiC;
    [~, idx] = min( abs( T-times ) ); 
    URf_GLR_HiC=mean(UR5f_HiC(idx(1):idx(2),j))
    % HiF
    T = Tm_HiF;
    [~, idx] = min( abs( T-times ) ); 
    URm_GLR_HiF=mean(UR5m_HiC(idx(1):idx(2),j))
    T = Tf_HiF;
    [~, idx] = min( abs( T-times ) ); 
    URf_GLR_HiF=mean(UR5f_HiC(idx(1):idx(2),j))
 
% % FFA
% j = 6; % FFA is 6
%     % HiC
%     T = Tm_HiC;
%     [~, idx] = min( abs( T-times ) ); 
%     URm_FFA_HiC=mean(UR5m_HiC(idx(1):idx(2),j));
%     T = Tf_HiC;
%     [~, idx] = min( abs( T-times ) ); 
%     URf_FFA_HiC=mean(UR5f_HiC(idx(1):idx(2),j));
%     % HiF
%     T = Tm_HiF;
%     [~, idx] = min( abs( T-times ) ); 
%     URm_FFA_HiF=mean(UR5m_HiC(idx(1):idx(2),j));
%     T = Tf_HiF;
%     [~, idx] = min( abs( T-times ) ); 
%     URf_FFA_HiF=mean(UR5f_HiC(idx(1):idx(2),j));

% FLUXES
% we are interested in the following fluxes 
% PYR --> GAP; flux 4
% GLY --> G6P; flux 8
% GLR --> GRP; flux 11
% GRP --> GAP; flux 13
% FFA --> ACOA; flux 17
% Sex differences exist for F4, F13, F8
j = [4 13 8];
    % HiC
    disp('High carbohydrate meal')
    T = Tm_HiC;
    [~, idx] = min( abs( T-times ) ); 
    F5m_GLR_HiC=mean(F5m_HiC(idx(1):idx(2),j), 1)
    T = Tf_HiC;
    [~, idx] = min( abs( T-times ) ); 
    F5f_GLR_HiC=mean(F5f_HiC(idx(1):idx(2),j), 1)

    % HiF
    disp('High fat meal')
    T = Tm_HiF;
    [~, idx] = min( abs( T-times ) ); 
    F5m_GLR_HiF=mean(F5m_HiF(idx(1):idx(2),j), 1)
    T = Tf_HiF;
    [~, idx] = min( abs( T-times ) ); 
    F5f_GLR_HiF=mean(F5f_HiF(idx(1):idx(2),j), 1)

% bars
% order is F4, F13, F8; M then F
% bars_gng_gyg = [0.9817 0.1093 0.5649; 0.9902 0.1209 0.5251; 0 0 0; % HIC
%                 0.8468 0.1086 0.5361; 0.8580 0.1189 0.4967; 0 0 0]; %HIF

% compute male to female ratio
bars_gng_gyg = [0.9902/0.9817 0.1209/0.1093 0.5251/0.5649; 0 0 0; % HIC
                0.8580/0.8468 0.1189/0.1086 0.4967/0.5361; 0 0 0]; %HIF
% convert to percent difference F/M
bars_gng_gyg([1,3],:) = (bars_gng_gyg([1,3],:)-1)*100;

fg = figure(6);
cols = [0.93,0.97,0.69; 0.50,0.80,0.73; 0.17,0.50,0.72];
colororder(cols)
b = bar(round(bars_gng_gyg,2,'significant'),'stacked'); 
b(1).BarWidth = 0.9;
legend('gluconeogenesis (pyruvate)', 'gluconeogenesis (glycerol)', 'glycogenolysis', ...
    'Location','Northeast','Orientation','vertical','NumColumns',1)
ylabel({'Endogenous glucose production'; '%\Delta F/M'})
set(gca,'XTickLabel', {'HiC' '' 'HiF' ''});
set(gca,'FontSize',18,'FontName','Times New Roman')
