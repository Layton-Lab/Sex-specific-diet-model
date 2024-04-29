%% CALIBRATION AND VALIDATION
clear;
clc;
addpath( './helper_functions' );

%% RUN SIMULATIONS FOR ALL 4 MEAL COMPOSITIONS
t0=0; 
WR=0; % NO EXERCISE
tstart=0; % unused, only for exercise
tend=0; % unused, only for exercise
ndays=1; % number of days of IF
lag = 24; % single meal
% SET 1 - Frayn1993
carbs=96; fats=33; % amounts in grams per day 
sex = 0; % male
    [Tm1,Ym1,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm]=call_IFMod(t0,sex,WR,carbs,fats,ndays,lag,tstart,tend);
    Tm1=(Tm1-tstart)/60;
sex = 1; % female
    [Tf1,Yf1,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm]=call_IFMod(t0,sex,WR,carbs,fats,ndays,lag,tstart,tend);
    Tf1=(Tf1-tstart)/60;

% SET 2 - Taylor1996
carbs=139; fats=17; % amounts in grams per day 
sex = 0; % male
    [Tm2,Ym2,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm]=call_IFMod(t0,sex,WR,carbs,fats,ndays,lag,tstart,tend);
    Tm2=(Tm2-tstart)/60;
sex = 1; % female
    [Tf2,Yf2,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm]=call_IFMod(t0,sex,WR,carbs,fats,ndays,lag,tstart,tend);
    Tf2=(Tf2-tstart)/60;

% SET 3 - Hansen2011
carbs=58; fats=27.7; % amounts in grams per day 
sex = 0; % male
    [Tm3,Ym3,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm]=call_IFMod(t0,sex,WR,carbs,fats,ndays,lag,tstart,tend);
    Tm3=(Tm3-tstart)/60;
sex = 1; % female
    [Tf3,Yf3,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm]=call_IFMod(t0,sex,WR,carbs,fats,ndays,lag,tstart,tend);
    Tf3=(Tf3-tstart)/60;

% SET 4 - Taylor1993
carbs=289; fats=45; % amounts in grams per day 
sex = 0; % male
    [Tm4,Ym4,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm]=call_IFMod(t0,sex,WR,carbs,fats,ndays,lag,tstart,tend);
    Tm4=(Tm4-tstart)/60;
sex = 1; % female
    [Tf4,Yf4,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm]=call_IFMod(t0,sex,WR,carbs,fats,ndays,lag,tstart,tend);
    Tf4=(Tf4-tstart)/60;

%% FIGURES

%--------------------------------------------------------------------------
% FIGURE 2
% Time profile of plasma insulin and glucose concentrations after 
% an overnight fast and following a single meal.
%--------------------------------------------------------------------------
N = 4;
C = brewermap(N,'BrBG'); 
cols = C(4,:); 
colororder(cols)

fg = figure(1);
% Row 1
% Frayn1993
% Info: eight healthy subjects (including three females) were studied after an overnight fast. 
% The meal contained 96 g carbohydrate and 33 g fat. 
t = [0 30 60 90 120 180 240 300 360]';
t=t/60;
glucose = [5 8 8.2 7.5 6.2 5.2 5.1 5.2 4.9]';
glucose_err = [0.42 0.60 0.66 0.26 0.39 0.61 0.84 0.48 0.62]';
insulin = [55 500 490 470 330 180 180 110 50]';
insulin_err=[4.0 66.7 60.0 36.1 60.9 47.6 51.2 33.0 15.4]';
subplot(4, 2, 1)
errorbar(t,insulin,insulin_err,'ks','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k')
hold on
colororder(cols)
plot(Tm1,Ym1(:,156),'LineWidth',3)
hold off
xlabel('Time, hr')
ylabel('Insulin, pM')
xlim([0,12])
ylim([0 650])
xticks(0:4:12)
grid on
box on
legend('Experiment 1','Male model','Fontsize',16)
set(gca,'FontSize',16,'FontName','Times New Roman')
subplot(4, 2, 2)
errorbar(t,glucose,glucose_err,'ks','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k')
hold on
colororder(cols)
plot(Tm1,Ym1(:,1),'LineWidth',3)
hold off
xlabel('Time, hr')
ylabel('Glucose, mM')
xlim([0,12])
ylim([4 9.5])
xticks(0:4:12)
grid on
box on
legend('Experiment 1','Male model','Fontsize',16)
set(gca,'FontSize',16,'FontName','Times New Roman')

% Row 2
% Taylor1996
% Info: Healthy volunteers, aged 18–40 yr were recruited (males and females)
% liquid test meal (824 kcal; 67.3% carbohydrate as glucose, 18.5% fat, 14.2% protein) 
% was consumed over a 10-min period.This is roughly 139g of carbs, 17g of fat, and 29g of protein.
t3_2=[-10	0	10	20	30	40	50	60	70	80	90	100	110	120	180	240	300	360	420	480	540	600];
glucose=[5.03	5.00	5.76	6.85	7.90	8.44	8.60	8.57	8.25	7.61	7.54	7.45	7.19	7.03	6.84	6.42	5.23	4.91	4.94	4.77	4.89	4.79];
glucose_err=[0.19	0.13	0.45	0.54	0.19	0.35	0.57	0.70	0.61	0.61	0.48	0.51	0.45	0.57	0.38	0.41	0.61	0.35	0.13	0.35	0.19	0.26];
t3_2=t3_2/60; 
t3_22=[-10	0	10	20	30	40	60	90	120	180	240	300	360	420	480	540	600];	
insulin=[25.58	25.61	143.03	264.37	438.53	403.34	436.66	366.31	356.62	276.58	194.58	67.58	30.59	30.76	27.02	27.21	29.35];	
insulin_err=[3.92	3.92	54.78	86.09	78.26	54.79	68.48	48.91	66.52	60.65	41.09	23.48	7.83	5.87	7.83	7.82	5.87];	
t3_22=t3_22/60; 
subplot(4, 2, 3)
errorbar(t3_22,insulin,insulin_err,'k^','MarkerSize',9,'MarkerEdgeColor','k','MarkerFaceColor','k')
hold on
colororder(cols)
plot(Tm2,Ym2(:,156),'LineWidth',3)
hold off
xlabel('Time, hr')
ylabel('Insulin, pM')
xlim([0,12])
ylim([0 650])
xticks(0:4:12)
grid on
box on
legend('Experiment 2','Male model','Fontsize',16)
set(gca,'FontSize',16,'FontName','Times New Roman')
subplot(4, 2, 4)
errorbar(t3_2,glucose,glucose_err,'k^','MarkerSize',9,'MarkerEdgeColor','k','MarkerFaceColor','k')
hold on
colororder(cols)
plot(Tm2,Ym2(:,1),'LineWidth',3)
hold off
xlabel('Time, hr')
ylabel('Glucose, mM')
xlim([0,12])
ylim([4 9.5])
xticks(0:4:12)
grid on
box on
legend('Experiment 2','Male model','Fontsize',16)
set(gca,'FontSize',16,'FontName','Times New Roman')

% Row 3
% Hansen2011
% Info: Ten healthy young Caucasian men were studied after a 10h overnight fast
% The subjects ingested 100 g of NAN 1 [2200 KJ (520 kcal): 
% 58 g carbohydrate, 27.7 g fat, and 9.5 g protein] dissolved in 300 ml water over 5 min
t2=[-15 -10 0 10 20 30 40 50 60 75 90 120 150 180 240]; % in mins
glucose = [4.98	4.98	4.91	5.19	5.51	5.97	6.33	6.19	5.83	5.40	5.25	4.89	4.85	4.71	4.74]';
glucose_err = [0.18	0.18	0.12	0.14	0.11	0.18	0.14	0.32	0.36	0.28	0.36	0.28	0.21	0.14	0.14]';
insulin = [36.47	37.32	36.47	67.48	124.01	226.14	297.26	306.38	271.73	227.96	182.37	131.31	72.95	45.59	34.65]';
insulin_err=[7.30	5.11	3.77	29.18	41.94	74.78	71.12	94.83	102.13	87.53	52.89	38.30	20.06	12.77	9.12]';
t2=t2/60;
subplot(4, 2, 5)
errorbar(t2,insulin,insulin_err,'k^','MarkerSize',9,'MarkerEdgeColor','k','MarkerFaceColor','k')
hold on
colororder(cols)
plot(Tm3,Ym3(:,156),'LineWidth',3)
hold off
xlabel('Time, hr')
ylabel('Insulin, pM')
xlim([0,12])
ylim([0 650])
xticks(0:4:12)
grid on
box on
legend('Experiment 3','Male model','Fontsize',16)
set(gca,'FontSize',16,'FontName','Times New Roman')
subplot(4, 2, 6)
errorbar(t2,glucose,glucose_err,'k^','MarkerSize',9,'MarkerEdgeColor','k','MarkerFaceColor','k')
hold on
colororder(cols)
plot(Tm3,Ym3(:,1),'LineWidth',3)
hold off
xlabel('Time, hr')
ylabel('Glucose, mM')
xlim([0,12])
ylim([4 9.5])
xticks(0:4:12)
grid on
box on
legend('Experiment 3','Male model','Fontsize',16)
set(gca,'FontSize',16,'FontName','Times New Roman')

% Row 4
% Taylor 1993
% not using data for insulin and glucose
subplot(4, 2, 7)
plot(Tm4,Ym4(:,156),'Color',C(4,:),'LineWidth',3)
xlabel('Time, hr')
ylabel('Insulin, pM')
xlim([0,12])
ylim([0 650])
xticks(0:4:12)
grid on
box on
legend(['Male model' newline '(experiment 4)'],'Fontsize',16)
set(gca,'FontSize',16,'FontName','Times New Roman')
subplot(4, 2, 8)
plot(Tm4,Ym4(:,1),'Color',C(4,:),'LineWidth',3)
xlabel('Time, hr')
ylabel('Glucose, mM')
xlim([0,12])
ylim([3 11])
xticks(0:4:12)
grid on
box on
legend(['Male model' newline '(experiment 4)'],'Fontsize',16)
set(gca,'FontSize',16,'FontName','Times New Roman')
% add letters to captions
AddLetters2Plots(fg, {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)','(j)','(h)'},...
'HShift', -0.1, 'VShift', -0.06, 'Direction', 'LeftRight', ...
'FontSize', 18, 'FontWeight', 'normal')
% export to png
% % exportgraphics(gcf,'Fig2.png','resolution',600)

%--------------------------------------------------------------------------
% FIGURE 3
% Time profile of glycogen concentration in liver (left column) and 
% skeletal muscle (right column), relative to its initial value, after an 
% overnight fast and following a single meal.
%--------------------------------------------------------------------------
N = 4;
C = brewermap(N,'BrBG'); 
cols =  [C(2,:); C(4,:)]; 
colororder(cols)
fg = figure(2);
% Row 1
% Frayn1993
% no glycogen data
subplot(4, 2, 1)
plot(nan,nan,'w','LineWidth',2);
hold on
colororder(cols)
plot(Tm1,Ym1(:,121)./Ym1(1,121),'LineWidth',3)
plot(Tf1,Yf1(:,121)./Yf1(1,121),'LineWidth',3)
hold off
xlabel('Time, hr')
ylabel('Relative glycogen')
xlim([0,12])
ylim([0.9 2.2])
xticks(0:4:12)
grid on
box on
legend('Experiemnt 1', 'Male model', 'Female model','Fontsize',16)
set(gca,'FontSize',16,'FontName','Times New Roman')
subplot(4, 2, 2)
hold on
plot(nan,nan,'w','LineWidth',2);
colororder(cols)
plot(Tm1,Ym1(:,77)./Ym1(1,77),'LineWidth',3)
plot(Tf1,Yf1(:,77)./Yf1(1,77),'LineWidth',3)
hold off
xlabel('Time, hr')
ylabel('Relative glycogen')
xlim([0,12])
ylim([0.9 1.3])
xticks(0:4:12)
yticks(0.9:0.2:1.3)
grid on
box on
legend('Experiment 1','Fontsize',16)
set(gca,'FontSize',16,'FontName','Times New Roman')

% Row 2
% Taylor 1996 - Liver
% Info: Healthy volunteers, aged 18–40 yr were recruited (males and females)
% liquid test meal (824 kcal; 67.3% carbohydrate as glucose, 18.5% fat, 14.2% protein) 
% was consumed over a 10-min period.This is roughly 139g of carbs, 17g of fat, and 29g of protein.
t3 = [0 120 150 180 260 290 360 470 560]';
glycogen_L = [206.88 248.53 260.60 280.85 296.01 289.37 294.39 258.85 242.81]'/206.88;
glycogen_L_err = [22.22 18.71 17.93 20.27 23.00 19.49 19.10 21.83 21.83]'/206.88;
t3=t3/60; 
subplot(4, 2, 3)
errorbar(t3,glycogen_L,glycogen_L_err,'ks','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k')
hold on
colororder(cols)
plot(Tm2,Ym2(:,121)./Ym2(1,121),'LineWidth',3)
plot(Tf2,Yf2(:,121)./Yf2(1,121),'LineWidth',3)
hold off
xlabel('Time, hr')
ylabel('Relative glycogen')
xlim([0,12])
ylim([0.9 2.2])
xticks(0:4:12)
grid on
box on
legend('Experiment 2','Fontsize',16)
set(gca,'FontSize',16,'FontName','Times New Roman')
subplot(4, 2, 4)
hold on
plot(nan,nan,'w','LineWidth',2);
colororder(cols)
plot(Tm2,Ym2(:,77)./Ym2(1,77),'LineWidth',3)
plot(Tf2,Yf2(:,77)./Yf2(1,77),'LineWidth',3)
hold off
xlabel('Time, hr')
ylabel('Relative glycogen')
xlim([0,12])
ylim([0.9 1.3])
xticks(0:4:12)
yticks(0.9:0.2:1.3)
grid on
box on
legend('Experiment 2','Fontsize',16)
set(gca,'FontSize',16,'FontName','Times New Roman')

% Row 3
% Hansen 2011
% No glycogen data 
subplot(4, 2, 5)
plot(nan,nan,'w','LineWidth',2);
hold on
colororder(cols)
plot(Tm3,Ym3(:,121)./Ym3(1,121),'LineWidth',3)
plot(Tf3,Yf3(:,121)./Yf3(1,121),'LineWidth',3)
hold off
xlabel('Time, hr')
ylabel('Relative glycogen')
xlim([0,12])
ylim([0.9 2.2])
xticks(0:4:12)
grid on
box on
legend('Experiment 3','Fontsize',16)
set(gca,'FontSize',16,'FontName','Times New Roman')
subplot(4, 2, 6)
plot(nan,nan,'w','LineWidth',2);
hold on
colororder(cols)
plot(Tm3,Ym3(:,77)./Ym3(1,77),'LineWidth',3)
plot(Tf3,Yf3(:,77)./Yf3(1,77),'LineWidth',3)
hold off
xlabel('Time, hr')
ylabel('Relative glycogen')
xlim([0,12])
ylim([0.9 1.3])
xticks(0:4:12)
yticks(0.9:0.2:1.3)
grid on
box on
legend('Experiment 3','Fontsize',16)
set(gca,'FontSize',16,'FontName','Times New Roman')

% Row 4
% Taylor 1993 - Muscle
% Info: Eight healthy young volunteers (males and females)
% the test meal consisted of 289 g of carbs, 45g of fat, and 89g of protein.									
t4 = [-60 -40 -20 60 120 180 240 300 360 420]';
glycogen_M=[83.3	86.84	87.835	84.065	89.779	92.621	100.2	98.161	99.568	94.514]'/83.3;
glycogen_M_err=[5.224	7.064	6.344	6.583	6.344	5.944	6.7	6.6	6.52	5.063]'/83.3;
t4=t4/60; 
subplot(4, 2, 7)
plot(nan,nan,'w','LineWidth',2);
hold on
colororder(cols)
plot(Tm4,Ym4(:,121)./Ym4(1,121),'LineWidth',3)
plot(Tf4,Yf4(:,121)./Yf4(1,121),'LineWidth',3)
hold off
xlabel('Time, hr')
ylabel('Relative glycogen')
xlim([0,12])
ylim([0.9 2.2])
xticks(0:4:12)
grid on
box on
legend('Experiment 4','Fontsize',16,'location','southeast')
set(gca,'FontSize',16,'FontName','Times New Roman')
subplot(4, 2, 8)
errorbar(t4,glycogen_M,glycogen_M_err,'ks','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k')
hold on
colororder(cols)
plot(Tm4,Ym4(:,77)./Ym4(1,77),'LineWidth',3)
plot(Tf4,Yf4(:,77)./Yf4(1,77),'LineWidth',3)
hold off
xlabel('Time, hr')
ylabel('Relative glycogen')
xlim([0,12])
ylim([0.9 1.3])
xticks(0:4:12)
yticks(0.9:0.2:1.3)
grid on
box on
legend('Experiment 4','Fontsize',16)
set(gca,'FontSize',16,'FontName','Times New Roman')
% add letters to captions
AddLetters2Plots(fg, {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)',},...
'HShift', -0.1, 'VShift', -0.06, 'Direction', 'LeftRight', ...
'FontSize', 18, 'FontWeight', 'normal')
% export to png
% % exportgraphics(gcf,'Fig3.png','resolution',600)

%--------------------------------------------------------------------------
% FIGURE 4
% Time profile of plasma metabolite concentrations after an overnight 
% fast and following a single meal.
%--------------------------------------------------------------------------
N = 4;
C = brewermap(N,'BrBG'); 
cols = [C(2,:); C(4,:)];
colororder(cols)
fg = figure(3);
% Frayn1993 
% Info: eight healthy subjects (including three females) were studied after an overnight fast. 
% The meal contained 96 g carbohydrate and 33 g fat. 
t = [0 30 60 90 120 180 240 300 360]';
ffa = [0.52 0.39 0.07 0.06 0.06 0.1 0.2 0.32 0.5]'/0.52;
ffa_err = [0.11 0.08 0.05 0.06 0.06 0.08 0.06 0.06 0.09]'/0.52;
% Coppack1990
% Info: seven healthy subjects (including four females) were studied after an overnight fast. 
% The meal contained 96 g carbohydrate and 33 g fat. 
% t = [0 30 60 90 120 180 240 300 360]';
lactate = [0.52 0.93 1.64 1.28 0.97 0.71 0.65 0.54 0.51]'/0.52;
lactate_err = [0.10 0.11 0.21 0.10 0.05 0.07 0.07 0.04 0.05]'/0.52;
TG = [0.65 0.68 0.74 0.72 0.80 1.03 1.20 1.15 0.88]'/0.65;
TG_err = [0.11 0.10 0.10 0.10 0.11 0.15 0.15 0.15 0.08]'/0.65;
glycerol = [0.0491 0.0307 0.0232 0.0237 0.0207 0.0235 0.0360 0.0471 0.0583]'/0.0491;
glycerol_err = [0.0066 0.0066 0.0053 0.0044 0.0049 0.0049 0.0035 0.0084 0.0062]'/0.0491;
% convert experimental time from minutes to hours
t=t/60;
subplot(2, 2, 1)
errorbar(t,lactate,lactate_err,'ks','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k')
hold on
colororder(cols)
plot(Tm1,Ym1(:,3)./Ym1(1,3),'LineWidth',3)
plot(Tf1,Yf1(:,3)./Yf1(1,3),'LineWidth',3)
hold off
xlabel('Time, hr')
ylabel('Relative lactate')
xlim([0,12])
xticks(0:4:12)
grid on
box on
legend('Experiemnt 1', 'Male model', 'Female model','Fontsize',16)
set(gca,'FontSize',16,'FontName','Times New Roman')

subplot(2, 2, 2)
hold on
errorbar(t,TG,TG_err,'ks','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k')
colororder(cols)
plot(Tm1,Ym1(:,7)./Ym1(1,7),'LineWidth',3)
plot(Tf1,Yf1(:,7)./Yf1(1,7),'LineWidth',3)
hold off
xlabel('Time, hr')
ylabel('Relative TG')
xlim([0,12])
xticks(0:4:12)
grid on
box on
legend('Experiemnt 1', 'Male model', 'Female model','Fontsize',16)
set(gca,'FontSize',16,'FontName','Times New Roman')

subplot(2, 2, 3)
errorbar(t,ffa,ffa_err,'ks','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k')
hold on
colororder(cols)
plot(Tm2,Ym2(:,6)./Ym2(1,6),'LineWidth',3)
plot(Tf2,Yf2(:,6)./Yf2(1,6),'LineWidth',3)
hold off
xlabel('Time, hr')
ylabel('Relative FFA')
xlim([0,12])
xticks(0:4:12)
grid on
box on
legend('Experiemnt 1', 'Male model', 'Female model','Fontsize',16)
set(gca,'FontSize',16,'FontName','Times New Roman')
subplot(2, 2, 4)
hold on
errorbar(t,glycerol,glycerol_err,'ks','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k')
colororder(cols)
plot(Tm2,Ym2(:,5)./Ym2(1,5),'LineWidth',3)
plot(Tf2,Yf2(:,5)./Yf2(1,5),'LineWidth',3)
hold off
xlabel('Time, hr')
ylabel('Relative glycerol')
xlim([0,12])
xticks(0:4:12)
grid on
box on
legend('Experiemnt 1', 'Male model', 'Female model','Fontsize',16)
set(gca,'FontSize',16,'FontName','Times New Roman')
% add letters to captions
AddLetters2Plots(fg, {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)',},...
'HShift', -0.08, 'VShift', -0.04, 'Direction', 'LeftRight', ...
'FontSize', 18, 'FontWeight', 'normal')
% export to png
% % exportgraphics(gcf,'Fig4.png','resolution',600)