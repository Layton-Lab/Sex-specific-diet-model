%% [3] - MUSCLE COMPARTMENT
function [F3,P3,U3] = muscle(t,tstart,tend,trecovery,y3,y,sex,WR,MR,RMR,epi,epi0,c0,d0)

    %Initialize
    ptemp = zeros(22,1);
    utemp = zeros(22,1);

    %Assign parameters to estimate
    lambda1  = c0(19); alpha1  = c0(20);
    lambda8  = c0(21); alpha8  = c0(22);
    lambda14 = c0(23); alpha14 = c0(24);
    lambda17 = c0(25); alpha17 = c0(26);
    lambda19 = c0(27); alpha19 = c0(28);
    if sex == 0  % male
        lactate_muscle = 1; 
    else    % sex==1, female
        lactate_muscle = c0(48); %=0.85 <1;
    end
   

    %Vmax values
    Vmax(1) = 0.398;
    Vmax(2) = 0.66;
    Vmax(3) = 5.28;
    Vmax(4) = 0.0;
    Vmax(5) = 0.0;
    Vmax(6) = 0.0;
    Vmax(7) = 0.5;
    Vmax(8) = 1.0;
    Vmax(9) = 14.85*lactate_muscle;
    Vmax(10) = 12.51;
    Vmax(11) = 0.508;
    Vmax(12) = 0.0;
    Vmax(13) = 0.0;
    Vmax(14) = 0.08;
    Vmax(15) = 0.0;
    Vmax(16) = 2.745;
    Vmax(17) = 0.701;
    Vmax(18) = 0.0;
    Vmax(19) = 0.26;
    Vmax(20) = 3.048;
    Vmax(21) = 9.968;
    Vmax(22) = 14.68;
    Vmax(23) = 80.0;
    Vmax(24) = 80.0;
    Vmax(25) = MR(3)*2.0;

    % RMR = Relative Metabolic Rate --> ATP Hydrolysis Related to Work Rate
    if (t>tstart) && (t<tend+trecovery) && (WR>0)
        Vmax(1)=Vmax(1)*(1+lambda1*(epi-epi0)^2/(alpha1+(epi-epi0)^2));
        Vmax(2)=Vmax(2)*RMR;
        Vmax(3)=Vmax(3)*RMR;
        Vmax(8)=Vmax(8)*RMR*(1+lambda8*(epi-epi0)^2/(alpha8+(epi-epi0)^2));
        Vmax(14)=Vmax(14)*(1+lambda14*(epi-epi0)^2/(alpha14+(epi-epi0)^2));
        Vmax(16)=Vmax(16)*RMR;
        Vmax(17)=Vmax(17)*(1+lambda17*(epi-epi0)^2/(alpha17+(epi-epi0)^2));
        Vmax(19)=Vmax(19)*(1+lambda19*(epi-epi0)^2/(alpha19+(epi-epi0)^2));
        Vmax(21)=Vmax(21)*RMR; 
        Vmax(22)=Vmax(22)*RMR; 
        % Notes: 
        % 1) RMR slowly returns to basal levels after exercise ends --> consequently, Vmax([2,3,8,16,21,22]) do too
        % 2) epinephrine also slowly returns to basal levels after exercise ends --> Vmax([1,8,14,17,19]) do too
    else %(WR==0)        
        % glucose metabolism
        Vmax_ins_M = d0(7);  
        Km_ins_M   = d0(8);   
        n_ins_M    = d0(9);
        beta_M     = hill(y(156),Vmax_ins_M,Km_ins_M,n_ins_M);

        Vmax_ins_M2= d0(10); 
        Km_ins_M2  = d0(11);  
        n_ins_M2   = d0(12);
        beta_M2    = hill(y(156),Vmax_ins_M2,Km_ins_M2,n_ins_M2);

        Vmax_ins_inh_M = d0(13); 
        Km_ins_inh_M   = d0(14);  
        n_ins_inh_M    = d0(15);
        gamma_M        = hill(y(156),Vmax_ins_inh_M,Km_ins_inh_M,n_ins_inh_M);
        
        % fat metabolism
        rho = 1-1/(1+exp(-4*(y(1)-4.4))); % fat metabolism activation function

        Vmax_ins_M_fat = d0(16);  
        Km_ins_M_fat   = d0(17);   
        n_ins_M_fat    = d0(18);
        beta_M_fat     = rho*hill(y(156),Vmax_ins_M_fat,Km_ins_M_fat,n_ins_M_fat);

        Vmax_ins_M_fat2= d0(19); 
        Km_ins_M_fat2  = d0(20);  
        n_ins_M_fat2   = d0(21);
        beta_M_fat2    = rho*hill(y(156),Vmax_ins_M_fat2,Km_ins_M_fat2,n_ins_M_fat2);

        Vmax_ins_inh_M_fat = d0(22); 
        Km_ins_inh_M_fat   = d0(23);  
        n_ins_inh_M_fat    = d0(24);
        gamma_M_fat        = rho*hill(y(156),Vmax_ins_inh_M_fat,Km_ins_inh_M_fat,n_ins_inh_M_fat);

        Vmax(1)  = Vmax(1)*(1+beta_M);
        Vmax(2)  = Vmax(2)*(1+beta_M2); 
        Vmax(3)  = Vmax(3)*(1+beta_M2); 
        Vmax(7)  = Vmax(7)*(1+beta_M/2.5-gamma_M_fat); % Note: beta_M/2.5
        
        Vmax(16) = Vmax(16)*(1+beta_M_fat2);

        Vmax(17) = Vmax(17)*(1+beta_M_fat);
        Vmax(8)  = Vmax(8)*(1-gamma_M+beta_M_fat); 
        Vmax(14) = Vmax(14)*(1+beta_M_fat);
    end 

    %Km values
    Km(1) = 0.1;    %GLC
    Km(2) = 0.048;  %PYR
    Km(3) = 1.44;   %LAC
    Km(4) = 1.3;    %ALA
    Km(5) = 0.064;  %GLR
    Km(6) = 0.53;   %FFA
    if sex == 0 % male
        Km(7) = 14.8;   %TGL
    else  % female
        Km(7) = 14.8*1.28;   %TGL % Km is chosen to match the IC for TG; see y3(7) in call_IFMod
    end
    Km(8) = 7.0e-4; %O2
    Km(9) = 15.43;  %CO2
    Km(10) = 0.24;  %G6P
    Km(11) = 95.0;  %GLY
    Km(12) = 0.08;  %GAP
    Km(13) = 0.15;  %GRP
    Km(14) = 0.0022;%ACoA
    Km(15) = 0.018; %CoA
    Km(16) = 0.45;  %NAD+
    Km(17) = 0.05;  %NADH
    Km(18) = 6.15;  %ATP
    Km(19) = 0.02;  %ADP
    Km(20) = 2.70;  %Pi
    Km(21) = 20.1;  %PCR
    Km(22) = 10.45; %CR
    Km(23) = Km(16)/Km(17); %NAD/NADH
    Km(24) = Km(17)/Km(16); %NADH/NAD
    Km(25) = Km(18)/Km(19); %ATP/ADP
    Km(26) = Km(19)/Km(18); %ADP/ATP

    [F3,P3,U3] = callflux(y3,Vmax,Km,ptemp,utemp,3);

end
% Activation fucntion
function val = hill(x, Vmax, Km, n)
val=Vmax*x^n/(Km^n+x^n);
end