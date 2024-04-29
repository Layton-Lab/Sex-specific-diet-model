%% FUNCTION CALL -- MODEL ODEs
function [dydt,bflow,MR,P,U,uprel,F,epi,v_feed_Glc_B,v_feed_TG_B] = IFMod(t,y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend,sm,sa)

    %Initialize variables
    dydt = zeros(156,1);
    UR1  = zeros(1,22);
    UR2  = zeros(1,22);
    UR3  = zeros(1,22);
    UR4  = zeros(1,22);
    UR5  = zeros(1,22);
    UR6  = zeros(1,22);
    UR7  = zeros(1,22);
    
    %Cluster substrates per compartment
    Ca(1:22) = y(1:22); 
    y1(1:22) = y(23:44);
    y2(1:22) = y(45:66);
    y3(1:22) = y(67:88);
    y4(1:22) = y(89:110);
    y5(1:22) = y(111:132);
    y6(1:22) = y(133:154);
    GIR = (y(155)/y(156))-glu0/ins0;

    % organ weights [whole body = 70kg]
    if sex == 0  % male
        Wbody = 70;
        WB = 1.4;    % Brain
        WH = 0.331;  % Heart
        WM = Wbody*(0.4+sm)*0.8; % 28*0.8; %Skeletal muscles, excluding upper extremities which account for 18-20%
        WG = 2.0;    % GI tract 
        WL = 1.8;    % Liver 
        WA = Wbody*(0.16+sa); %11;     % Adipose tissue, 16% body fat content
    else  % female
        Wbody = 58;
        WB = 1.2;    % Brain
        WH = 0.253;  % Heart
        WM = Wbody*(0.3+sm)*0.8; %17*0.8; % Skeletal muscles, excluding upper extremities which account for 18-20%
        WG = 2.0;    % GI tract 
        WL = 1.4;    % Liver 
        WA = Wbody*(0.295+sa); %17.1;   % Adipose tissue, 29.5% body fat content
    end
    WO = Wbody - WB - WH -WM - WG - WL - WA; % Others

    %Incorporating an integral rein control
    %phi and psi are formulated to maintain an arterial glucose
    %concentration at 5mM 
    if(Ca(1)<2.5) 
        phi=1.0;
        psi=0.0;
    elseif((Ca(1)>=2.5) && (Ca(1)<=7.5)) 
        phi=1.0-(Ca(1)-2.5)^2.0/25.0;
        psi=1.0-(Ca(1)-7.5)^2.0/25.0;
    elseif((Ca(1)>7.5))
        phi=0.0;
        psi=1.0;
    end

    %--------------------- START OF EXERCISE ---------------------%
    % post-exercise recovery
    if t<tend
        recovery=1;
    else
        % time constant for recovery is \tau_recovery = 5;
        % after 10 mins, basal levels are recovered up to 90%; 
        % after 30 mins, basal levels are fully recovered.
        recovery=exp(-(t-tend)/5);
    end
    % Assign a mandatory recovery time of 30 mins
    % This will allow us to stagger exercise and feeding
    trecovery = 30;

    %blood flow and metabolic rate related parameters
    sig_heart  = c0(1);
    sig_muscle = c0(2);
    sig_gi     = c0(3);
    tau_Q      = c0(4); 
    tau_Qm     = c0(47); 
    gamma_m    = c0(5);
    gamma_h    = c0(48);

    %Glucagon-insulin related parameters
    h  = c0(6);
    k1 = c0(7); 
    k2 = c0(8);
    k3 = c0(9);
    k4 = c0(10);
    k5 = c0(11); 
    D  = c0(12);

    %Epinephrine
    if sex==0
        epi0 = 250.0;
    else
        epi0 = 150.0;
    end
    epi  = epi0;
    %Metabolic rate of ATP --> ADP
        MR(1) = 15.20;
        MR(2) = 7.33;
        MR(3) = 10.82;
        MR(4) = 3.04;
        MR(5) = 13.92;
        MR(6) = 2.74;
        %ATP hydrolysis related to Work Rate
        RMR_M=MR(3)/10.82;
        RMR_H=MR(2)/7.33;
   
    %Blood flow in tissues
    Q(1)=0.75;
    Q(2)=0.25;
    Q(3)=0.9;
    Q(4)=1.1;
    Q(6)=0.36;
    if sex==0
        Q(5)=0.4;  % chosen so that Q4 + Q5 = 1.5;
        Q(7)=1.74; % chosen to keep blood volume at 5.5L/min
    else
        Q(5)=0.25; % chosen so that Q4 + Q5 = 1.35;
        Q(7)=1.39; % chosen to keep blood volume at 5L/min
    end

    if sex==0
        gain=1100;
    else
        gain=330;
    end
    if (t>tstart) && (WR>0)
        epi=epi0+gain*(1.0-exp(-(t-tstart)/30))*recovery;
        Q(2)=0.25+sig_heart*(1-exp(-(t-tstart)/tau_Q))*recovery;
        Q(3)=0.90+sig_muscle*(1-exp(-(t-tstart)/tau_Qm))*recovery;
        Q(4)=1.10+sig_gi*(1-exp(-(t-tstart)/tau_Q))*recovery;
        MR(2)=7.33+WR*gamma_h*recovery;
        MR(3)=10.82+WR*gamma_m*recovery;
       %ATP hydrolysis related to Work Rate
        RMR_M=MR(3)/10.82;
        RMR_H=MR(2)/7.33;
    end 

%--------------------- END OF EXERCISE ----------------------%

%--------------------- START OF FEEDING ---------------------%

    % pancreas-related parameters (specific to feeding)
    k_glu_b = 0.5; 
    k_ins_b = 1.8;
    Vmax_glc_glu = d0(67);
    Vmax_glc_ins = d0(68); 
    Ki_glc_glu   = d0(69);
    Ki_ins_glu   = d0(70);
    Km_glc_ins   = d0(71);
    n_glc_glu    = d0(72);
    n_ins_glu    = d0(73); 
    n_glc_ins    = d0(74);
%   D            = 0.1; % keep this at 0.1, same as for exercise

% Adjust UPTAKE sensitivities based on blood levels of glucose and TG
% Note: val = sens(x,ic,a,slope)
    % detect hyperglycemia
    hyperGLC=sens(Ca(1),9,1,2); 
    % inflection point (9) is chosen so that values Ca(1)>6.1 elicit a response

    % detect hypertriglyceridemia
    hyperTG=sens(Ca(7),4,1,2);
    % inflection point (4) is chosen so that values Ca(7)>2 elicit a response

%---------------------  END OF FEEDING  ---------------------%

% Uptake-release rates for the Brain (B) compartment - y1(x)
    % Partition coefficients
    if sex == 0  % male
        s(1) = 4.0119;  %GLC
        s(2) = 0.4533;  %PYR
        s(3) = 0.4828;  %LAC
        s(4) = 0.0;     %ALA
        s(5) = 0.0;     %GLR
        s(6) = 0.0;     %FFA
        s(7) = 0.0;     %TGL
        s(8) = 183.704; %O2
        s(9) = 1.6034;  %CO2
    else  % female
        s(1) = 3.9315;  %GLC vs ----> male value = 4.0119
        s(2) = 0.4533;  %PYR
        s(3) = 0.4828;  %LAC
        s(4) = 0.0;     %ALA
        s(5) = 0.0;     %GLR
        s(6) = 0.0;     %FFA
        s(7) = 0.0;     %TGL
        s(8) = 183.704; %O2
        s(9) = 1.6034;  %CO2
    end
    s(10:22) = 0.0;
    
    % Uptake/Release rates
    UR1(1) = Q(1) * (Ca(1) - s(1)*y1(1));
    UR1(2) = Q(1) * (Ca(2) - s(2)*y1(2));
    UR1(3) = Q(1) * (Ca(3) - s(3)*y1(3));
    UR1(4:7) = 0.0;
    UR1(8) = Q(1) * (Ca(8) - s(8)*y1(8));
    UR1(9) = Q(1) * (Ca(9) - s(9)*y1(9));
    UR1(10:22) = 0.0;

% Uptake-release rates for the Heart (H) compartment - y2(x)
    % Partition coefficients
    if sex == 0  % male
        s(1) = 4.84;    %GLC
        s(2) = 0.34;    %PYR
        s(3) = 0.1392;  %LAC
        s(4) = 0.0;     %ALA
        s(5) = 4.6667;  %GLR
        s(6) = 24.7619; %FFA
        s(7) = 0.3173;  %TGL
        s(8) = 3.3483;  %O2
        s(9) = 1.269;   %CO2
    else  % female
        s(1) = 4.75;    %GLC vs ----> male value = 4.84
        s(2) = 0.34;    %PYR
        s(3) = 0.1392;  %LAC
        s(4) = 0.0;     %ALA
        s(5) = 4.6667;  %GLR
        s(6) = 29.5238; %FFA vs ----> male value = 24.7619
        s(7) = 0.2981;  %TGL vs ----> male value = 0.3173
        s(8) = 3.3483;  %O2
        s(9) = 1.269;   %CO2
    end
    s(10:22) = 0.0;

    % Uptake/Release rates
    UR2(1) = Q(2) * (Ca(1) - s(1)*y2(1));
    UR2(2) = Q(2) * (Ca(2) - s(2)*y2(2));
    UR2(3) = Q(2) * (Ca(3) - s(3)*y2(3));
    UR2(4) = 0.0;
    UR2(5) = Q(2) * (Ca(5) - s(5)*y2(5));
    UR2(6) = Q(2) * (Ca(6) - s(6)*y2(6));
    UR2(7) = 0.0;
    UR2(8) = Q(2) * (Ca(8) - s(8)*y2(8));
    UR2(9) = Q(2) * (Ca(9) - s(9)*y2(9));
    UR2(10:22) = 0.0;

% Uptake-release rates for the Muscle (M) compartment - y3(x)
    % Partition coefficients
    if sex == 0  % male
        s(1) = 10.03472*hyperGLC; %GLC
        s(2) = 1.3009;  %PYR
        s(3) = 0.5725;  %LAC
        s(4) = 0.1819;  %ALA
        s(5) = 1.1458;  %GLR
        s(6) = 1.1488;  %FFA
        s(7) = 0.06667; %TGL
        s(8) = 12.1723; %O2
        s(9) = 1.5092;  %CO2
    else  % female
        s(1) = 9.8472*hyperGLC;  %GLC
        s(2) = 1.3530;  %PYR
        s(3) = 0.5596;  %LAC
        s(4) = 0.1819;  %ALA
        s(5) = 1.1458;  %GLR
        s(6) = 1.3375;  %FFA vs ----> male value = 1.1488
        s(7) = 0.0489;  %TGL vs ----> male value = 0.06667
        s(8) = 12.1723; %O2
        s(9) = 1.5092;  %C02;
    end
    s(10:22) = 0.0;
    
    %Uptake/Release rates
    UR3(1) = Q(3) * (Ca(1) - s(1)*y3(1));
    UR3(2) = Q(3) * (Ca(2) - s(2)*y3(2));
    UR3(3) = Q(3) * (Ca(3) - s(3)*y3(3));
    UR3(4) = Q(3) * (Ca(4) - s(4)*y3(4));
    UR3(5) = Q(3) * (Ca(5) - s(5)*y3(5));
    UR3(6) = Q(3) * (Ca(6) - s(6)*y3(6));
    UR3(7) = Q(3) * (Ca(7) - s(7)*y3(7));
    UR3(8) = Q(3) * (Ca(8) - s(8)*y3(8));
    UR3(9) = Q(3) * (Ca(9) - s(9)*y3(9));
    UR3(10:22) = 0.0;
 
% Uptake-release rates for the GI tract (G) compartment - y4(x)  
    % Blood/tissue partition coefficients
    if sex == 0  % male
        s(1) = 4.9309;   %GLC
        s(2) = 0.34;     %PYR
        s(3) = 0.1804;   %LAC
        s(4) = 0.0;      %ALA
        s(5) = 7.0909;   %GLR
        s(6) = 36.6234;  %FFA
        s(7) = 0.002187; %TGL
        s(8) = 15.4805;  %O2
        s(9) = 1.4332;   %CO2
    else  % female
        s(1) = 4.8409;   %GLC vs ----> male value = 4.9309
        s(2) = 0.34;     %PYR
        s(3) = 0.1804;   %LAC
        s(4) = 0.0;      %ALA
        s(5) = 7.0909;   %GLR
        s(6) = 41.3853;  %FFA vs ----> male value = 36.6234
        s(7) = 0.0020545;%TGL vs ----> male value = 0.002187; 
        s(8) = 15.4805;  %O2
        s(9) = 1.4332;   %CO2
    end
    s(10:22) = 0.0;
   
    % Uptake/Release rates
    UR4(1) = Q(4) * (Ca(1) - s(1)*y4(1));
    UR4(2) = Q(4) * (Ca(2) - s(2)*y4(2));
    UR4(3) = Q(4) * (Ca(3) - s(3)*y4(3));
    UR4(4) = 0.0;
    UR4(5) = Q(4) * (Ca(5) - s(5)*y4(5));
    UR4(6) = Q(4) * (Ca(6) - s(6)*y4(6));
    UR4(7) = Q(4) * (Ca(7) - s(7)*y4(7));
    UR4(8) = Q(4) * (Ca(8) - s(8)*y4(8));
    UR4(9) = Q(4) * (Ca(9) - s(9)*y4(9));
    UR4(10:22) = 0.0;

%Uptake-release rates for the Liver (L) compartment - y5(x)  
    % Blood/tissue partition coefficients
    if sex == 0  % male
        s(1) = 0.6796*hyperGLC; %GLC
        s(2) = 0.1838; %PYR
        s(3) = 0.6341; %LAC
        s(4) = -0.10;  %ALA
        s(5) = 0.0481; %GLR
        s(6) = 1.0526; %FA
        s(7) = 0.3430; %TGL
        s(8) = 215.90; %O2
        s(9) = 1.5104; %CO2
    else  % female
        s(1) = 0.6744*hyperGLC; %GLC vs ----> male value = 0.6796
        s(2) = 0.1838; %PYR
        s(3) = 0.6098; %LAC vs ----> male value = 0.6341
        s(4) = -0.1958;%ALA vs ----> male value = -0.10;
        s(5) = 0.0481; %GLR vs ----> male value = 0.0481 
        s(6) = 1.2164; %FFA vs ----> male value = 1.0526
        s(7) = 0.3232; %TGL vs ----> male value = 0.3430; 
        s(8) = 206.9738; %O2
        s(9) = 1.5219; %CO2
    end
    s(10:22) = 0.0;

    %Uptake/Release rates
    %Relative to the liver
    %The blood input to the liver comes from the hepatic artery AND venous blood of the GI tract.
    Cp(1:3)=(Q(4)*(Ca(1:3) - UR4(1:3)/Q(4))+Q(5)*Ca(1:3))/(Q(4)+Q(5));
    Cp(4)=Ca(4);
    Cp(5:7)=(Q(4)*(Ca(5:7) - UR4(5:7)/Q(4))+Q(5)*Ca(5:7))/(Q(4)+Q(5));
    Cp(8:9)=(Q(4)*(Ca(8:9) - UR4(8:9)/Q(4))+Q(5)*Ca(8:9))/(Q(4)+Q(5));
    Cp(10:22)=Ca(10:22);
    
    %Uptake/Release rates
    UR5(1) = (Q(4)+Q(5)).*(Cp(1) - s(1).*y5(1));
    UR5(2) = (Q(4)+Q(5)).*(Cp(2) - s(2).*y5(2));
    UR5(3) = (Q(4)+Q(5)).*(Cp(3) - s(3).*y5(3));
    UR5(4) = (Q(4)+Q(5)).*(Cp(4) - s(4).*y5(4));
    UR5(5) = (Q(4)+Q(5)).*(Cp(5) - s(5).*y5(5));
    UR5(6) = (Q(4)+Q(5)).*(Cp(6) - s(6).*y5(6));
    UR5(7) = (Q(4)+Q(5)).*(Cp(7) - s(7).*y5(7));
    UR5(8) = (Q(4)+Q(5)).*(Cp(8) - s(8).*y5(8));
    UR5(9) = (Q(4)+Q(5)).*(Cp(9) - s(9).*y5(9));
    UR5(10:22) = 0.0;

%Uptake-release rates for the Adipose (A) compartment - y6(x)  
    %Blood/tissue partition coefficients
    if sex == 0  % male
        s(1) = 1.9269;    %GLC
        s(2) = 0.1838;    %PYR
        s(3) = 1.0434;    %LAC
        s(4) = 0.0;       %ALA
        s(5) = 1.5429;    %GLR
        s(6) = 2.1861;    %FFA
        s(7) = 0.0009439*hyperTG; %TG
        s(8) = 249.9761;  %O2
        s(9) = 1.4640;    %CO2
    else  % female
        s(1) = 1.8915;    %GLC
        s(2) = 0.1838;    %PYR
        s(3) = 1.0434;    %LAC
        s(4) = 0.0;       %ALA
        s(5) = 1.5429;    %GLR
        s(6) = 2.3616;    %FFA   vs ----> male value = 2.1861
        s(7) = 8.8328e-4*hyperTG; %TGL vs ----> male value = 0.0009439; 
        s(8) = 249.9761; %O2
        s(9) = 1.4560;  %CO2
    end
    s(10:22) = 0.0;

    %Uptake/Release rates
    UR6(1) = Q(6) * (Ca(1) - s(1)*y6(1));
    UR6(2) = Q(6) * (Ca(2) - s(2)*y6(2));
    UR6(3) = Q(6) * (Ca(3) - s(3)*y6(3));
    UR6(4) = 0.0;
    UR6(5) = Q(6) * (Ca(5) - s(5)*y6(5));
    UR6(6) = Q(6) * (Ca(6) - s(6)*y6(6));
    UR6(7) = Q(6) * (Ca(7) - s(7)*y6(7));
    UR6(8) = Q(6) * (Ca(8) - s(8)*y6(8));
    UR6(9) = Q(6) * (Ca(9) - s(9)*y6(9));
    UR6(10:22) = 0.0;

%Uptake-release rates for the Other tissues (O) compartment - y7(x) 
    UR7(1) = 0.062;
    UR7(2) = -0.005;
    UR7(3) = -0.142 + max(0,Q(7)*(Ca(3)-2)); 
    UR7(4) = -0.280;
    UR7(5) = 0.0 + max(0,Q(7)*(Ca(5)-0.33)); 
    UR7(6) = 0.04; 
    UR7(7) = 0.0;
    UR7(8) = 2.146;
    UR7(9) = -1.572; 
    UR7(10:22) = 0.0; 

    %Call tissue-specific functions
    [F1,P1,U1] = brain(t,y1,MR);
    P(1,1:22)=P1(1:22); U(1,1:22)=U1(1:22);
    
    [F2,P2,U2] = heart(t,tstart,tend,trecovery,y2,y,WR,epi,epi0,MR,RMR_H,c0,d0);
    P(2,1:22)=P2(1:22); U(2,1:22)=U2(1:22); 
    
    [F3,P3,U3] = muscle(t,tstart,tend,trecovery,y3,y,sex,WR,MR,RMR_M,epi,epi0,c0,d0);
    P(3,1:22)=P3(1:22); U(3,1:22)=U3(1:22); 
    
    [F4,P4,U4] = gi(t,tstart,tend,trecovery,y4,y,WR,MR,GIR,epi,epi0,c0,d0);
    P(4,1:22)=P4(1:22); U(4,1:22)=U4(1:22); 
    
    [F5,P5,U5] = liver(t,tstart,tend,trecovery,y5,y,sex,WR,GIR,MR,c0,d0);
    P(5,1:22)=P5(1:22); U(5,1:22)=U5(1:22); 
    
    [F6,P6,U6] = adipose(t,tstart,tend,trecovery,y6,y,sex,WR,MR,GIR,epi,epi0,c0,d0);
    P(6,1:22)=P6(1:22); U(6,1:22)=U6(1:22);
    
    %Tissue volumes, convert Kg to L using tissue density
    V(1)=WB*1.04;
    V(2)=WH*1.05;
    V(3)=WM*1.05;
    V(4)=WG*1.06;
    V(5)=WL*1.08;
    V(6)=WA*0.923;
    V(7)=WO*1.05;
    if sex == 0  % male
        Vb = 5.5; % Blood volume
    else  % female
        Vb = 5; % Blood volume
    end

    %% store outputs
    UR(1,1:22) = UR1;
    UR(2,1:22) = UR2;
    UR(3,1:22) = UR3;
    UR(4,1:22) = UR4;
    UR(5,1:22) = UR5;
    UR(6,1:22) = UR6;
    uprel = [UR;UR7];
    bflow = Q;
    F = [F1; F2; F3; F4; F5; F6]; 

    %% Begin odes
    % Dietary input fluxes
    v_feed_Glc_B = zeros(1,nmealsTot);
    v_feed_TG_B  = zeros(1,nmealsTot);
    if (tot_Glc_B>0 || tot_TG_B>0)

        % time delays for gastric emptying
        t_delay_Glc_B = 30+hill(tot_Glc_B,60,754,2); % Frayn, page 174, GLC peaks 30-60 mins after a meal
        t_delay_TG_B  = 60+hill(tot_TG_B,120,129,2); % Frayn, page 185, TG peaks 3-4 hours after a meal 
        t1=t;
        
        % GI to blood influx functions
        i = 1:1:nmealsTot;
        v_feed_Glc_B(1,:) = influx_Glc_B(t1,i,lagm,tot_Glc_B,t_delay_Glc_B);  
        v_feed_TG_B(1,:) = influx_TG_B(t1,i,lagm,tot_TG_B,t_delay_TG_B); 
    else
        v_feed_Glc_B(1,:) = 0;
        v_feed_TG_B(1,:) = 0;
    end
    % v_feed_[] are row vectors, thus summing over the columns produces a single value
    % We add the different fluxes in case the meals are so close to each other 
    % that gastric emptying is still in progress when a new meal arrives, 
    % so we must consider the combined flux.
    v_feed_Glc_B = sum(max(0,v_feed_Glc_B),2);
    v_feed_TG_B  = sum(max(0,v_feed_TG_B),2);

    % Blood compartment
    netflux     = -(UR1+UR2+UR3+UR4+UR5+UR6+UR7);
    dydt(1,1)   = (netflux(1)+v_feed_Glc_B)/Vb;   % add glucose uptake from meal, if any 
    dydt(2:6,1) = netflux(2:6)/Vb;
    dydt(7,1)   = (netflux(7)+v_feed_TG_B)/Vb;    % add TG uptake from meal, if any 
    
    % Tissues B,H,M,G,L,A
    for k=1:6 
        dydt((1+22*k):(22*(k+1)),1) = (P(k,:)-U(k,:)+UR(k,:))/V(k);
    end  

    % Insulin and glucagon
    if  (t>tstart) && (t<tend+trecovery) && (WR>0) 
        dydt(155,1)=y(155)*(phi*(h-k1*(y(155)-glu0)-k2*(y(156)-ins0))-D); % Glucagon 
        dydt(156,1)=y(156)*(psi*(h-k3*(y(155)-glu0)-k4*(y(156)-ins0)-k5*(epi-epi0))-D); % Insulin
    else
        dydt(155,1) = k_glu_b+phi*(hill_inh(Ca(1),Vmax_glc_glu,Ki_glc_glu,n_glc_glu)*hill_inh(y(156),1,Ki_ins_glu,n_ins_glu))-D*y(155); % Glucagon
        dydt(156,1) = k_ins_b+psi*(hill(Ca(1),Vmax_glc_ins,Km_glc_ins,n_glc_ins))-D*y(156); % Insulin 
    end

    % end odes
end

% Activation fucntion
function val = hill(x, Vmax, Km, n)
val=Vmax*x^n/(Km^n+x^n);
end

% Inhibition fucntion
function val = hill_inh(x, Vmax, Ki, n)
val=Vmax*Ki^n/(Ki^n+x^n);
end

% uptake/release sensitivity function
% allows organs to store more substrate when available aplenty
function val = sens(x,ic,a,slope) % ic is initial concentration of substrate x
    val = 1-a/(1+exp(-slope*(x-ic)));
end

% IF functions  
function val = influx_Glc_B(t,k,lag,tot_Glc_B,t_delay_Glc_B)  
    val = tot_Glc_B .* (t-(k-1)*lag)/(t_delay_Glc_B^2) .* exp(-(t-(k-1)*lag).^2./(2*t_delay_Glc_B^2));
end

function val = influx_TG_B(t,k,lag,tot_TG_B,t_delay_TG_B)  
    val = tot_TG_B  .* (t-(k-1)*lag)/(t_delay_TG_B^2)  .* exp(-((t-(k-1)*lag)-200).^2./(2*t_delay_TG_B^2));
end   