%% [5] - LIVER COMPARTMENT
function [F5,P5,U5] = liver(t,tstart,tend,trecovery,y5,y,sex,WR,GIR,MR,c0,d0)
    
    %Initialize
    ptemp = zeros(22,1);
    utemp = zeros(22,1);

    lambda4  = c0(33); alpha4  = c0(34);
    lambda5  = c0(35); alpha5  = c0(36);
    lambda6  = c0(37); alpha6  = c0(38);
    lambda8  = c0(39); alpha8  = c0(40);
    lambda15 = c0(41); alpha15 = c0(42);
    
    %Vmax values
    Vmax(1) = 0.765;
    Vmax(2) = 0.68;
    Vmax(3) = 5.44;
    Vmax(4) = 7.44;
    Vmax(5) = 2.08;
    Vmax(6) = 1.80;
    Vmax(7) = 0.40;
    Vmax(8) = 3.84;
    Vmax(9) = 0.84;
    Vmax(10) = 1.92;
    Vmax(11) = 0.576;
    Vmax(12) = 0.0;
    Vmax(13) = 0.444;
    Vmax(14) = 0.0;
    Vmax(15) = 0.64;
    Vmax(16) = 0.01;
    Vmax(17) = 1.088;
    Vmax(18) = 0.896;
    Vmax(19) = 0.008;
    Vmax(20) = 0.8;
    Vmax(21) = 15.62;
    Vmax(22) = 22.18;
    Vmax(23) = 0.0;
    Vmax(24) = 0.0;
    Vmax(25) = MR(5)*2.0;
    
    if (t>tstart) && (t<tend+trecovery) && (WR>0)
        Vmax(4)=Vmax(4)*(1+lambda4*GIR^2/(alpha4+GIR^2));
        Vmax(5)=Vmax(5)*(1+lambda5*GIR^2/(alpha5+GIR^2));
        Vmax(6)=Vmax(6)*(1+lambda6*GIR^2/(alpha6+GIR^2));
        Vmax(8)=Vmax(8)*(1+lambda8*GIR^2/(alpha8+GIR^2));
        Vmax(15)=Vmax(15)*(1+lambda15*GIR^2/(alpha15+GIR^2));
    else %if (WR==0)
        % glucose metabolism
        Vmax_ins_L  = d0(34);
        Km_ins_L    = d0(35);
        n_ins_L     = d0(36);
        beta_L      = hill(y(156),Vmax_ins_L,Km_ins_L,n_ins_L);

        Vmax_ins_L2 = d0(37);
        Km_ins_L2   = d0(38);
        n_ins_L2    = d0(39);
        beta_L2     = hill(y(156),Vmax_ins_L2,Km_ins_L2,n_ins_L2);

        Vmax_ins_L3 = d0(40);
        Km_ins_L3   = d0(41);
        n_ins_L3    = d0(42);
        beta_L3     = hill(y(156),Vmax_ins_L3,Km_ins_L3,n_ins_L3);

        Vmax_ins_inh_L = d0(43);
        Km_ins_inh_L   = d0(44);
        n_ins_inh_L    = d0(45);
        gamma_L        = hill(y(156),Vmax_ins_inh_L,Km_ins_inh_L,n_ins_inh_L);
        
        % fat metabolism     
        Vmax_ins_L_fat = d0(46);  
        Km_ins_L_fat   = d0(47);   
        n_ins_L_fat    = d0(48);
        beta_L_fat     = hill(GIR^2,Vmax_ins_L_fat,Km_ins_L_fat,n_ins_L_fat);

        Vmax_ins_L_fat2= d0(49); 
        Km_ins_L_fat2  = d0(50);  
        n_ins_L_fat2   = d0(51);
        beta_L_fat2    = hill(GIR^2,Vmax_ins_L_fat2,Km_ins_L_fat2,n_ins_L_fat2);

        Vmax_ins_inh_L_fat = d0(52); 
        Km_ins_inh_L_fat   = d0(53);  
        n_ins_inh_L_fat    = d0(54);
        gamma_L_fat        = hill(GIR^2,Vmax_ins_inh_L_fat,Km_ins_inh_L_fat,n_ins_inh_L_fat);

        Vmax(1)     = Vmax(1)*(1+beta_L2);
        Vmax(2)     = Vmax(2)*(1+beta_L3); 
        Vmax(3)     = Vmax(3)*(1+beta_L3);
        Vmax(7)     = Vmax(7)*(1+beta_L2*4-gamma_L_fat); % Note: beta_L2*4 

        Vmax(16)    = Vmax(16)*(1+beta_L_fat2);
        Vmax(10)    = Vmax(10)*(1+beta_L_fat); 
        Vmax(9)     = Vmax(9)*(1-gamma_L_fat); 

        Vmax(18)    = Vmax(18)*(1+beta_L);      
        Vmax(20)    = Vmax(20)*(1+beta_L);

        Vmax(4)     = Vmax(4)*(1+beta_L_fat-gamma_L);
        Vmax(5)     = Vmax(5)*(1+beta_L_fat-gamma_L);
        Vmax(6)     = Vmax(6)*(1+beta_L_fat-gamma_L);
        Vmax(8)     = Vmax(8)*(1+beta_L_fat-gamma_L); 
        Vmax(15)    = Vmax(15)*(1+beta_L_fat);
        Vmax(17)    = Vmax(17)*(1+beta_L_fat);  
   end 

    %Km values
    Km(1) = 10.0;  %GLC
    Km(2) = 0.37;  %PYR
    Km(3) = 0.82;  %LAC
    Km(4) = 0.23;  %ALA
    Km(5) = 0.07;  %GLR
    Km(6) = 0.57;  %FFA
    Km(7) = 2.93;  %TGL
    Km(8) = 7.0e-4;%O2
    Km(9) = 15.43; %CO2
    Km(10) = 0.2;  %G6P
    Km(11) = 417;  %GLY
    Km(12) = 0.11; %GAP
    Km(13) = 0.24; %GRP
    Km(14) = 0.035;%ACoA
    Km(15) = 0.14; %CoA
    Km(16) = 0.45; %NAD+
    Km(17) = 0.05; %NADH
    Km(18) = 2.74; %ATP
    Km(19) = 1.22; %ADP
    Km(20) = 4.60; %Pi
    Km(21) = 0.1;  %PCR
    Km(22) = 0.1;  %CR
    Km(23) = Km(16)/Km(17); %NAD/NADH
    Km(24) = Km(17)/Km(16); %NADH/NAD
    Km(25) = Km(18)/Km(19); %ATP/ADP
    Km(26) = Km(19)/Km(18); %ADP/ATP 

    [F5,P5,U5] = callflux(y5,Vmax,Km,ptemp,utemp,5);

end
% Activation fucntion
function val = hill(x, Vmax, Km, n)
val=Vmax*x^n/(Km^n+x^n);
end