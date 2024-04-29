%% [6] - ADIPOSE TISSUE COMPARTMENT
function [F6,P6,U6] = adipose(t,tstart,tend,trecovery,y6,y,sex,WR,MR,GIR,epi,epi0,c0,d0)
   
    
    %Initialize
    ptemp = zeros(22,1);
    utemp = zeros(22,1);

    %Assign parameters to estimate 
    lambda19_gir = c0(43); alpha19_gir = c0(44);
    lambda19_epi = c0(45); alpha19_epi = c0(46);    
    if sex==0 % male
        ffa_adipose = 1;
    else % female
        ffa_adipose =  c0(49); %=1.2 >1
    end

    %Vmax values
    Vmax(1) = 0.079;
    Vmax(2) = 0.152;
    Vmax(3) = 0.896;
    Vmax(4) = 0.0;
    Vmax(5) = 0.0;
    Vmax(6) = 0.0;
    Vmax(7) = 0.0;
    Vmax(8) = 0.0;
    Vmax(9) = 0.144;
    Vmax(10) = 0.04;
    Vmax(11) = 0.0;
    Vmax(12) = 0.08;
    Vmax(13) = 0.0;
    Vmax(14) = 0.0;
    Vmax(15) = 0.0;
    Vmax(16) = 0.24;
    Vmax(17) = 0.16;
    Vmax(18) = 0.0;
    Vmax(19) = 0.19*ffa_adipose;
    Vmax(20) = 0.48;
    Vmax(21) = 1.28;
    Vmax(22) = 2.05;
    Vmax(23) = 0.0;
    Vmax(24) = 0.0;
    Vmax(25) = MR(6)*2.0;

      if(t>tstart) && (t<tend+trecovery) && (WR>0)
            Vmax(19) = Vmax(19)*(1+(lambda19_gir*GIR^2/(alpha19_gir+GIR^2)) + (lambda19_epi*(epi-epi0)^2/(alpha19_epi+(epi-epi0)^2))); 
      else %if (WR==0)
        Vmax_ins_A = d0(55);
        Km_ins_A   = d0(56);
        n_ins_A    = d0(57); 
        beta_A     = hill(y(156),Vmax_ins_A,Km_ins_A,n_ins_A);

        Vmax_ins_A2 = d0(58);
        Km_ins_A2   = d0(59);
        n_ins_A2    = d0(60); 
        beta_A2     = hill(y(156),Vmax_ins_A2,Km_ins_A2,n_ins_A2);

        Vmax_ins_inh_A = d0(61); 
        Km_ins_inh_A   = d0(62); 
        n_ins_inh_A    = d0(63); 
        gamma_A        = hill(y(156),Vmax_ins_inh_A,Km_ins_inh_A,n_ins_inh_A);
       
%         Vmax_ins_A_fat = d0(64);
%         Km_ins_A_fat   = d0(65);
%         n_ins_A_fat    = d0(66); 
%         beta_A_fat     = hill(GIR^2,Vmax_ins_A_fat,Km_ins_A_fat,n_ins_A_fat);
       
%-----------------------------------------------------          
        %de novo lipogenesis via glycolysis
        Vmax(1)=Vmax(1)*(1+beta_A);
        Vmax(2)=Vmax(2)*(1+beta_A);
        Vmax(3)=Vmax(3)*(1+beta_A);
        Vmax(16)=Vmax(16)*(1+beta_A);
        Vmax(20)=Vmax(20)*(1+beta_A2);
        % lipolysis
        Vmax(19)=Vmax(19)*(1-gamma_A);
%         % TG synthesis
%         Vmax(20)=Vmax(20)*(1+beta_A_fat);
     end

    %Km values^
    Km(1) = 0.1;   %GLC
    Km(2) = 0.37;  %PYR
    Km(3) = 0.82;  %LAC
    Km(4) = 0.1;   %ALA
    Km(5) = 0.22;  %GLR
    Km(6) = 0.57;  %FFA
    Km(7) = 990.0; %TGL
    Km(8) = 7.0e-4;%O2
    Km(9) = 15.43; %CO2
    Km(10) = 0.2;  %G6P
    Km(11) = 0.1;  %GLY
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

    [F6,P6,U6] = callflux(y6,Vmax,Km,ptemp,utemp,6);

end
% Activation fucntion
function val = hill(x, Vmax, Km, n)
val=Vmax*x^n/(Km^n+x^n);
end