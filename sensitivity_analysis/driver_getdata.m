%% Save the data
if isempty(my_field) == 1
    Y1(:,1:22)=Y(:,23:44);   %brain
    Y2(:,1:22)=Y(:,45:66);   %heart
    Y3(:,1:22)=Y(:,67:88);   %muscle
    Y4(:,1:22)=Y(:,89:110);  %GI
    Y5(:,1:22)=Y(:,111:132); %liver
    Y6(:,1:22)=Y(:,133:154); %adipose
    arterial = Y(:,1:22);    %arterial substrates
    hormone = Y(:,155:156);  %glucagon and insulin
    
    [dydt,bflowi,MRi,P,U,uprel,F,epi,v_feed_Glc_B,v_feed_TG_B] = cellfun(@(t,y) IFMod(t,y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend,sm,sa),num2cell(T),num2cell(Y,2),'uni',0);
    
    bflow = cell2mat(bflowi); %blood flow
    MR = cell2mat(MRi);       %metabolic rate of ATP -> ADP
    epinephrine = cell2mat(epi);
    v_feed_Glc = cell2mat(v_feed_Glc_B);
    v_feed_TG = cell2mat(v_feed_TG_B);
    
    
    %production
    P  = cell2mat(P);
    P1 = P(1:6:end,:);
    P2 = P(2:6:end,:);
    P3 = P(3:6:end,:);
    P4 = P(4:6:end,:);
    P5 = P(5:6:end,:);
    P6 = P(6:6:end,:);
    P = {P1 P2 P3 P4 P5 P6}; % override P

    %utilization
    U  = cell2mat(U);
    U1 = U(1:6:end,:);
    U2 = U(2:6:end,:);
    U3 = U(3:6:end,:);
    U4 = U(4:6:end,:);
    U5 = U(5:6:end,:);
    U6 = U(6:6:end,:);
    U = {U1 U2 U3 U4 U5 U6}; % override U
    
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
else
    Y1.(my_field)(:,1:22)=Y.(my_field)(:,23:44);   %brain
    Y2.(my_field)(:,1:22)=Y.(my_field)(:,45:66);   %heart
    Y3.(my_field)(:,1:22)=Y.(my_field)(:,67:88);   %muscle
    Y4.(my_field)(:,1:22)=Y.(my_field)(:,89:110);  %GI
    Y5.(my_field)(:,1:22)=Y.(my_field)(:,111:132); %liver
    Y6.(my_field)(:,1:22)=Y.(my_field)(:,133:154); %adipose
    arterial.(my_field) = Y.(my_field)(:,1:22);    %arterial substrates
    hormone.(my_field) = Y.(my_field)(:,155:156);  %glucagon and insulin
    
    [dydt,bflowi,MRi,P,U,uprel,F,epi,v_feed_Glc_B,v_feed_TG_B] = cellfun(@(t,y) IFMod(t,y,sex,WR,ins0,glu0,c0,d0,tot_Glc_B,tot_TG_B,nmealsTot,lagm,tstart,tend,sm,sa),num2cell(T.(my_field)),num2cell(Y.(my_field),2),'uni',0);
    
    bflow.(my_field) = cell2mat(bflowi); %blood flow
    MR.(my_field) = cell2mat(MRi);       %metabolic rate of ATP -> ADP
    epinephrine.(my_field) = cell2mat(epi);
    v_feed_Glc.(my_field) = cell2mat(v_feed_Glc_B);
    v_feed_TG.(my_field) = cell2mat(v_feed_TG_B);
    
    
    %production
    P  = cell2mat(P);
    P1.(my_field) = P(1:6:end,:);
    P2.(my_field) = P(2:6:end,:);
    P3.(my_field) = P(3:6:end,:);
    P4.(my_field) = P(4:6:end,:);
    P5.(my_field) = P(5:6:end,:);
    P6.(my_field) = P(6:6:end,:);
    
    %utilization
    U  = cell2mat(U);
    U1.(my_field) = U(1:6:end,:);
    U2.(my_field) = U(2:6:end,:);
    U3.(my_field) = U(3:6:end,:);
    U4.(my_field) = U(4:6:end,:);
    U5.(my_field) = U(5:6:end,:);
    U6.(my_field) = U(6:6:end,:);
    
    %uptake-release rates
    uprel = cell2mat(uprel);
    UR1.(my_field) = uprel(1:7:end,:);
    UR2.(my_field) = uprel(2:7:end,:);
    UR3.(my_field) = uprel(3:7:end,:);
    UR4.(my_field) = uprel(4:7:end,:);
    UR5.(my_field) = uprel(5:7:end,:);
    UR6.(my_field) = uprel(6:7:end,:);
    UR7.(my_field) = uprel(7:7:end,:); % other tissues
    
    %metabolic fluxes
    F  = cell2mat(F);
    F1.(my_field) = F(1:6:end,:);
    F2.(my_field) = F(2:6:end,:);
    F3.(my_field) = F(3:6:end,:);
    F4.(my_field) = F(4:6:end,:);
    F5.(my_field) = F(5:6:end,:);
    F6.(my_field) = F(6:6:end,:);
end
 