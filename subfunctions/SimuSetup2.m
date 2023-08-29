function [Feeder,Ybus,SlackBus,Base,puTrue] = SimuSetup2(SystemN,iScenario)
switch SystemN
    case 1
        load('IEEE37.mat')
        SlackBus=1;
        Base.Sbase = abs(P_node(SlackBus,1));
        Base.Vbase = abs(V_node(SlackBus,1));
        Base.Ibase=Base.Sbase/Base.Vbase;
        Base.Ybase =  Base.Sbase/Base.Vbase^2;
     case 2
%         load('IEEE123_2.mat') %% ADDED 2022
        load(strcat('IEEE123_',num2str(iScenario),'.mat')) %% ADDED 2022
        SlackBus=1;  %% Modified 8.20 117 for 150 node 1 for regulator, 124 for no regulator
        Base.Sbase = abs(P_node(SlackBus,1));
        Base.Vbase = abs(V_node(SlackBus,1));
        Base.Ibase=Base.Sbase/Base.Vbase;
        Base.Ybase =  Base.Sbase/Base.Vbase^2;
%         V_node(end,:) = V_node(end,:)*9/10;

%         Base.Sbase = 1;
%         Base.Vbase = 1;
%         Base.Ibase=1;
%         Base.Ybase =  1;

        case 3
        load('SHE215-3_temp.mat') %% ADDED 2022
        SlackBus=67;
        Base.Sbase = abs(P_node(SlackBus,1));
        Base.Vbase = abs(V_node(SlackBus,1));
        Base.Ibase=Base.Sbase/Base.Vbase;
        Base.Ybase =  Base.Sbase/Base.Vbase^2;
end


puTrue.VOL=V_node/Base.Vbase;%% complex voltage
puTrue.PN=real(P_node)/Base.Sbase;
puTrue.QN=imag(P_node)/Base.Sbase;
puTrue.IN=abs(I_node)/Base.Ibase;
puTrue.PL=real(P_line)/Base.Sbase;
puTrue.QL=imag(P_line)/Base.Sbase;
puTrue.IL=abs(I_line)/Base.Ibase;
puTrue.Ve=real(puTrue.VOL);
puTrue.Vf=imag(puTrue.VOL); 
puTrue.VN=abs(puTrue.VOL);
puTrue.AN=angle(puTrue.VOL);

% data preprocessing
puTrue.PN(abs(puTrue.PN)<1e-8)=0;
puTrue.QN(abs(puTrue.QN)<1e-8)=0;
puTrue.PL(abs(puTrue.PL)<1e-8)=0;
puTrue.QL(abs(puTrue.QL)<1e-8)=0;
puTrue.IN(abs(puTrue.IN)<1e-8)=0;
puTrue.IL(abs(puTrue.IL)<1e-8)=0;




