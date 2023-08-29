%%% Main_DSE
%%% 3 phase unbalanced distribution state estimation toolbox
%%% Author: Jiaojiao Dong
function  seResult = DSSE_demo7(SCADA_VN, measure_VN, NodeList, Ybig, Pbig, Qbig, vckt, Lines, zeroinj_V)
%% SCADA_VN is the location of Vn (size N*2, column is location, row is phase), VN is the measurements of Vn (size N*1)
%% Ybig is the system Y, Lines, NodeList is the name of vckt in OpenDSS
%% Pbig, Qbig is the pseduo measurements
%% vckt is the rough estimate of voltage, column variable

% ##system('taskkill /F /IM EXCEL.EXE'); 
% clear;clc;warning off;
addpath('subfunctions');global Base;
% tic

%% ADDED 2022
addpath('PF')
nScenario = 1; %% ADDED 8.9
% ratio = 1.05;
isPlotPOW = 0;
isPlotFeeder = 0;
capNumber = 4;  %% ADDED 8.9
for iScenario = 1: nScenario
%     iScenario = 2;
    %% ADDED 8.9
    load('loadPVInfo.mat');
    
    PV_temp = P_PV(iScenario,Loc_PV);
    PV2 = [Loc_PV, -PV_temp.'];
    ratio = ratio_load(iScenario);
%     ratio = 1;
    
%     SEData('IEEE123',iScenario,ratio,P_PV,capNumber)
    %% END ADDED 8.9
%     processOpenDSS_v3_2(Pbig, Qbig, vckt, Ybig, Lines, NodeList)
    processOpenDSS_v3_1(Pbig, Qbig, Ybig, Lines, NodeList, zeroinj_V)
    nSimu = 1;
    for iSimu = 1:nSimu
%% END ADDED 2022    
        % %% PARAMETERS SETTING
        SystemN         =    2;   % 1 for 37 node system, 2 for IEEE 123 node 
        Flag.ini        =    1;   % 1 for (1,0) initial value. 2 for true value, 3 for slack node voltage
        Flag.noise      =    2;   % 1 and 11 for gauss   noise , 2 for no noise, 3 and 31 for uniform noise 4 for extreme noise
        Flag.zeroinj    =    1;   % 1 for considering zero injection, 0 for not considering zero injection
        Flag.error      =    1;   % 1 for absolute error, 2 for relative error
        iterMax         =    100; % Max Iteration in Newton method
        tol             =    1e-4;% tolerance of the Newton Method
%         [Feeder,Ybus,SlackBus,Base,puTrue] = SimuSetup(SystemN,iScenario);%% ADDED 2022 
        [Feeder,Ybus,SlackBus,Base,puTrue] = SimuSetup(SystemN);%% ADDED 2022 
        
    
        %% Sensor Error
        SensorError.PN      = 0.001;
        SensorError.QN      = 0.001;
        SensorError.VN      = 0.001;
        SensorError.AN      = 0.001; %% ADDED 2022
        SensorError.PL      = 0.001;
        SensorError.QL      = 0.001;
    
        %% Numbers
        Num.Node                        = Feeder.NumN;
        Num.Line                        = Feeder.NumL;
        Num.StateVar                    = 2; % for each node
        [Num.Zeroinj, Loc.Zeroinj]      = Zero_inj(puTrue.IN); %Zero Injection Num and location
        [Num.Zerovol, Loc.Zerovol]      = Zero_inj(puTrue.VN); 
%         Loc.Line                        = Feeder.Topology(:,1:2);
    
        %% Creat measure with noise
        Noise = Add_noise(Flag,Num,SensorError);%% ADDED 2022
        
        %% Index of measure locations 
%         MeaIdx.PN = [1:Num.Node]';
        MeaIdx.PN = [1:size(NodeList,1)]';
        MeaIdx.QN = MeaIdx.PN;
    
%         SlackIdx = [SlackBus:Num.Node:3*Num.Node]';
% 
%         MeaIdx.PN = setdiff(MeaIdx.PN,SlackIdx);
%         MeaIdx.QN = setdiff(MeaIdx.QN,SlackIdx);
%% ADDED 2022
        PMU = []; 
        POW = [];  
%         SCADA_VN = [37 1; 37 2; 37 3;
%                  96 1; 96 2; 96 3;
%                  83 1; 83 2; 83 3;];


        [VN_POW, AN_POW] =  processPowData();
        

        % if isPlotFeeder
        %     figure('Renderer', 'painters', 'Position', [100 100 900 500])
        %     hold on; box on
        %     PlotFeeder(Feeder)
        %     set(gca, 'Fontname', 'Times New Roman', 'Fontsize', 12);

        %     isPlotFeeder = 0;
        % end        

        % if isPlotPOW
        %     load("powData.mat")
        %     figure('Renderer', 'painters', 'Position', [100 100 600 300])
        %     hold on; box on
        %     plot(t,powData60A,'b-o', 'MarkerFaceColor','k', 'MarkerEdgeColor','k','MarkerSize',2)
        %     xlabel('Time (second)')
        %     ylabel('Voltage (p.u.)')
        %     set(gca, 'Fontname', 'Times New Roman', 'Fontsize', 12);

        %     isPlotPOW = 0;
        % end
        
    
%         MeaIdx.VN = pha2idx([PMU;POW;SCADA_VN],Num.Node);
        MeaIdx.VN = [PMU;POW;SCADA_VN];
        MeaIdx.AN =pha2idx([PMU;POW],Num.Node);
    
    %% END ADDED 2022

        MeaIdx.PL = [];
        MeaIdx.QL = [];
        MeaPha.PN = idx2pha(MeaIdx.PN,Num.Node);
        MeaPha.QN = idx2pha(MeaIdx.QN,Num.Node);
        MeaPha.VN = idx2pha(MeaIdx.VN,Num.Node);
        MeaPha.AN = idx2pha(MeaIdx.AN,Num.Node);%% ADDED 2022
        MeaPha.PL = idx2pha(MeaIdx.PL,Num.Line);
        MeaPha.QL = idx2pha(MeaIdx.QL,Num.Line);
    
        Weight = WeightMatrix(MeaIdx,SensorError);%% ADDED 2022
        Weight = sparse(Weight);
    
        if Flag.noise == 2; Weight = eye(size(Weight,1)); end
        %% creat measure with noise
        Mea.PN  = puTrue.PN .* Noise.PN; 
        Mea.QN  = puTrue.QN .* Noise.QN;
        Mea.VN  = puTrue.VN.^2 .* Noise.VN;
        Mea.AN  = puTrue.AN.* Noise.AN;%% ADDED 2022

%         Mea.VN(60,:) = VN_POW;%% ADDED 2022
%         Mea.AN(60,:) = AN_POW;%% ADDED 2022

        Mea.PL  = puTrue.PL .* Noise.PL;
        Mea.QL  = puTrue.QL .* Noise.QL;

%         Mea.PN = reshape(Mea.PN,3*Num.Node,1);
%         Mea.QN = reshape(Mea.QN,3*Num.Node,1);
%         Mea.VN = reshape(Mea.VN,3*Num.Node,1);
%         Mea.AN = reshape(Mea.AN,3*Num.Node,1);%% ADDED 2022
%         Mea.PL = reshape(Mea.PL,3*Num.Line,1);
%         Mea.QL = reshape(Mea.QL,3*Num.Line,1);
        measure = [Mea.PN(MeaIdx.PN); Mea.QN(MeaIdx.QN);measure_VN;Mea.AN(MeaIdx.AN);...
          Mea.PL(MeaIdx.PL);Mea.QL(MeaIdx.QL)];%% ADDED 2022
    
        %% Newton
        display(strcat('started solving DSSE for t',num2str(iScenario),' in the time series.'))
        display('==============Voltage==============Angle (degree)============')
        display('==== Phase A  Phase B  Phase C == Phase A  Phase B  Phase C==')        
        [StateVar,~] = Ini_val(Flag,Num,puTrue,SlackBus);
        G = real(Ybus)/Base.Ybase; B = imag(Ybus)/Base.Ybase;
        tol1 = [];
        for i = 1:iterMax
          % Jacobian matrix
          tic
          
          joc_vn          = vnjoc(Num,StateVar,MeaPha);
          [joc_pn,joc_qn] = pnjoc(Num,StateVar,MeaPha,B,G);
%           joc_pn([1,1+132,1+2*132],:) = [];
%           joc_qn([1,1+132,1+2*132],:) = [];

          
%           joc_an          = anjoc(Num,StateVar,MeaPha);%% ADDED 2022
%           [joc_pl,joc_ql] = pljoc(Num,StateVar,MeaPha,B,G,Loc);
%           [joc_c,joc_d] = cjoc(Num,B,G,Loc);
          joc_an = [];
          joc_pl =[];
          joc_ql =[];
          joc_c = [];
          joc_d = [];
          

          H = [joc_pn;joc_qn;joc_vn;joc_an; joc_pl;joc_ql];%% ADDED 2022
          C = [joc_c;joc_d];
          H(:,2*SlackBus:2*Num.Node:end)=0;   
%           C(:,2*SlackBus:2*Num.Node:end)=0;
          H(:,2*SlackBus-1:2*Num.Node:end)=0; 
%           C(:,2*SlackBus-1:2*Num.Node:end)=0;
          % MeaEqu
          [equ_pn,equ_qn] = pnequ(Num,StateVar,MeaPha,B,G);
          equ_vn = vnequ(Num,StateVar,MeaPha);
          
          
%           equ_an = anequ(Num,StateVar,MeaPha);%% ADDED 2022
%           [equ_pl,equ_ql] = plequ(Num,StateVar,MeaPha,B,G,Loc);
%           [equ_c,equ_d] = cequ(Num,StateVar,B,G,Loc);
          equ_an = [];
          equ_pl = [];
          equ_ql = [];
          equ_c = [];
          equ_d = [];

          
          
          MeaEqu = [equ_pn;equ_qn;equ_vn;equ_an; equ_pl;equ_ql];%% ADDED 2022
          delta_z = measure-MeaEqu;
          delta_c = [-equ_c;-equ_d];

          clear joc_pn joc_qn joc_vn joc_an 
          clear joc_pl joc_ql joc_c joc_d  C 
          clear equ_pn equ_qn equ_vn equ_an
          clear equ_pl equ_ql equ_c equ_d MeaEqu delta_c
          
          G2 = H'/Weight*delta_z;
          G3 = H'/Weight*H;    
 
          %% ADDED 10.9.2022
          [zeroRow,zeroCol] = computeZero(G3);
          sizeG3 = size(G3,1);
          G3(zeroRow,:)=[];
          G3(:,zeroCol)=[];

          %% END ADDED 10.9.2022

          clear H delta_z
          
          Ginv = pinv(G3);
          
          %% ADDED 10.9.2022
          nonindex = setdiff([1:sizeG3],zeroRow);
          G3_temp = zeros(sizeG3, sizeG3);
          G3_temp(nonindex,nonindex) = Ginv;
          Ginv = G3_temp;

          %% END ADDED 10.9.2022
          

          
%           C_temp1 = C*Ginv*C';
%           
%           %% ADDED 10.9.2022
%           [zeroRow,zeroCol] = computeZero(C_temp1);
%           sizeG3 = size(C_temp1,1);
%           C_temp1(zeroRow,:)=[];
%           C_temp1(:,zeroCol)=[];
% 
%           %% END ADDED 10.9.2022
%           
%           C_temp2 = C_temp1;
%           
%                     %% ADDED 10.9.2022
%           nonindex = setdiff([1:sizeG3],zeroRow);
%           C_temp3 = zeros(sizeG3, sizeG3);
%           C_temp3(nonindex,nonindex) = C_temp2;
%           C_temp2 = C_temp3;

%           deltaX = Ginv * G2 - Flag.zeroinj * Ginv*C'*(C_temp2*(C*Ginv*G2-delta_c));      
          deltaX = Ginv * G2; 
          tol_temp = max(abs(deltaX));
          if tol_temp < tol
              tol1 = [tol1; tol_temp];
              break;
          end
          StateVar = StateVar + deltaX;
          [i tol_temp]
          tol1 = [tol1; tol_temp];

          clear G2 G3 Ginv G3_temp
          clear C_temp1 C_temp2 C_temp3 deltaX  
          toc
        end
        
        tempVN=StateVar(1:Num.Node);
        tempVA=StateVar(Num.Node+1:end);
%         seResult.Ve=reshape(tempVe,Num.Node,3);
%         seResult.Vf=reshape(tempVf,Num.Node,3);

        seResult.Ve=tempVN.*cos(tempVA);
        seResult.Vf=tempVN.*sin(tempVA);
        seResult.VOL=seResult.Ve+1i*seResult.Vf;
        seResult.VN=tempVN;
        seResult.AN=tempVA*180/pi;
        seResult.AN(seResult.VN<1e-4) = 0;
        DATA=[seResult.VN,seResult.AN];
       %% ADDED 2022
%         display(strcat('DSSE results for the time series t',num2str(iScenario),' is below:'))
% 
%         
% %         display(DATA);
%     % % %     [seResult.AN  puTrue.AN*180/pi]
%         error.VN = seResult.VN - puTrue.VN;
%         error.AN = seResult.AN - puTrue.AN*180/pi;
%     
%     
%         ERROR_VN(iSimu,:) = [max(abs(error.VN))] ;
%         ERROR_AN(iSimu,:) = [max(abs(error.AN))] ; % in degree
%         Flag_ini = Flag.ini;
%         Flag_noise = Flag.noise;
%         Flag_zeroinj = Flag.zeroinj;
%         SensorError_PN = SensorError.PN;
% %         T = table(SystemN,Flag_ini,Flag_noise, Flag_zeroinj,SensorError_PN,i, tol_temp,ERROR_VN(end,:),ERROR_AN(end,:) );
%         
%         display(strcat('completed solving DSSE for t',num2str(iScenario),' in the time series.'))
%         seResultsScenarios{iScenario}{1} = seResult;
%         seResultsScenarios{iScenario}{2} = DATA;
%         seResultsScenarios{iScenario}{3} = T;
%         seResultsScenarios{iScenario}{4} = puTrue;
%         if iScenario ==nScenario
%             display(strcat('completed solving DSSE for the whole time series.'))
%         end
%     end
    A = seResult.VN;
    A(A>1.1) = NaN;
    A(A<0.6) = NaN;
%     plot(A)
%     xlabel('Node');
%     ylabel('Voltage magnitude (p.u.)');
% figure
%         B = puTrue.VN;
%     B(B<0.01) = NaN;
%     plot(B)
%     xlabel('Node');
%     ylabel('Voltage magnitude (p.u.)');
%  figure
%     plot(tol1,'LineWidth', 1.5, 'Marker','o')
%     xlabel('Iteration');
%     ylabel('Tolerance');
    %%
% %% %% ADDED 8.9 Save Data to Excel
%         % 1 dataShare.Loads
%         % 2 dataShare.PV
%         % 3 dataShare.pfSolVN
%         % 4 dataShare.pfSolAN
%         % 5 dataShare.meaVN
%         % 6 dataShare.meaAN
%         % 7 dataShare.meaPN
%         % 8 dataShare.meaQN
%         % 9 dataShare.meaPL
%         % 10 dataShare.meaQL
%         % 11 dataShare.locPMU
%         % 12 dataShare.locPOW
%         % 13 dataShare.locSCADAVN
%         % 14 dataShare.locPseudo
%         % 15 dataShare.locVirtual
%         % 16 dataShare.locSlack   
%         
%         dataShare.Loads   = Feeder.Loads1(1:end-4,[1,4:9]);
%         dataShare.PV = PV2;
%         dataShare.pfSolVN = puTrue.VN;
%         dataShare.pfSolAN = puTrue.AN;
%         
%         dataShare.LoadsTitle   = {'Node','A-P','A-Q','B-P','B-Q','C-P','C-Q'};
%         dataShare.PVTitle = {'NodeID','PV output'};
%         dataShare.pfSolVNTitle = {'A','B','C'};
%         dataShare.pfSolANTitle = {'A','B','C'};  
% 
%         dataShare.meaVN = [idx2pha(MeaIdx.VN,Num.Node),Mea.VN(MeaIdx.VN)];
%         dataShare.meaAN = [idx2pha(MeaIdx.AN,Num.Node),Mea.AN(MeaIdx.AN)];
%         dataShare.meaPN = [idx2pha(MeaIdx.PN,Num.Node),Base.Sbase*Mea.PN(MeaIdx.PN)];
%         dataShare.meaQN = [idx2pha(MeaIdx.QN,Num.Node),Base.Sbase*Mea.QN(MeaIdx.QN)];
%         dataShare.meaPL = [idx2pha(MeaIdx.PL,Num.Line),Base.Sbase*Mea.PL(MeaIdx.PL)];
%         dataShare.meaQL = [idx2pha(MeaIdx.QL,Num.Line),Base.Sbase*Mea.QL(MeaIdx.QL)];
%        
%         dataShare.locPMU = [PMU];
%         dataShare.locPOW = [POW];
%         dataShare.locSCADAVN = [SCADA_VN];
%         dataShare.locPseudo = [idx2pha(MeaIdx.PN,Num.Node)];
%         dataShare.locVirtual = [Loc.Zeroinj];
% %         dataShare.locSlack = [SlackBus 1; SlackBus 2; SlackBus 3];
% 
%         dataShare.meaVN(:,1) = Feeder.Nodes_ID(dataShare.meaVN(:,1));
%         dataShare.meaAN(:,1) = Feeder.Nodes_ID(dataShare.meaAN(:,1));
%         dataShare.meaPN(:,1) = Feeder.Nodes_ID(dataShare.meaPN(:,1));
%         dataShare.meaQN(:,1) = Feeder.Nodes_ID(dataShare.meaQN(:,1));
%         dataShare.meaPL(:,1) = Feeder.Nodes_ID(dataShare.meaPL(:,1));
%         dataShare.meaQL(:,1) = Feeder.Nodes_ID(dataShare.meaQL(:,1));
%         dataShare.locPMU(:,1) = Feeder.Nodes_ID(dataShare.locPMU(:,1));
%         dataShare.locPOW(:,1) = Feeder.Nodes_ID(dataShare.locPOW(:,1));
%         dataShare.locSCADAVN(:,1) = Feeder.Nodes_ID(dataShare.locSCADAVN(:,1));
%         dataShare.locPseudo(:,1) = Feeder.Nodes_ID(dataShare.locPseudo(:,1));
%         dataShare.locVirtual(:,1) = Feeder.Nodes_ID(dataShare.locVirtual(:,1));
% 
%         dataShare.Ybus = num2cell(Ybus);
%         dataShare.Ybus = cellfun(@num2str , dataShare.Ybus, 'UniformOutput', false);
% 
% 
% %         dataShare.locSlack(:,1) = Feeder.Nodes_ID(dataShare.locSlack(:,1));
% 
%         dataShare.meaVNTitle = {'Node','Phase','VN'}; 
%         dataShare.meaANTitle = {'Node','Phase','AN'}; 
%         dataShare.meaPNTitle = {'Node','Phase','PN'}; 
%         dataShare.meaQNTitle = {'Node','Phase','QN'}; 
%         dataShare.meaPLTitle = {'Line','Phase','PL'}; 
%         dataShare.meaQLTitle = {'Line','Phase','QL'}; 
%         dataShare.locPMUTitle = {'Node','Phase'};
%         dataShare.locPOWTitle = {'Node','Phase'};
%         dataShare.locSCADAVNTitle = {'Node','Phase'};
%         dataShare.locPseudoTitle = {'Node','Phase'};
%         dataShare.locVirtualTitle = {'Node','Phase'};
%         dataShare.locSlackTitle = {'Node','Phase'};
    
    %%
%     xlsFileName = strcat('t',num2str(iScenario),'.xls');
%     writecell(dataShare.LoadsTitle,xlsFileName,'Sheet',1,'Range','A1')
%     writematrix(dataShare.Loads,xlsFileName,'Sheet',1,'Range','A2') 
% 
%     writecell(dataShare.PVTitle,xlsFileName,'Sheet',2,'Range','A1')
%     writematrix(dataShare.PV,xlsFileName,'Sheet',2,'Range','A2') 
% 
%     writecell(dataShare.pfSolVNTitle,xlsFileName,'Sheet',3,'Range','A1')
%     writematrix(dataShare.pfSolVN,xlsFileName,'Sheet',3,'Range','A2')  
% 
%     writecell(dataShare.pfSolANTitle,xlsFileName,'Sheet',4,'Range','A1')
%     writematrix(dataShare.pfSolAN,xlsFileName,'Sheet',4,'Range','A2')  
% 
%     writecell(dataShare.meaVNTitle,xlsFileName,'Sheet',5,'Range','A1')
%     writematrix(dataShare.meaVN,xlsFileName,'Sheet',5,'Range','A2')  
% 
%     writecell(dataShare.meaANTitle,xlsFileName,'Sheet',6,'Range','A1')
%     writematrix(dataShare.meaAN,xlsFileName,'Sheet',6,'Range','A2')  
% 
%     writecell(dataShare.meaPNTitle,xlsFileName,'Sheet',7,'Range','A1')
%     writematrix(dataShare.meaPN,xlsFileName,'Sheet',7,'Range','A2')  
% 
%     writecell(dataShare.meaQNTitle,xlsFileName,'Sheet',8,'Range','A1')
%     writematrix(dataShare.meaQN,xlsFileName,'Sheet',8,'Range','A2')  
% 
%     writecell(dataShare.meaPLTitle,xlsFileName,'Sheet',9,'Range','A1')
%     writematrix(dataShare.meaPL,xlsFileName,'Sheet',9,'Range','A2')  
% 
%     writecell(dataShare.meaQLTitle,xlsFileName,'Sheet',10,'Range','A1')
%     writematrix(dataShare.meaQL,xlsFileName,'Sheet',10,'Range','A2')  
% 
% 
% %     writecell(dataShare.locPMUTitle,xlsFileName,'Sheet',11,'Range','A1')
% %     writematrix(dataShare.locPMU,xlsFileName,'Sheet',11,'Range','A2')  
% % 
% %     writecell(dataShare.locPOWTitle,xlsFileName,'Sheet',12,'Range','A1')
% %     writematrix(dataShare.locPOW,xlsFileName,'Sheet',12,'Range','A2')  
% % 
% %     writecell(dataShare.locSCADAVNTitle,xlsFileName,'Sheet',13,'Range','A1')
% %     writematrix(dataShare.locSCADAVN,xlsFileName,'Sheet',13,'Range','A2')  
% 
% %     writecell(dataShare.locPseudoTitle,xlsFileName,'Sheet',14,'Range','A1')
% %     writematrix(dataShare.locPseudo,xlsFileName,'Sheet',14,'Range','A2')  
% 
% %     writecell(dataShare.locVirtualTitle,xlsFileName,'Sheet',15,'Range','A1')
% %     writematrix(dataShare.locVirtual,xlsFileName,'Sheet',15,'Range','A2')  
% 
% %     writecell(dataShare.locSlackTitle,xlsFileName,'Sheet',16,'Range','A1')
% %     writematrix(dataShare.locSlack,xlsFileName,'Sheet',16,'Range','A2')     
% 
% %     writematrix(dataShare.Ybus,xlsFileName,'Sheet',16,'Range','A1')
%     
%     writecell(dataShare.Ybus(:,1:100),xlsFileName,'Sheet',11,'Range','A1')  
%     writecell(dataShare.Ybus(:,101:200),xlsFileName,'Sheet',12,'Range','A1')  
%     writecell(dataShare.Ybus(:,201:300),xlsFileName,'Sheet',13,'Range','A1')  
%     writecell(dataShare.Ybus(:,301:end),xlsFileName,'Sheet',14,'Range','A1')  
% 
%     e = actxserver('Excel.Application'); 
%     ewb = e.Workbooks.Open(strcat('C:\Users\jjdon\Desktop\DSSE V2\',xlsFileName)); 
%     ewb.Worksheets.Item(1).Name = 'Loads'; 
%     ewb.Worksheets.Item(2).Name = 'PV'; 
%     ewb.Worksheets.Item(3).Name = 'pfSolVN';
%     ewb.Worksheets.Item(4).Name = 'pfSolAN';
%     ewb.Worksheets.Item(5).Name = 'meaVN';
%     ewb.Worksheets.Item(6).Name = 'meaAN';
%     ewb.Worksheets.Item(7).Name = 'meaPN';
%     ewb.Worksheets.Item(8).Name = 'meaQN';
%     ewb.Worksheets.Item(9).Name = 'meaPL';
%     ewb.Worksheets.Item(10).Name = 'meaQL';
% %     ewb.Worksheets.Item(11).Name = 'locPMU';
% %     ewb.Worksheets.Item(12).Name = 'locPOW';
% %     ewb.Worksheets.Item(13).Name = 'locSCADAVN';
% %     ewb.Worksheets.Item(14).Name = 'locPseudo';
% %     ewb.Worksheets.Item(15).Name = 'locVirtual';
% %     ewb.Worksheets.Item(16).Name = 'locSlack';
%     ewb.Worksheets.Item(11).Name = 'YbusCol100';
%     ewb.Worksheets.Item(12).Name = 'YbusCol200';
%     ewb.Worksheets.Item(13).Name = 'YbusCol300';
%     ewb.Worksheets.Item(14).Name = 'YbusCol400';
%     ewb.Save 
%     ewb.Close(false)
%     e.Quit
%% END ADDED 8.9
end

%% END ADDED 2022
% toc
end


%% ADDED 2022


%% END ADDED 2022
