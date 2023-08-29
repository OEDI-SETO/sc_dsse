function R = WeightMatrix(MeaIdx,SensorError)
WPN = SensorError.PN*ones(size(MeaIdx.PN));
WQN = SensorError.QN*ones(size(MeaIdx.QN));
WVN = SensorError.VN*ones(size(MeaIdx.VN));
WAN = SensorError.AN*ones(size(MeaIdx.AN));%% ADDED 2022
WPL = SensorError.PL*ones(size(MeaIdx.PL));
WQL = SensorError.QL*ones(size(MeaIdx.QL));

R =diag([WPN.^2;WQN.^2;WVN.^2;WAN.^2;WPL.^2;WQL.^2]);%% ADDED 2022
% RPN = diag(WPN.^2);RQN = diag(WQN.^2);RVN = diag(WVN.^2);
% RPL = diag(WPL.^2);RQL = diag(WQL.^2);
% RPN(row,row)= SensorError.PN^2;
% RQN(row,row)= SensorError.QN^2;
% RVN(row,row)= SensorError.VN^2;
% RPL(row,row)= SensorError.PL^2;
% RQL(row,row)= SensorError.QL^2;

    
    