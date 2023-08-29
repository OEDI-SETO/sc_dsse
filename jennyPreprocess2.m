function output = jennyPreprocess2()
load('saved_vars.mat');
power_P1 = zeros(1,size(Y,1));
power_Q1 = zeros(1,size(Y,1));

for i = 1:size(ids_PQ,1)
%     a = strfind(string(NodeIDs),ids_PQ(i,:));
    indices = findCharacterLocation(NodeIDs, ids_PQ(i,:));
%     for j = 1:size(Y,1)
%         if a{j}==1
%             index = j;
            power_P1(indices) = power_P(i);
            power_Q1(indices) = power_Q(i);
%         end
%     end
end
power_P1(1:3) = -sum(power_P1(4:end))/3;
power_Q1(1:3) = -sum(power_Q1(4:end))/3;

SCADA_VN1 = [];
for i = 1:size(ids_voltage,1)
%     b = strfind(string(NodeIDs),ids_voltage(i,:));
    indices = findCharacterLocation(NodeIDs, ids_voltage(i,:));
%     for j = 1:size(Y,1)
%         if b{j}==1
            SCADA_VN1 = [SCADA_VN1,indices];
%         end
%     end
end



SCADA_VN = SCADA_VN1';
measure_VN = voltage';
vckt = []';

Pbig = power_P1'/base_power;
Qbig = power_Q1'/base_power;
Ybig = Y;
% zeroinj_V = Node_ZeroInjec;
zeroinj_V = [];
% vckt =[];
NodeList = ids_PQ;
Lines = [];


output = DSSE_demo7(SCADA_VN, measure_VN, NodeList, Ybig, Pbig, Qbig, vckt, Lines, zeroinj_V);

end 