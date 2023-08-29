function [VN, AN] = processPowData

load("powData.mat")

VN(1) = (max(powData60A)-min(powData60A)) /2;
VN(2) = (max(powData60B)-min(powData60B)) /2;
VN(3) = (max(powData60C)-min(powData60C)) /2;
AN(1) = acos(powData60A(1)/VN(1));
AN(2) = -acos(powData60B(1)/VN(2));
AN(3) = acos(powData60C(1)/VN(3));

% AN(1) = acos(powData60A(2)/sqrt(2)/VN(1))-2*pi*60*1/1000;
return
