%%% input: index: eg,1,2,...,10,11,12,...20,21,22,...30;
%%% output: byphase: eg, a1;a2;...a10;b1;b2;...b10;c1;c2;...c10
function byphase = idx2pha(index,n)
a = index/n;
a1=find(a>0 & a<=1);
a2=find(a>1 & a<=2);
a3=find(a>2 & a<=3);
byphase = zeros(size(index,1),2);
byphase(a1,1) = index(a1);
byphase(a1,2) = 1;
byphase(a2,1) = index(a2)-n;
byphase(a2,2) = 2;
byphase(a3,1) = index(a3)-2*n;
byphase(a3,2) = 3;

%% test this function
% clear;clc;
% index = [1:30]';n=10;
% a = index/10;
% a1=find(a>0 & a<=1);
% a2=find(a>1 & a<=2);
% a3=find(a>2 & a<=3);
% byphase = zeros(size(index,1),2);
% byphase(a1,1) = index(a1);
% byphase(a1,2) = 1;
% byphase(a2,1) = index(a2)-n;
% byphase(a2,2) = 2;
% byphase(a3,1) = index(a3)-2*n;
% byphase(a3,2) = 3;































