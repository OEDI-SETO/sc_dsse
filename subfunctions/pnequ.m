function [equ_pn,equ_qn] = pnequ(Num,StateVar,MeaPha,B,G)
%   byphase = MeaPha.PN;    
%   row_p = size(byphase,1);
%   equ_pn = zeros(row_p,1); equ_qn = equ_pn;
%   %% Power Injection Measure Equ
%   for row = 1:row_p
%     i = byphase(row,1) ; d = byphase(row,2);  
%     x=2*Num.Node*(d-1)+2*i-1;y=x+1;
%     for k = 1:Num.Node
%       for t = 1:3
%         xx=2*Num.Node*(t-1)+2*k-1;yy=xx+1;
%         eid = StateVar(x);
%         fid = StateVar(y);
%         ekt = StateVar(xx);
%         fkt = StateVar(yy);                     
%         Bikdt = B(y/2,yy/2);
%         Gikdt = G(y/2,yy/2);
%         Term1 = fkt*Gikdt+ekt*Bikdt;       %fG+eB 
%         Term2 = ekt*Gikdt-fkt*Bikdt;       %eG-fB
%         equ_pn(row) = equ_pn(row) + fid*Term1 + eid*Term2;
%         equ_qn(row) = equ_qn(row) + fid*Term2 - eid*Term1;
%       end
%     end    
%   end
npi = size(MeaPha.PN,1);
nqi = size(MeaPha.QN,1);
nbus = Num.Node;
V = StateVar(1:nbus); % Initialize the bus voltages..
del = StateVar(nbus+1:end); % Initialize the bus angles..

ppi = MeaPha.PN(:,1);
qi  = MeaPha.QN(:,1);
    h2 = zeros(npi,1);
    h3 = zeros(nqi,1);
% 
    for i = 1:npi
        m = ppi(i);
        for k = 1:nbus
            h2(i) = h2(i) + V(m)*V(k)*(G(m,k)*cos(del(m)-del(k)) + B(m,k)*sin(del(m)-del(k)));
        end
    end
    
    for i = 1:nqi
        m = qi(i);
        for k = 1:nbus
            h3(i) = h3(i) + V(m)*V(k)*(G(m,k)*sin(del(m)-del(k)) - B(m,k)*cos(del(m)-del(k)));
        end
    end

b=-B; %线路电导矩阵
g=-G; %线路电纳矩阵
P=zeros(3,1); %初始化，节点注入功率
Q=zeros(3,1);
PP=zeros(3,3); %线路注入功率
QQ=PP;

for i=1:npi
    P_P=0;
    Q_Q=0;
    for j=1:npi
        if(j~=i)
        P_P=P_P+V(i)*V(j)*(G(i,j)*cos(del(i)-del(j))+B(i,j)*sin(del(i)-del(j)));
        Q_Q=Q_Q+V(i)*V(j)*(G(i,j)*sin(del(i)-del(j))-B(i,j)*cos(del(i)-del(j)));
        PP(i,j)=(V(i)^2)*g(i,j)-V(i)*V(j)*(g(i,j)*cos(del(i)-del(j))...
+b(i,j)*sin(del(i)-del(j)));
        QQ(i,j)=-(V(i)^2)*b(i,j)-V(i)*V(j)*(g(i,j)*sin(del(i)-del(j))...
-b(i,j)*cos(del(i)-del(j)));
        end
    end
    P(i)=(V(i)^2)*G(i,i)+P_P;
    Q(i)=-(V(i)^2)*B(i,i)+Q_Q;
end

    equ_pn = h2;
    equ_qn  = h3;