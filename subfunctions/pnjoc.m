function [joc_pn,joc_qn] = pnjoc(Num,StateVar,MeaPha,B,G)
%     byphase = MeaPha.PN;    
%     row_p = size(byphase,1);
%     joc_pn = zeros(row_p, Num.StateVar*Num.Node); joc_qn = joc_pn;
    
    %% Jacobian for Power injection
%     TermA = zeros(row_p,1); TermB = zeros(row_p,1);
%     for row = 1:row_p
%       i = byphase(row,1);d = byphase(row,2);
%       x=2*Num.Node*(d-1)+2*i-1;y=x+1;
%         for t = 1:3
%             for k = 1:Num.Node 
%                  eid = StateVar(x);
%                  fid = StateVar(y);
%                  xx=2*Num.Node*(t-1)+2*k-1;yy=xx+1;
% 
%                  ekt = StateVar(xx);
%                  fkt = StateVar(yy);
% 
%                  Bikdt = B(y/2,yy/2);
%                  Gikdt = G(y/2,yy/2);                     
% 
%                  joc_pn(y/2,xx) = (fid*Bikdt+eid*Gikdt);  %P_e=fB+eG
%                  joc_pn(y/2,yy) =   (fid*Gikdt-eid*Bikdt);  %P_f=fG-eB
%                  joc_qn(y/2,xx) =   (fid*Gikdt-eid*Bikdt);  %Q_e=fG-eB
%                  joc_qn(y/2,yy) =     (-fid*Bikdt-eid*Gikdt); %Q_f=-fB-eG
%                  TermA(row) =  TermA(row) + Gikdt*ekt - Bikdt*fkt;
%                  TermB(row) =  TermB(row) + Gikdt*fkt + Bikdt*ekt;     
%             end
%         end
%         
%         % Modify H1 because of diag elements
%         joc_pn(y/2,x)   =  joc_pn(y/2,x) + TermA(row); 
%         joc_pn(y/2,y)   =  joc_pn(y/2,y)   + TermB(row);  
%         joc_qn(y/2,x)   =  joc_qn(y/2,x)   - TermB(row);   
%         joc_qn(y/2,y)   =  joc_qn(y/2,y)      + TermA(row);   

% 
npi = size(MeaPha.PN,1);
nqi = size(MeaPha.QN,1);
nbus = Num.Node;
% fbus = Feeder.Topology(:,1); % From bus..
% tbus = Feeder.Topology(:,2); % To bus..
V = StateVar(1:nbus); % Initialize the bus voltages..
del = StateVar(nbus+1:end); % Initialize the bus angles..

% H21 - Derivative of Real Power Injections with Angles..
    H21 = zeros(npi,nbus);
    for i = 1:npi
        m = i;
        for k = 1:(nbus)
            if k == m
                for n = 1:nbus
                    H21(i,k) = H21(i,k) + V(m)* V(n)*(-G(m,n)*sin(del(m)-del(n)) + B(m,n)*cos(del(m)-del(n)));
                end
                H21(i,k) = H21(i,k) - V(m)^2*B(m,m);
            else
                H21(i,k) = V(m)* V(k)*(G(m,k)*sin(del(m)-del(k)) - B(m,k)*cos(del(m)-del(k)));
            end
        end
    end
    
    % H22 - Derivative of Real Power Injections with V..
    H22 = zeros(npi,nbus);
    for i = 1:npi
        m = i;
        for k = 1:(nbus)
            if k == m
                for n = 1:nbus
                    H22(i,k) = H22(i,k) + V(n)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
                end
                H22(i,k) = H22(i,k) + V(m)*G(m,m);
            else
                H22(i,k) = V(m)*(G(m,k)*cos(del(m)-del(k)) + B(m,k)*sin(del(m)-del(k)));
            end
        end
    end
    
    % H31 - Derivative of Reactive Power Injections with Angles..
    H31 = zeros(nqi,nbus);
    for i = 1:nqi
        m = i;
        for k = 1:(nbus)
            if k == m
                for n = 1:nbus
                    H31(i,k) = H31(i,k) + V(m)* V(n)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
                end
                H31(i,k) = H31(i,k) - V(m)^2*G(m,m);
            else
                H31(i,k) = V(m)* V(k)*(-G(m,k)*cos(del(m)-del(k)) - B(m,k)*sin(del(m)-del(k)));
            end
        end
    end
    
    % H32 - Derivative of Reactive Power Injections with V..
    H32 = zeros(nqi,nbus);
    for i = 1:nqi
        m = i;
        for k = 1:(nbus)
            if k == m
                for n = 1:nbus
                    H32(i,k) = H32(i,k) + V(n)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
                end
                H32(i,k) = H32(i,k) - V(m)*B(m,m);
            else
                H32(i,k) = V(m)*(G(m,k)*sin(del(m)-del(k)) - B(m,k)*cos(del(m)-del(k)));
            end
        end
    end

   joc_pn = [H22 H21];
   joc_qn = [H32 H31];
end  
     
     