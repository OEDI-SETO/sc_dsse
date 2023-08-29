function joc_vn = vnjoc(Num,StateVar,MeaPha)
%   byphase = MeaPha.VN;
%   row_v = size(byphase,1);
%   joc_vn = zeros(row_v, Num.StateVar*Num.Node);
% 
%   %% Jacobian for Node Voltage
%   for row = 1:row_v
%     i=byphase(row,1);d=byphase(row,2);
%     x=2*Num.Node*(d-1)+2*i-1;y=x+1;
%     joc_vn(row,x) =  2 * StateVar(x);
%     joc_vn(row,y)   =  2 * StateVar(y);
% 
%     
%   end

nvi = size(MeaPha.VN,1);
nbus = Num.Node;

    % Jacobian..
    % H11 - Derivative of V with respect to angles.. All Zeros
    H11 = zeros(nvi,nbus);

    % H12 - Derivative of V with respect to V.. 
    H12 = zeros(nvi,nbus);
    for k = 1:nvi
        for n = 1:nbus
            if n == k
                H12(k,n) = 1;
            end
        end
    end

    joc_vn = [H12 H11];
