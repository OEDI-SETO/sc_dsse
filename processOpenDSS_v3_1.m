function processOpenDSS_v3_1(Pbig, Qbig, Ybig, Lines, NodeList, zeroinj_V) 
%%%%%% Bus information
% vckt = DSSCircuit.YNodeVarray;
vckt = ones(size(Ybig,1),1);
vckt_len = length(vckt);
Vbig = vckt;
% Vbig = (vckt(1:2:vckt_len)-1i*vckt(2:2:vckt_len))';
% NodeList = DSSCircuit.YNodeOrder;R
% Vbig(1:256,1) = Vbig1(4:259);
% Vbig(257:259,1) = Vbig1(1:3);
% Vbig(260:269,1) = Vbig1(260:end);

%% get bus information
% ## Buses = getbusinfo();

% DSSText.Command = 'vsource.source.enabled=no';
% DSSText.command = 'batchedit load..* enabled=no';
% DSSText.Command = 'solve';
% 
% Ybus = DSSCircuit.SystemY;
% NodeList1 = DSSCircuit.YNodeOrder;
% nb_node = length(NodeList1);
% Y = reshape(Ybus,nb_node*2,nb_node)';
% % create complex polyphase nodal Y
% Ybig = Y(:,1:2:(nb_node*2)) + 1i*Y(:,2:2:nb_node*2);
% Y_dss = Ybig;

Ibig=Ybig*Vbig;
% Pbig=real(Vbig.*conj(Ibig));
% Qbig=imag(Vbig.*conj(Ibig));

Ybus = Ybig;
V_node = Vbig;
I_node = Ibig;
P_node = Pbig + 1i*Qbig;

for i=1:size(V_node,1)
    for d=1:1
            if d==1
                Ve_temp(i,d)=1;Vf_temp(i,d)=0;
                V_node(i,d) = Ve_temp(i,d) + j*Vf_temp(i,d);
            end
            if d==2
                Ve_temp(i,d)=-0.5;Vf_temp(i,d)=-sqrt(3)/2;
                V_node(i,d) = Ve_temp(i,d) + j*Vf_temp(i,d);
            end
            if d==3
                Ve_temp(i,d)=-0.5;Vf_temp(i,d)=sqrt(3)/2;
                V_node(i,d) = Ve_temp(i,d) + j*Vf_temp(i,d);
            end
    end
end
% if isempty(zeroinj_V) == 0
%     row_Modify = zeroinj_V(:,1);
%     col_Modify = zeroinj_V(:,2);
%     num_Modify = size(zeroinj_V,1);
%     for i = 1:num_Modify
%         V_node(row_Modify(i),col_Modify(i)) = 0;
%     end
% end

% V_node = 7415.5* V_node; % large system
V_node = 2.401694479938087e+03* V_node; % small system


Feeder.NumN = size(Ybig,1);
Feeder.NumL = size(Ybig,1);
% Feeder.NumL = size(Lines,1);

P_line = zeros(Feeder.NumL,1);
I_line = zeros(Feeder.NumL,1);

Ybus = sparse(Ybus);

FileName = 'IEEE';


save(strcat(FileName,'_temp.mat'),...
    'Feeder', 'Ybus', 'P_node', 'V_node',...
    'I_node', 'I_line', 'P_line')
return

%% 
a = NodeList;
b1 = P_node;
b2 = I_node;
b3 = V_node;
c  = Ybus;

[a_new,index,numDim] = calcIndex(a);
[b1_new] = funConvertArrayToThreePhaseMatrix(b1,index,numDim);
[b2_new] = funConvertArrayToThreePhaseMatrix(b2,index,numDim);
[b3_new] = funConvertArrayToThreePhaseMatrix(b3,index,numDim);
[c_new] = funConvertYbusToThreePhaseMatrix(c,index,numDim);

% % [a_new,b1_new] = funConvertArrayToThreePhaseMatrix(a, b1);
% % [~,b2_new] = funConvertArrayToThreePhaseMatrix(a, b2);
% % [~,b3_new] = funConvertArrayToThreePhaseMatrix(a, b3);
% % [~,c_new] = funConvertYbusToThreePhaseMatrix(a, c);
% % return
% % 
P_node = b1_new;
I_node = b2_new;
V_node = b3_new;

for i=1:size(V_node,1)
    for d=1:3
            if d==1
                Ve_temp(i,d)=1;Vf_temp(i,d)=0;
                V_node(i,d) = Ve_temp(i,d) + j*Vf_temp(i,d);
            end
            if d==2
                Ve_temp(i,d)=-0.5;Vf_temp(i,d)=-sqrt(3)/2;
                V_node(i,d) = Ve_temp(i,d) + j*Vf_temp(i,d);
            end
            if d==3
                Ve_temp(i,d)=-0.5;Vf_temp(i,d)=sqrt(3)/2;
                V_node(i,d) = Ve_temp(i,d) + j*Vf_temp(i,d);
            end
    end
end
if isempty(zeroinj_V) == 0
    row_Modify = zeroinj_V(:,1);
    col_Modify = zeroinj_V(:,2);
    num_Modify = size(zeroinj_V,1);
    for i = 1:num_Modify
        V_node(row_Modify(i),col_Modify(i)) = 0;
    end
end

% V_node = 7415.5* V_node; % large system
V_node = 2.401694479938087e+03* V_node; % small system

% for i =1:6
%     subplot(2,3,i)
%     compass(V_node(i,:))
% end

Ybus = c_new;
% % % % % load('Ybus.mat');
NodeList = a_new;
% % % 
% % % % P_node(132,:) = [];
% % % % I_node(132,:) = [];
% % % % V_node(132,:) = [];
% % % % tempIdx = [132,132+132,132+132*2];
% % % % Ybus(tempIdx,:) = [];
% % % % Ybus(:,tempIdx) = [];
% % % %%
% % % % x = double(string(NodeList));
% % % % [xSorted,idxSorted] = sort(x);
% % % % 
% % % % P_nodeNew = P_node(idxSorted,:);
% % % % I_nodeNew = I_node(idxSorted,:);
% % % % V_nodeNew = V_node(idxSorted,:);
% % % 
% % % 
% % % %%
% % % % V_node(end,:) = V_node(end,:)*9;
% % % % I_node(end,:) = I_node(end,:)*9;
% % % % P_node(end,:) = P_node(end,:)*81;
%%%%% Feeder information %%%%%%

% load('matlab.mat')
Feeder.NumN = size(NodeList,1);
Feeder.NumL = size(Lines,1);

%%%%% Line information %%%%%%

% P_temp = V_temp.*conj(Res.Ilines);
newTopology = zeros(Feeder.NumL,2);
for iRow = 1:Feeder.NumL
%     iRow = 8;
    oldTopology1 = Lines(iRow).bus1;
    oldTopology2 = Lines(iRow).bus2;
%     temp1 = split(oldTopology1,'.');
%     temp2 = split(oldTopology2,'.');

    lenOfList = size(oldTopology1,1);
    temp1 = [];
    for itemp = 1:lenOfList
        temp1 = [temp1;strsplit(oldTopology1(itemp,:),'.')];
    end
    
    lenOfList = size(oldTopology2,1);
    temp2 = [];
    for itemp = 1:lenOfList
        temp2 = [temp2;strsplit(oldTopology2(itemp,:),'.')];
    end
    
    idxFromNode = [];
    idxToNode   = [];
    for jjj = 1:length(a_new)
        
        if strcmpi(temp1{1},a_new(jjj)) % string compare, case insensitive, i.e., "yes"=="Yes
            idxFromNode = jjj;
        end
        if strcmpi(temp2{1},a_new(jjj))
            idxToNode = jjj;
        end
    end
    if isempty(idxFromNode) || isempty(idxToNode)
        newTopology(iRow,1:2) = NaN;
    else   
        newTopology(iRow,1) = idxFromNode;
        newTopology(iRow,2) = idxToNode;
    end
end
idxNaN = find(isnan(newTopology(:,1)));
newTopology(idxNaN,:) = [];
Feeder.Topology(:,1:2) = newTopology;




P_line = zeros(Feeder.NumL,3);
I_line = zeros(Feeder.NumL,3);

% for k = 1:Feeder.NumL
%       N1 = Feeder.Topology(k,1);
%       N2 = Feeder.Topology(k,2);
%       P_line(k,:) = P_temp(N2,:);
%       I_line(k,:) = Res.Ilines(N2,:);
% end

Ybus = sparse(Ybus);

FileName = 'IEEE';


save(strcat(FileName,'_temp.mat'),...
    'Feeder', 'Ybus', 'P_node', 'V_node',...
    'I_node', 'I_line', 'P_line')
end

%%
% module 1: calculate the index of the 269 elements in the 384 elements
% module 2: define a zero vector 384*1, assign the 269 elements into the
% 384 locations using the index in module 1 (for Pnode, Inode, Vnode)
% module 3: define a zero matrix 384*384, assign the 269*269 elements into
% the 384*384 locations using the index in module 1 (for Ybus)

%


%%





