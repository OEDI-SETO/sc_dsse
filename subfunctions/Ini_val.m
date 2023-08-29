function [InitialValue,StateVar_true]=...
    Ini_val(Flag,Num,puTrue,SlackBus)

Ve1 = reshape(puTrue.Ve,Num.Node,1);
Vf1 = reshape(puTrue.Vf,Num.Node,1);
VV1 = [Ve1,Vf1];
StateVar_true= reshape(VV1',Num.StateVar * Num.Node,1);%a:e1f1e2f2...;b:e1f1e2f2...;...

switch Flag.ini
  
    case 1  
        
Ve_temp = zeros(Num.Node,1);
Vf_temp = zeros(Num.Node,1);

for i=1:Num.Node
    for d=1:1
        if puTrue.Ve(i,d) == 0
            Ve_temp(i,d)=0;Vf_temp(i,d)=0;
        else
            if d==1
                Ve_temp(i,d)=1;Vf_temp(i,d)=0;
            end
            if d==2
                Ve_temp(i,d)=-0.5;Vf_temp(i,d)=-sqrt(3)/2;
            end
            if d==3
                Ve_temp(i,d)=-0.5;Vf_temp(i,d)=sqrt(3)/2;
            end
        end
    end
end

% Ve1p = reshape(Ve_temp,3*Num.Node,1);
% Vf1p = reshape(Vf_temp,3*Num.Node,1);
Ve1p = Ve_temp;
Vf1p = Vf_temp;
VV1p = [Ve1p;Vf1p];
InitialValue = reshape(VV1p',Num.StateVar*Num.Node,1);%a:e1f1e2f2...;b:e1f1e2f2...;..

% InitialValue(2*SlackBus:2*Num.Node:4*Num.Node+2*SlackBus)=...
%     StateVar_true(2*SlackBus:2*Num.Node:4*Num.Node+2*SlackBus);
% InitialValue(2*SlackBus-1:2*Num.Node:4*Num.Node+2*SlackBus-1)=...
%     StateVar_true(2*SlackBus-1:2*Num.Node:4*Num.Node+2*SlackBus-1);


%%
    case 2
        Ve = StateVar_true(1:2:end);
        Vf = StateVar_true(2:2:end);
        VV = Ve +1i*Vf;
        VM = abs(VV);
        VA = angle(VV);
        InitialValue = [VM;VA];
        

    case 3
Ve_temp = zeros(Num.Node,3);
Vf_temp = zeros(Num.Node,3);

for i=1:Num.Node
    for d=1:3
        if puTrue.Ve(i,d) == 0
            Ve_temp(i,d)= 0 ;Vf_temp(i,d)= 0 ;
        else
            if d==1
%                 Ve_temp(i,d)=puTrue.Ve(SlackBus,1);Vf_temp(i,d)=puTrue.Vf(SlackBus,1);
                Ve_temp(i,d)=puTrue.VN(SlackBus,1)*cos(puTrue.AN(SlackBus,1));
                Vf_temp(i,d)=puTrue.VN(SlackBus,1)*sin(puTrue.AN(SlackBus,1));
            end
            if d==2
%                 Ve_temp(i,d)=puTrue.Ve(SlackBus,2);Vf_temp(i,d)=puTrue.Vf(SlackBus,2);
                Ve_temp(i,d)=puTrue.VN(SlackBus,2)*cos(puTrue.AN(SlackBus,2));
                Vf_temp(i,d)=puTrue.VN(SlackBus,2)*sin(puTrue.AN(SlackBus,2));                
            end
            if d==3
%                 Ve_temp(i,d)=puTrue.Ve(SlackBus,3);Vf_temp(i,d)=puTrue.Vf(SlackBus,3);
                Ve_temp(i,d)=puTrue.VN(SlackBus,3)*cos(puTrue.AN(SlackBus,3));
                Vf_temp(i,d)=puTrue.VN(SlackBus,3)*sin(puTrue.AN(SlackBus,3));                
            end
        end
    end
end

Ve1p = reshape(Ve_temp,3*Num.Node,1);
Vf1p = reshape(Vf_temp,3*Num.Node,1);
VV1p = [Ve1p,Vf1p];
InitialValue = reshape(VV1p',Num.StateVar*Num.Node,1);%a:e1f1e2f2...;b:e1f1e2f2...;..

InitialValue(2*SlackBus:2*Num.Node:4*Num.Node+2*SlackBus)=...
    StateVar_true(2*SlackBus:2*Num.Node:4*Num.Node+2*SlackBus);
InitialValue(2*SlackBus-1:2*Num.Node:4*Num.Node+2*SlackBus-1)=...
    StateVar_true(2*SlackBus-1:2*Num.Node:4*Num.Node+2*SlackBus-1);

        
end



