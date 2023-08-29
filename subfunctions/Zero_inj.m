function [Num_ZeroInjec,Node_ZeroInjec] = Zero_inj(A)
% [minvalue,maxvalue] = Minvalue(A);
Num_ZeroInjec=0;
Node_ZeroInjec=[];
for i =1:size(A,1)
    for d= 1:size(A,2)
        if abs(A(i,d))<10^-3
            Num_ZeroInjec=Num_ZeroInjec+1;
            Node_ZeroInjec=[Node_ZeroInjec;i,d];
        end
    end
end