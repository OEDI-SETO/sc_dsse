function [zeroRow,zeroCol] = computeZero(A)
N = size(A,1);
zeroRow = [];
zeroCol = [];
for i = 1:N
    if sum(abs(A(i,:))<1e-6) == N
        zeroRow = [zeroRow;i];
    end
    if sum(abs(A(:,i))<1e-6) == N
        zeroCol = [zeroCol;i];
    end    
end
end