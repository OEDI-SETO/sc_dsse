function Noise = Add_noise(Flag,Num,SensorError)
switch Flag.noise

    
    case 1
%% Gauss noise, three phase the same
spn = SensorError.PN;
sqn = SensorError.QN;
svn = SensorError.VN;
san = SensorError.AN;%% ADDED 2022
spl = SensorError.PL;
sql = SensorError.QL;


for i=1:Num.Node
    sspn(i,1)=TruncatedGaussian(spn, [- 3 * spn 3 * spn]);
    ssqn(i,1)=TruncatedGaussian(sqn, [- 3 * sqn 3 * sqn]);
    ssvn(i,1)=TruncatedGaussian(svn, [- 3 * svn 3 * svn]);
    ssan(i,1)=TruncatedGaussian(san, [- 3 * san 3 * san]);%% ADDED 2022
end

for i=1:Num.Line
    sspl(i,1)=TruncatedGaussian(spl, [- 3 * spl 3 * spl]);
    ssql(i,1)=TruncatedGaussian(sql, [- 3 * sql 3 * sql]);
end


Noise.PN = [1-sspn , 1-sspn , 1-sspn];
Noise.QN = [1-ssqn , 1-ssqn , 1-ssqn];
Noise.VN = [1-ssvn , 1-ssvn , 1-ssvn];
Noise.AN = [1-ssan , 1-ssan , 1-ssan];%% ADDED 2022
Noise.PL = [1-sspl , 1-sspl , 1-sspl];
Noise.QL = [1-ssql , 1-ssql , 1-ssql];

    case 11
%% Gauss noise, three phase different
spn = SensorError.PN;
sqn = SensorError.QN;
svn = SensorError.VN;
san = SensorError.AN;%% ADDED 2022
spl = SensorError.PL;
sql = SensorError.QL;


for i=1:Num.Node
    sspn1(i,1)=TruncatedGaussian(spn, [- 3 * spn 3 * spn]);
    ssqn1(i,1)=TruncatedGaussian(sqn, [- 3 * sqn 3 * sqn]);
    ssvn1(i,1)=TruncatedGaussian(svn, [- 3 * svn 3 * svn]);
    ssan1(i,1)=TruncatedGaussian(san, [- 3 * san 3 * san]);%% ADDED 2022
    
    sspn2(i,1)=TruncatedGaussian(spn, [- 3 * spn 3 * spn]);
    ssqn2(i,1)=TruncatedGaussian(sqn, [- 3 * sqn 3 * sqn]);
    ssvn2(i,1)=TruncatedGaussian(svn, [- 3 * svn 3 * svn]);
    ssan2(i,1)=TruncatedGaussian(san, [- 3 * san 3 * san]);%% ADDED 2022
    
    sspn3(i,1)=TruncatedGaussian(spn, [- 3 * spn 3 * spn]);
    ssqn3(i,1)=TruncatedGaussian(sqn, [- 3 * sqn 3 * sqn]);
    ssvn3(i,1)=TruncatedGaussian(svn, [- 3 * svn 3 * svn]);
    ssan3(i,1)=TruncatedGaussian(san, [- 3 * san 3 * san]);%% ADDED 2022  
end

for i=1:Num.Line
    sspl1(i,1)=TruncatedGaussian(spl, [- 3 * spl 3 * spl]);
    ssql1(i,1)=TruncatedGaussian(sql, [- 3 * sql 3 * sql]);
    
    sspl2(i,1)=TruncatedGaussian(spl, [- 3 * spl 3 * spl]);
    ssql2(i,1)=TruncatedGaussian(sql, [- 3 * sql 3 * sql]); 
    
    sspl3(i,1)=TruncatedGaussian(spl, [- 3 * spl 3 * spl]);
    ssql3(i,1)=TruncatedGaussian(sql, [- 3 * sql 3 * sql]);    
end


Noise.PN = [1-sspn1 , 1-sspn2 , 1-sspn3];
Noise.QN = [1-ssqn1 , 1-ssqn2 , 1-ssqn3];
Noise.VN = [1-ssvn1 , 1-ssvn2 , 1-ssvn3];
Noise.AN = [1-ssan1 , 1-ssan2 , 1-ssan3];%% ADDED 2022
Noise.PL = [1-sspl1 , 1-sspl2 , 1-sspl3];
Noise.QL = [1-ssql1 , 1-ssql2 , 1-ssql3];
    
    
    case 2
%% No noise        
Noise.PN = 1;
Noise.QN = 1;
Noise.VN = 1;
Noise.AN = 1;%% ADDED 2022
Noise.PL = 1;
Noise.QL = 1;



end

end



function y = f(x,n)
    for jj = 1:n
        if rand(1) > 0.5 
            y(jj,1:3) = x; 
        else
            y(jj,1:3) = -x; 
        end
    end
end