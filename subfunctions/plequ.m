function [equ_pl,equ_ql] = plequ(Num,StateVar,MeaPha,B,G,Loc)
  byphase = MeaPha.PL;    
  row_p   = size(byphase,1);
  equ_pl = zeros(row_p,1); equ_ql = equ_pl;
  TermC  = zeros(row_p,1); TermD  = zeros(row_p,1);
  %%
for row = 1: row_p
  j= byphase(row,1); d = byphase(row,2);
  i=Loc.Line(j,1);k=Loc.Line(j,2);
  for t=1:3
    Gikdt = G(Num.Node*(d-1)+i,Num.Node*(t-1)+k);
    Bikdt = B(Num.Node*(d-1)+i,Num.Node*(t-1)+k);
    eit = StateVar(2*Num.Node*(t-1)+2*i-1);
    fit = StateVar(2*Num.Node*(t-1)+2*i);  
    ekt = StateVar(2*Num.Node*(t-1)+2*k-1);
    fkt = StateVar(2*Num.Node*(t-1)+2*k);                      
    TermC(row)=TermC(row)+(Gikdt*eit-Bikdt*fit)-(Gikdt*ekt-Bikdt*fkt);
    TermD(row)=TermD(row)+(Bikdt*eit+Gikdt*fit)-(Bikdt*ekt+Gikdt*fkt);
  end
  
  ekd = StateVar(2*Num.Node*(d-1)+2*k-1);
  fkd = StateVar(2*Num.Node*(d-1)+2*k);
  Vkd=ekd+1i*fkd;
  Ijd=TermC(row)+1i*TermD(row);
  Sflow=Vkd.*conj(Ijd);
  equ_pl(row)=-real(Sflow);
  equ_ql(row)=-imag(Sflow);
end
 

% for row = 1: row_p
%   j= byphase(row,1); d = byphase(row,2); k=Loc.Line(j,2);
%   ekd = StateVar(2*Num.Node*(d-1)+2*k-1);
%   fkd = StateVar(2*Num.Node*(d-1)+2*k);
%   Vkd=ekd+1i*fkd;
%   Ijd=TermC(row)+1i*TermD(row);
%   Sflow=Vkd.*conj(Ijd);
%   equ_pl(row)=-real(Sflow);
%   equ_ql(row)=-imag(Sflow);
% end