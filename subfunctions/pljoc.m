function [joc_pl,joc_ql] = pljoc(Num,StateVar,MeaPha,B,G,Loc)
    byphase = MeaPha.PL;    
    row_p   = size(byphase,1);
    joc_pl  = zeros(row_p, Num.StateVar*Num.Node); joc_ql = joc_pl;   
    TermE=zeros(row_p,1); TermF=zeros(row_p,1);
            %% Jacobian for Line Power flow(P and Q)
    for row = 1:row_p
      j = byphase(row,1);d = byphase(row,2);
        for t=1:3
            i=Loc.Line(j,1);k=Loc.Line(j,2);
            Gikdt = G(Num.Node*(d-1)+i,Num.Node*(t-1)+k);
            Bikdt = B(Num.Node*(d-1)+i,Num.Node*(t-1)+k);
            ekd = StateVar(2*Num.Node*(d-1)+2*k-1);
            fkd = StateVar(2*Num.Node*(d-1)+2*k); 
            joc_pl(row,2*Num.Node*(t-1)+2*i-1)=ekd*Gikdt+fkd*Bikdt; 
            joc_pl(row,2*Num.Node*(t-1)+2*i)=-ekd*Bikdt+fkd*Gikdt; 
            joc_pl(row,2*Num.Node*(t-1)+2*k-1)=-ekd*Gikdt-fkd*Bikdt;
            joc_pl(row,2*Num.Node*(t-1)+2*k)=ekd*Bikdt-fkd*Gikdt;     

            joc_ql(row,2*Num.Node*(t-1)+2*i-1)=fkd*Gikdt-ekd*Bikdt; 
            joc_ql(row,2*Num.Node*(t-1)+2*i)=-fkd*Bikdt-ekd*Gikdt; 
            joc_ql(row,2*Num.Node*(t-1)+2*k-1)=ekd*Bikdt-fkd*Gikdt;
            joc_ql(row,2*Num.Node*(t-1)+2*k)=ekd*Gikdt+fkd*Bikdt;                        
        end
        %modify the diagonal elemens
        t=d; i=Loc.Line(j,1); k=Loc.Line(j,2);
        Gikdt = G(Num.Node*(d-1)+i,Num.Node*(t-1)+k);
        Bikdt = B(Num.Node*(d-1)+i,Num.Node*(t-1)+k);
        eit = StateVar(2*Num.Node*(t-1)+2*i-1);
        fit = StateVar(2*Num.Node*(t-1)+2*i);  
        ekt = StateVar(2*Num.Node*(t-1)+2*k-1);
        fkt = StateVar(2*Num.Node*(t-1)+2*k);                     
        TermE(row)=(Gikdt*eit-Bikdt*fit)-(Gikdt*ekt-Bikdt*fkt);
        TermF(row)=(Bikdt*eit+Gikdt*fit)-(Bikdt*ekt+Gikdt*fkt);

        joc_pl(row,2*Num.Node*(d-1)+2*k-1)=...
          joc_pl(row,2*Num.Node*(d-1)+2*k-1)+TermE(row);
        joc_pl(row,2*Num.Node*(d-1)+2*k)=...
          joc_pl(row,2*Num.Node*(d-1)+2*k)+TermF(row);
        joc_ql(row,2*Num.Node*(d-1)+2*k-1)=...
          joc_ql(row,2*Num.Node*(d-1)+2*k-1)-TermF(row);
        joc_ql(row,2*Num.Node*(d-1)+2*k)=...
          joc_ql(row,2*Num.Node*(d-1)+2*k)+TermE(row);    
    end    
%  row <-- Num.Line*(d-1)+j
    joc_pl = - joc_pl;
    joc_ql = - joc_ql;
       