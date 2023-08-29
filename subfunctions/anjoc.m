function joc_an = anjoc(Num,StateVar,MeaPha)
  byphase = MeaPha.AN;
  row_v = size(byphase,1);
  joc_an = zeros(row_v, Num.StateVar*Num.Node);

  %% Jacobian for Node Voltage
  for row = 1:row_v
    i=byphase(row,1);d=byphase(row,2);
    x=2*Num.Node*(d-1)+2*i-1;y=x+1;
    
    vx = StateVar(x);
    vy = StateVar(y);
    joc_an(row,x)   =  1./(1+(vy./vx).^2).*(-vy./(vx.^2));
    joc_an(row,y)   =  1./(1+(vy./vx).^2).*1./vx;
  end