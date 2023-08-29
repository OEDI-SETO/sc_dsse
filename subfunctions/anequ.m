function equ_an = anequ(Num,StateVar,MeaPha)
  byphase = MeaPha.AN;    
  row_v = size(byphase,1);
  equ_an = zeros(row_v,1);

  for row = 1: row_v
    i = byphase(row,1); d = byphase(row,2);
    vx = StateVar(2*Num.Node*(d-1)+2*i-1);
    vy = StateVar(2*Num.Node*(d-1)+2*i);
    switch d
        case 1
            equ_an(row) = atan2(vy,vx);
        case 2
            equ_an(row) = atan2(vy,vx);
        case 3
%             temp = atan(vy./vx);
%             if temp < 0
%               equ_an(row) = temp + pi;
%             else
%                equ_an(row) = temp; 
%             end
            equ_an(row) = atan2(vy,vx);
    end 
  end
