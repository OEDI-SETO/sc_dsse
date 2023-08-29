%%% input: byphase: eg, a1;a2;...a10;b1;b2;...b10;c1;c2;...c10;
%%% output: index: eg,1,2,...,10,11,12,...20,21,22,...30;
function index = pha2idx(byphase,n)
if isempty(byphase)
    index = [];
else
    node  = byphase(:,1);
    phase = byphase(:,2);
    row   = size(byphase,1);
    index = zeros(row,1);
    for i = 1:row
    switch phase(i)
      case 1
        index(i) = node(i);
      case 2
        index(i) = node(i) + n;
      case 3
        index(i) = node(i) + 2*n;
    end
end
end