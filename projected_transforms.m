function [edgemat] = projected_transforms( numEle, numNode, Tt, Pt, P0 )

%mat = zeros(numEle,3,3);
edgemat = zeros(9:numEle);

for i=1:numEle
    index = Tt(:,i);
    Dm = [  P0(:,index(1,1))- P0(:,index(4,1)),  P0(:,index(2,1))- P0(:,index(4,1)), P0(:,index(3,1))- P0(:,index(4,1)) ];
    Ds = [  Pt(:,index(1,1))- Pt(:,index(4,1)),  Pt(:,index(2,1))- Pt(:,index(4,1)), Pt(:,index(3,1))- Pt(:,index(4,1)) ];
    
    Bm = inv(Dm);
    
    F = Ds * Bm;
    [R S] = poldec(F);
    
    %mat(i,:) = R;
    edgemat(:,i) = reshape(R * Dm, 9,1);
end

%disp(edgemat);

end