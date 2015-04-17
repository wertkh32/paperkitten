function [edgemat Fs Rs Ss Dms] = projected_transforms( numEle, numNode, Tt, Pt, P0 )

%mat = zeros(numEle,3,3);
edgemat = zeros(9,numEle);

Fs = zeros(3,3,numEle);
Rs = zeros(3,3,numEle);
Ss = zeros(3,3,numEle);
Dms = zeros(3,3,numEle);

for i=1:numEle
    index = Tt(:,i);
    Dm = [  P0(:,index(1))- P0(:,index(4)),  P0(:,index(2))- P0(:,index(4)), P0(:,index(3))- P0(:,index(4)) ];
    Ds = [  Pt(:,index(1))- Pt(:,index(4)),  Pt(:,index(2))- Pt(:,index(4)), Pt(:,index(3))- Pt(:,index(4)) ];
    
    Bm = inv(Dm);
    
    F = Ds * Bm;
    [R S] = poldec(F);
    
    %mat(i,:) = R;
%disp(inner(S - eye(3),S - eye(3)));
    edgemat(:,i) = reshape(R * Dm, 9,1);
    
    Fs(:,:,i) = F;
    Rs(:,:,i) = R;
    Ss(:,:,i) = S;
    Dms(:,:,i) = Dm;
    
end

%disp(edgemat);

end