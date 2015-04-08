function [ A, vol ] = ProjDyn_Init( numEle, numNode, Tt, P0, M, W, FPS )
%UNTITLED Summary of this function goes here
%   numNode = n
%   mass is n-by-1

vol = zeros(numEle);

for i=1:numEle
    index = Tt(:,i);
    Dm = [  P0(:,index(1,1))- P0(:,index(4,1)),  P0(:,index(2,1))- P0(:,index(4,1)), P0(:,index(3,1))- P0(:,index(4,1)) ];
    vol(i) = (1.0/6.0) * abs(det(Dm));
end

LTL = generate_Laplacian(numEle,numNode,Tt,vol,W);

A = inv(M * (FPS * FPS) + LTL);

end

