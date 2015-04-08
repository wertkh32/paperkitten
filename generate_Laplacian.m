function [ mat ] = generate_Laplacian( numEle, numNode, Tt, Vol, W )

mat = zeros(numNode * 3, numNode * 3);

for i=1:numEle
    index = Tt(:,i);
    S = zeros(3,numNode);
    
     S(1,index(1)) = 1;
     S(2,index(2)) = 1;
     S(3,index(3)) = 1;
            
     S(1,index(4)) = -1;
     S(2,index(4)) = -1;
     S(3,index(4)) = -1;
            
     S = kron(S, eye(3));
    
     mat = mat + S' * S * Vol(i) * W(i);
end


end