function [ ret retvel ] = ProjDyn_timestep( numEle, numNode, Tt, Pt, P0, M, A, Vol, W, fixed, h, vel, fext )
    
    qn = reshape(Pt,numNode * 3,1);
    q = qn;
    sn = qn + h * vel + h * h * M\fext;
    
    for n=1:10
        b = M * (1.0 / (h*h)) * sn;
        
        edges = projected_transforms(numEle, numNode,Tt, reshape(q,3,numNode), P0);

        for i=1:numEle
            index = Tt(:,i);
            L = zeros(3,numNode);
            
            L(1,index(1)) = 1;
            L(2,index(2)) = 1;
            L(3,index(3)) = 1;
            
            L(1,index(4)) = -1;
            L(2,index(4)) = -1;
            L(3,index(4)) = -1;
            
            L = kron(L, eye(3));
            
            b = b + L' * edges(:,i) * Vol(i) * W(i); 
        end
        
        q = A * b; 
        q = reshape(q,3,numNode);
        
        for j=1:numNode
            if(fixed(j,1) == 1)
                q(:,j) = P0(:,j);
            end
        end
        
         q = reshape(q,numNode * 3,1);
    end
   
    ret = q;
    retvel = (q - qn)/h;
    
end

