function [en, derv] = ProjDyn_Energy(  numEle, numNode, Tt, Pt, Pold, P0, M, A, Vol, W, fixed, h, vel, fext )
%GRADIENT NEEDS VELOCITY ADDED

 q = reshape(Pt,numNode * 3,1);
 qn = reshape(Pold,numNode * 3,1);
 q0 = reshape(P0, numNode * 3,1);
 sn = qn + h * vel + h * h * M\fext;
 
 %disp(vel);
 
 [edges, F, R, S, Dm] = projected_transforms(numEle, numNode,Tt, reshape(q,3,numNode), P0);

 en = (q - sn)' * M * (q - sn) * (1.0/(2 * h * h));
 derv = (1/(h * h)) * M * ( q - sn );
 
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
    
    edgediff = L * q - edges(:,i);
    
	en = en + edgediff' * edgediff * Vol(i) * W(i)/2; 
    
    ele_derv = W(i) * Vol(i) * (L' * L) * q...
              -(L' * blkdiag(R(:,:,i),R(:,:,i),R(:,:,i)) * L) * q0 * W(i) * Vol(i);
    
    temp = zeros(numNode * 3, 1);
    
    for j=1:4

        for k=1:3
            dF = deform_grad_t_derv(k, j, Dm(:,:,i));
            dR = rotation_derv(R(:,:,i),S(:,:,i),dF);
            temp( (index(j) - 1) * 3 + k, 1 ) = -q' * L' * blkdiag(dR,dR,dR) * L * q0 * W(i) * Vol(i);
            
        end
    end
    
    ele_derv = ele_derv + temp;
    derv = derv + ele_derv;
    
 end
   
   % disp(derv);
 
end

function [dP] = piola_stress_Derv(i,j,F,R,S,Dm,miu,lambda)
    dF = deform_grad_derv(i,j,F,Dm);
    dR = rotation_derv(R,S,dF);
    dP = 2 * miu * dF + lambda * trace(R' * dF) * R + (lambda * trace(S - eye(3)) - 2 * miu) * dR;
end

function [ dF ] = deform_grad_derv(i, j, F, Dm)
    dF = -F * Dm_derv(i,j) * inv(Dm);
end

function [ dF ] = deform_grad_t_derv(i, j, Dm)
    dF = Dm_derv(i,j) * inv(Dm);
end

function [dR] = rotation_derv( R, S, dF )
    W = R' * dF;
    w = [W(2,3) - W(3,2);W(3,1) - W(1,3);W(1,2) - W(2,1)];
    r = (trace(S) * eye(3) - S) \ w;
    rx = [0 r(3) -r(2);...
          -r(3) 0 r(1);...
          r(2) -r(1) 0];
    dR = R * rx;
end

function [dV] = volume_derv(i,j,Dm)
    d = det(Dm);
    dVdDm = (1.0/6.0) * (d/abs(d)) * d * inv(Dm)';
    dV = inner(dVdDm, Dm_derv(i,j));
end

function [dDmInvT] = DmInvT_derv(i,j,Dm)
    Bm = inv(Dm);
    dDmInvT = (-Bm * Dm_derv(i,j) * Bm)';
end

function [dDm] = Dm_derv(i,j)
    mat = zeros(3,3);
    if(j < 4)
        mat(i,j) = 1;
    else
        mat(i,:) = -ones(1,3);
    end
    dDm = mat;
end