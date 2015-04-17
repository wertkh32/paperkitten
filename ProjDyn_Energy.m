function [en, derv, Hess] = ProjDyn_Energy(  numEle, numNode, Tt, Pt, Pold, P0, M, A, Vol, W, fixed, h, vel, fext, evalH )

Pt = reshape(Pt,3,numNode);

 for i=1:numNode
     if(fixed(i,1) == 1)
         Pt(:,i) = P0(:,i);
     end
 end

 q = reshape(Pt,numNode * 3,1);
 qn = reshape(Pold,numNode * 3,1);
 q0 = reshape(P0, numNode * 3,1);
 sn = qn + h * vel + h * h * M\fext;
 
 %disp(vel);
 
 [edges, F, R, S, Dm] = projected_transforms(numEle, numNode,Tt, reshape(q,3,numNode), P0);

 en = (q - sn)' * M * (q - sn) * (1.0/(2 * h * h));
 derv = (1/(h * h)) * M * ( q - sn );
 Hess = 1/(h * h) * M;
 
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
    
    H = L' * L;
    
    for j=1:4

        for k=1:3
            dF_jk = deform_grad_Ds_derv(k, j, Dm(:,:,i));
            dR_jk = rotation_derv(R(:,:,i),S(:,:,i),dF_jk);
            temp( (index(j) - 1) * 3 + k, 1 ) = -q' * L' * blkdiag(dR_jk,dR_jk,dR_jk) * L * q0 * W(i) * Vol(i);
            
            
            
            
            if(evalH)
                m = (index(j)-1) * 3 + k;

                %H(:,m) = H(:,m) - 2 * ( L' * blkdiag(dR_jk,dR_jk,dR_jk) * L * q0 );

                for a=1:4
                    for b=1:3
                         dF_ab = deform_grad_Ds_derv(b, a, Dm(:,:,i));
                         dR_ab = rotation_derv(R(:,:,i),S(:,:,i),dF_ab);
                         dS_ab = shear_scale_derv (R(:,:,i), S(:,:,i), dR_ab, dF_ab);
                         d2R = rotation_second_order_derv(R(:,:,i), dR_ab, S(:,:,i), dS_ab, dF_jk);


                         n = (index(a)-1) * 3 + b;
                         s = zeros(numNode * 3,1);
                         s1 = s;
                         s2 = s;
                         s1(n,1) = 1;
                         s2(m,1) = 1;
                         H(m,n) = H(m,n) - q' * L' * blkdiag(d2R,d2R,d2R) * L * q0  - (s1' * L' * blkdiag(dR_jk,dR_jk,dR_jk) * L * q0 ) - (s2' * L' * blkdiag(dR_ab,dR_ab,dR_ab) * L * q0 );

                    end
                end
            end
            
        end
    end
    
    H = H * W(i) * Vol(i);
    %disp(H);
    Hess = Hess + H;
    
    ele_derv = ele_derv + temp;
    derv = derv + ele_derv;
    
 end
   
  for i=1:numNode
     if(fixed(i,1) == 1)
         derv((i-1) * 3 + 1:(i-1) * 3 + 3) = 0;
         Hess((i-1) * 3 + 1:(i-1) * 3 + 3, :) = 0;
         Hess(:, (i-1) * 3 + 1:(i-1) * 3 + 3) = 0;
     end
  end

   % disp(Hess);
 
   % disp(derv);
 
end

function [dP] = piola_stress_Derv(i,j,F,R,S,Dm,miu,lambda)
    dF = deform_grad_derv(i,j,F,Dm);
    dR = rotation_derv(R,S,dF);
    dP = 2 * miu * dF + lambda * trace(R' * dF) * R + (lambda * trace(S - eye(3)) - 2 * miu) * dR;
end

function [ dF ] = deform_grad_Dm_derv(i, j, F, Dm)
    dF = -F * Dm_derv(i,j) * inv(Dm);
end

function [ dF ] = deform_grad_Ds_derv(i, j, Dm)
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

function [dS] = shear_scale_derv (R, S, dR, dF)
    dS = R' * dF - (R' * dR) * S;
end

function [d2R] = rotation_second_order_derv(R, dR_lk, S, dS_lk, dF_ij)
    W = R' * dF_ij;
    w = [W(2,3) - W(3,2);W(3,1) - W(1,3);W(1,2) - W(2,1)];
    
    M = (trace(S) * eye(3) - S);
    dM = (trace(dS_lk) * eye(3) - dS_lk);
    Minv = inv(M);
    dMinv = -Minv * dM * Minv;
    
    N = w;
    dN = dR_lk * dF_ij;
    dN = -[dN(2,3) - dN(3,2);dN(3,1) - dN(1,3);dN(1,2) - dN(2,1)];
    
    r = (trace(S) * eye(3) - S) \ w;
    Q =  [0 r(3) -r(2);...
          -r(3) 0 r(1);...
          r(2) -r(1) 0];
     
    dQ = dMinv * N + Minv * dN;
    dQ =  [0 dQ(3) -dQ(2);...
          -dQ(3) 0 dQ(1);...
          dQ(2) -dQ(1) 0];
    
    
    d2R = dR_lk * Q + R * dQ;
      
   
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

%works for Ds too hahahaha
function [dDm] = Dm_derv(i,j)
    mat = zeros(3,3);
    if(j < 4)
        mat(i,j) = 1;
    else
        mat(i,:) = -ones(1,3);
    end
    dDm = mat;
end