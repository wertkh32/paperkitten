function [dR d2R] = rotation_derv_test(Tt, Pt, P0,i,j,x)
    index = Tt;
    
    P = Pt;
    P(i,index(j)) = x;
    
    Dm = [  P0(:,index(1))- P0(:,index(4)),  P0(:,index(2))- P0(:,index(4)), P0(:,index(3))- P0(:,index(4)) ];
    Ds = [  P(:,index(1))- P(:,index(4)),  P(:,index(2))- P(:,index(4)), P(:,index(3))- P(:,index(4)) ];
    
    Bm = inv(Dm);
    
    F = Ds * Bm;
    [R S] = poldec(F);

   
    
    dF = deform_grad_Ds_derv(i, j, Dm);
    dR = rotation_derv( R, S, dF );
    dS = shear_scale_derv (R, S, dR, dF);
    
    [d2R Q dQ] = rotation_second_order_derv(R, dR, S, dS, dF);
    
     R = reshape(R,9,1);
     dR = reshape(dR,9,1);
     d2R = reshape(d2R,9,1);
     %Q = reshape(Q,9,1);
     %dQ = reshape(dQ,9,1);
    %disp(reshape(dR,9,1));

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

function [d2R Q dQ] = rotation_second_order_derv(R, dR_lk, S, dS_lk, dF_ij)
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