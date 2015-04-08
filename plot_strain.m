function [ h ] = plot_strain( numEle,Tt,P0, P )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

Bm = zeros(3,3, numEle);
for i=1:numEle
    index = Tt(:,i);
    tmp = [  P0(:,index(1,1))- P0(:,index(4,1)),  P0(:,index(2,1))- P0(:,index(4,1)), P0(:,index(3,1))- P0(:,index(4,1)) ];
    Bm(:,:,i) = inv(tmp);
end

E = 1e6;
nu = 0.45;
M = [ 1-nu nu nu 0 0 0;...
    nu 1-nu nu 0 0 0;... 
    nu nu 1-nu 0 0 0;...  
    0 0 0 1-2*nu 0 0;... 
    0 0 0 0 1-2*nu 0;... 
    0 0 0 0 0 1-2*nu];

%% deformaton gradient, strain, stress
Ds = zeros(3,3, numEle);
S_Von_Mises = zeros(numEle, 1);
for i=1:numEle
    index = Tt(:,i);
    Ds(:,:, i) = [  P(:,index(1,1))- P(:,index(4,1)),  P(:,index(2,1))- P(:,index(4,1)), P(:,index(3,1))- P(:,index(4,1)) ];
    
    F = Ds(:,:,i)*Bm(:,:,i);
    Strain = 0.5 * (F' * F - eye(3));
    %disp(Strain);
    S = M * [Strain(1,1) Strain(2,2) Strain(3,3) Strain(1,2) Strain(2,3) Strain(3,1)]';
    S_Von_Mises(i,1) = sqrt(0.5 * ( (S(1,1)-S(2,1))^2 + (S(2,1)-S(3,1))^2 + (S(3,1)-S(1,1))^2 + 6*(S(4,1)^2+S(5,1)^2+S(6,1)^2) ));
end


S_Von_Mises_face = zeros(4*numEle,1);
for i=1:size(Tt,2)
    TSum((i-1)*4+1, :) = [Tt(1,i) Tt(2,i) Tt(3,i)];
    TSum((i-1)*4+2, :) = [Tt(2,i) Tt(3,i) Tt(4,i)];
    TSum((i-1)*4+3, :) = [Tt(3,i) Tt(4,i) Tt(1,i)];
    TSum((i-1)*4+4, :) = [Tt(4,i) Tt(1,i) Tt(2,i)];
    S_Von_Mises_face((i-1)*4+1:i*4, 1) = S_Von_Mises(i);
end
figure();
view(21, 0);
caxis([0, 0.02]);
h = patch('Faces', TSum, 'Vertices', P', 'FaceColor' ,'Flat','FaceVertexCData',S_Von_Mises_face);

end

