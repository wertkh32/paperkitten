function [] = test_ProjDyn()

P0 = [0 10 0 0 10;...
      0 0 10 0 10; ...
      0 0 0 10 10];
Pt = [0 20 0  0  10;...
      0 0  10 0  10; ...
      0 0  0 20  10];  
     
Tt = [1 2; 2 3; 3 4; 4 5];

fixed = [0; 0; 1; 0; 1];
mass = [1;1;1;1;1];
vel = zeros(15,1);
fext = zeros(15,1);
W = [1;1];

numEle = 2;
numNode = 5;
damping = 0.98;


Ptt = Pt;

handle = plot_strain(numEle,Tt,Pt,Pt);
axis([-10,20,-20,20,-10,20]);

M = diag(kron(mass,[1;1;1]));

[A vol] = ProjDyn_Init(numEle,numNode,Tt,P0,M,W,60);

%qn = reshape(Pt, numNode * 3, 1);
h = (1.0/60);

%%Alternating Optimization Approach%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:1000
    Pold = Ptt;
    [Ptt, vel] = ProjDyn_timestep( numEle, numNode, Tt, Ptt, P0, M, A, vol, W, fixed, h, vel, fext );
    Ptt = reshape(Ptt,3,5);
    %disp(vel);
    vel = vel * damping;
    plot_timestep( handle, numEle ,Tt, Ptt );
    pause(1/60);
%disp(ProjDyn_Energy(numEle, numNode, Tt, Ptt, Pold,P0, M, A, vol, W, fixed, h, vel, fext));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%disp(A);
%disp(vol);


end