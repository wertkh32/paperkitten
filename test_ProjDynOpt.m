function [] = test_ProjDynOpt()

P0 = [0 10 0 0 10;...
      0 0 10 0 10; ...
      0 0 0 10 10];
Pt = [0 20 0  0  10;...
      0 0  10 0  10; ...
      0 0  0 20  10];  
     
Tt = [1 2; 2 3; 3 4; 4 5];

fixed = [1; 0; 1; 0; 1];
mass = [1;1;1;1;1];
vel = zeros(15,1);
fext = zeros(15,1);
W = [1;1];

numEle = 2;
numNode = 5;
damping = 0.999;


Ptt = Pt;

handle = plot_strain(numEle,Tt,Pt,Pt);
axis([-10,20,-20,20,-10,20]);

M = diag(kron(mass,[1;1;1]));

[A vol] = ProjDyn_Init(numEle,numNode,Tt,P0,M,W,60);

%qn = reshape(Pt, numNode * 3, 1);
h = (1.0/60);



%%OPTMIZATION%%APPROACH%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pold = Ptt;

disp(ProjDyn_Energy(numEle, numNode, Tt, Ptt, Pold,P0, M, A, vol, W, fixed, h, vel, fext));

options=optimset('LargeScale','off','GradObj','on','DerivativeCheck','off','Display','iter','MaxIter',20,'TolFun',1e-5,'TolX',1e-5);
[Pout, objval] = fminunc(@(x)(ProjDyn_Energy(numEle, numNode, Tt, x, Pold, P0, M, A, vol, W, fixed, h, vel, fext)),Ptt,options);

for i=1:1000
    Pold = Ptt;
    [Ptt, objval] = fminunc(@(x)(ProjDyn_Energy(numEle, numNode, Tt, x, Pold, P0, M, A, vol, W, fixed, h, vel, fext)),Ptt,options);
    vel = reshape( (Ptt - Pold) / h, numNode * 3,1);
    disp(vel);
    %Ptt = reshape(Ptt,numEle,numNode);
    %disp(vel);
    vel = vel * damping;
    plot_timestep( handle, numEle ,Tt, Ptt );
    pause(1/60);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end