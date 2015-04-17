P0 = [0 10 0 0;...
      0 0 10 0; ...
      0 0 0 10];
Pt = [0 20 0  0;...
      0 0  10 0; ...
      0 0  0 20];  
     
Tt = [1; 2; 3; 4];


disp(jacobianest(@(x)(rotation_derv_test(Tt, Pt, P0,1,1,x)),0));
%disp(hessian(@(x)(rotation_derv_test(Tt, Pt, P0,1,1,x)),0));
[r, dr] = rotation_derv_test(Tt, Pt, P0,1,1,0);
disp(dr);
%disp(d2r);