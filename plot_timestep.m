function [ ] = plot_timestep( handle, numEle,Tt, P )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
TSum = zeros(4*numEle,3);
dummy = ones(4*numEle,1);
for i=1:size(Tt,2)
    TSum((i-1)*4+1, :) = [Tt(1,i) Tt(2,i) Tt(3,i)];
    TSum((i-1)*4+2, :) = [Tt(2,i) Tt(3,i) Tt(4,i)];
    TSum((i-1)*4+3, :) = [Tt(3,i) Tt(4,i) Tt(1,i)];
    TSum((i-1)*4+4, :) = [Tt(4,i) Tt(1,i) Tt(2,i)];
end
%figure();
%patch('Faces', TSum, 'Vertices', P', 'FaceColor' ,'Flat','FaceVertexCData',dummy);

set(handle,'Faces',TSum);
set(handle,'Vertices',P');
set(handle,'FaceColor','Flat');
set(handle,'FaceVertexCData',dummy);
drawnow;

end

