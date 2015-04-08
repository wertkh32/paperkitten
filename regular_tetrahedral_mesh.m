function [T,V, Edge] = regular_tetrahedral_mesh(nx,ny,nz)
% REGULAR_TETRAHEDRAL_MESH
%
% [T,V] = regular_tetrahedral_mesh(nx,ny,nz)
%
% Generates a regular tetrahedral mesh with dimensions (nx,ny,nz)
%
% Input:
%   nx  number of points in x direction on grid
%   ny  number of points in y direction on grid
%   nz  number of points in z direction on grid
% Output:
%   T  tetrahedra list of indices into V
%   V  list of vertex coordinates in 3D
%
% Example
%    [T,V] = regular_tetrahedral_mesh(3,3,3);
%    tetramesh(T,V); %also try tetramesh(T,V,'FaceAlpha',0);
%
% See also delaunayn, tetramesh
%

% Create a grid
[x,y,z] = meshgrid(linspace(0,1,nx),linspace(0,1,ny),linspace(0,1,nz));
V = [x(:) y(:) z(:)];
% meshgrid flips x and y ordering
idx = reshape(1:prod([ny,nx,nz]),[ny,nx,nz]);
v1 = idx(1:end-1,1:end-1,1:end-1);v1=v1(:);
v2 = idx(1:end-1,2:end,1:end-1);v2=v2(:);
v3 = idx(2:end,1:end-1,1:end-1);v3=v3(:);
v4 = idx(2:end,2:end,1:end-1);v4=v4(:);
v5 = idx(1:end-1,1:end-1,2:end);v5=v5(:);
v6 = idx(1:end-1,2:end,2:end);v6=v6(:);
v7 = idx(2:end,1:end-1,2:end);v7=v7(:);
v8 = idx(2:end,2:end,2:end);v8=v8(:);
T = [ ...
    v1  v3  v8  v7; ...
    v1  v8  v5  v7; ...
    v1  v3  v4  v8; ...
    v1  v4  v2  v8; ...
    v1  v6  v5  v8; ...
    v1  v2  v6  v8];

E = [ ...
    T(:,1) T(:,2); ...
    T(:,2) T(:,3); ...
    T(:,3) T(:,4); ...
    T(:,4) T(:,1); ...
    T(:,1) T(:,3); ...
    T(:,2) T(:,4)];

% sort rows so that edge are reorder in ascending order of indices
sortedE1 = sort(E,2); 
sortedE = sort(sortedE1,1);

% determine uniqueness of faces
[Edge,m,n] = unique(sortedE,'rows');

% construct sparse matrix such that N(i,j) := 1 only if i and j are neighbors
% otherwise 0
% N = sparse( ...
%     [sortedE(:,1)], ...
%     [sortedE(:,2)], ...
%     ones(size(sortedE,1),1), ...
%     size(V,1), ...
%     size(V,1)) > 0;
end

