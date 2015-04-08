function [ ret ] = project_tri_and_get_area( n, verts )

verts_p = zeros(3,3);
n = n/norm(n);

for i=1:3
a = cross(n,cross(verts(:,i),n));
if(norm(a) < 1e-20)
    ret = 0;
    return;
end

a = a/norm(a);
verts_p(:,i) = dot(a,verts(:,i)) * a;
end

ret = 0.5 * norm(cross(verts_p(:,1),verts_p(:,2)));

end

