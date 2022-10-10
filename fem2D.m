% clear all;
% pdetool;

ib = unique([e(1,:),e(2,:)]);
inside_b = unique([e(1,e(6,:)==1),e(2,e(6,:)==1)]);
% ib = unique(e(5,e(5,:)==1));

ne = size(t,2);
np = size(p,2);
ip = size(e,2);

a= zeros(ne,1);
A = zeros(np);
b = zeros(np,1);
for k =1:ne
    i1 = t(1, k);
    i2 = t(2, k);
    i3 = t(3, k);
    
    i123 = [i1,i2,i3];
    
    p1 = p(:, i1);
    p2 = p(:, i2);
    p3 = p(:, i3);
    
    a(k) = abs(det([(p1-p2)'; (p1-p3)']))./2;
    g1 = [(p1-p2)'; (p1-p3)'] \ [1; 1];
    g2 = [(p2-p3)'; (p2-p1)'] \ [1; 1];
    g3 = [(p3-p1)'; (p3-p2)'] \ [1; 1];
    
    g = [g1, g2, g3];
    
    A(i1, i1) = A(i1, i1) + a(k) * dot(g1, g1);
    A(i2, i2) = A(i2, i2) + a(k) * dot(g2, g2);
    A(i3, i3) = A(i3, i3) + a(k) * dot(g3, g3);
    A(i1, i2) = A(i1, i2) + a(k) * dot(g1, g2);
    A(i1, i3) = A(i1, i3) + a(k) * dot(g1, g3);
    A(i2, i1) = A(i2, i1) + a(k) * dot(g2, g1);
    A(i2, i3) = A(i2, i3) + a(k) * dot(g2, g3);
    A(i3, i1) = A(i3, i1) + a(k) * dot(g3, g1);
    A(i3, i2) = A(i3, i2) + a(k) * dot(g3, g2);
%     for i = 1:3
%         for j = 1:3
%             A(i123(i), i123(j)) = A(i123(i), i123(j)) + ...
%                 a(k) * dot(g(:, i), g(:, j));
%         end
%     end
    b(i1) = b(i1) + a(k) ./ 3;
    b(i2) = b(i2) + a(k) ./ 3;
    b(i3) = b(i3) + a(k) ./ 3;    
end


interior = ones(np,1);
interior(ib) = 0;
interior = logical(interior);

inside_boundary = zeros(np, 1);
inside_boundary(inside_b) = 1;
inside_boundary = logical(inside_boundary);
% A = A(interior,interior);
% b = b(interior);
u = zeros(np, 1);
% u(inside_boundary) = 5;
u(interior) = A(interior,interior) \b(interior);



uu2 = zeros(np, 1);
% uu2(inside_boundary) = 5;
A_sp = sparse(A(interior,interior));
[V,D] = eigs(A_sp, length(A(interior,interior)),'sm');

uu2(interior) = V(:,1);

uu = zeros(np, 1);
uu(interior) = 1;
trisurf(t(1:3, :)', p(1,:), p(2,:), u)