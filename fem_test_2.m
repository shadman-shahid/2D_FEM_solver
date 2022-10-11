% clear all;
% pdetool;
%% Load the geometry of our desired structure
load geometry_description.mat

%% Get [p,e,t] mesh data
d = decsg(gd, sf, ns);
figure
pdegplot(d,'EdgeLabels','on','SubdomainLabels','on')
axis equal
figure
[p,e,t] = initmesh(d);
pdemesh(p,e,t,'NodeLabels','on')

boundary_value_1 = 1;

% Defining the boundaries and known node values. (Dirichelet Condition)
dirichelet_elem = t(4,:)==1;
metal_nodes = unique(t(1:3, dirichelet_elem));

% boundary = unique([e(1,e(5,:)==13),e(2,e(5,:)==13)]);
dom_bound = e(6:7,:)~=[1;2] & e(6:7,:)~=[2;1];
dom_bound = dom_bound(1,:) | dom_bound(2,:);
% dom_bound = dom_bound & e(6,:)~=0;

inside_b = unique([e(1,e(6,:)~=0 & dom_bound),e(2,e(6,:)~=0 & dom_bound)]);
outside_b = unique([e(1,e(6,:)==0),e(2,e(6,:)==0)]);
% domain_b
% ib = unique(e(5,e(5,:)==1));

n_dom = max(unique(t(4,:)));
ne = size(t,2);
np = size(p,2);
ip = size(e,2);

a= zeros(ne,1);
A = zeros(np);
b = zeros(np,1);
eps = ones(ne, 2);
phi = zeros(np,1);

inside_boundary = zeros(np, 1);
inside_boundary(inside_b) = 1;
inside_boundary = logical(inside_boundary);

eps(t(4,:) == 1,:) = 122.9;
eps(t(4,:) == 2,:) = 1;

for sub_domain = 1:n_dom
    [~,id] = find(t(4,:)==sub_domain);
%     p_sub_dom = unique(t(1:3,id));
    for k = 1:length(id)
        i1 = t(1, id(k));
        i2 = t(2, id(k));
        i3 = t(3, id(k));

        i123 = [i1,i2,i3];

        p1 = p(:, i1);
        p2 = p(:, i2);
        p3 = p(:, i3);

        a(id(k)) = abs(det([(p1-p2)'; (p1-p3)']))./2;
        g1 = [(p1-p2)'; (p1-p3)'] \ [1; 1];
        g2 = [(p2-p3)'; (p2-p1)'] \ [1; 1];
        g3 = [(p3-p1)'; (p3-p2)'] \ [1; 1];

        g = [g1, g2, g3];

        A(i1, i1) = A(i1, i1) + a(id(k)) * dot(eps(k,:)'.*g1, g1);
        A(i2, i2) = A(i2, i2) + a(id(k)) * dot(eps(k,:)'.*g2, g2);
        A(i3, i3) = A(i3, i3) + a(id(k)) * dot(eps(k,:)'.*g3, g3);
        A(i1, i2) = A(i1, i2) + a(id(k)) * dot(eps(k,:)'.*g1, g2);
        A(i1, i3) = A(i1, i3) + a(id(k)) * dot(eps(k,:)'.*g1, g3);
        A(i2, i1) = A(i2, i1) + a(id(k)) * dot(eps(k,:)'.*g2, g1);
        A(i2, i3) = A(i2, i3) + a(id(k)) * dot(eps(k,:)'.*g2, g3);
        A(i3, i1) = A(i3, i1) + a(id(k)) * dot(eps(k,:)'.*g3, g1);
        A(i3, i2) = A(i3, i2) + a(id(k)) * dot(eps(k,:)'.*g3, g2);
        % imposing non-zero dirichelet condition on A matrix
    %     A(inside_b, :) = 0;
    %     A(:, inside_b) = 0;
    %     A(inside_b, inside_b) = eye(length(inside_b));
    %     for i = 1:3
    %         for j = 1:3
    %             A(i123(i), i123(j)) = A(i123(i), i123(j)) + ...
    %                 a(k) * dot(g(:, i), g(:, j));
    %         end
    %     end
        b(i1) = b(i1) + a(id(k)) ./ 3;
        b(i2) = b(i2) + a(id(k)) ./ 3;
        b(i3) = b(i3) + a(id(k)) ./ 3;

    end
%     A(p_sub_dom, p_sub_dom) =  A(p_sub_dom, p_sub_dom).* eps(sub_domain);
end

interior = ones(np,1);
% interior(dom_bound) = 0;
interior(outside_b) = 0;
interior(inside_b) = 0;
% interior(boundary) = 0;
interior = logical(interior);

%% b(interior) = b(interior) - A(interior,dirichelet) * U(dirichelet)
b_modifier = A(interior,[inside_b,outside_b]) * ...
    [ones(length(inside_b),1) .* boundary_value_1; zeros(length(outside_b),1)];
b(interior) = b(interior) - b_modifier;

%% 
% Continuity at a boundary

u = zeros(np, 1);
%% Imposing Dirichelet BC on U
u(inside_b) = boundary_value_1;
%% u(inside_boundary) = 1.5;
u(interior) = A(interior,interior) \ b(interior);



% uu2 = zeros(np, 1);
% % uu2(inside_boundary) = 5;
% A_sp = sparse(A(interior,interior));
% [V,D] = eigs(A_sp, length(A(interior,interior)),'sm');
% 
% uu2(interior) = V(:,1);

uu = zeros(np, 1);
uu(interior) = 1;
trisurf(t(1:3, :)', p(1,:), p(2,:), u)