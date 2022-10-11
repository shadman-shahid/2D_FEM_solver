clear;
close all;
R1 = [3;4; 0;0.01;0.01;0; 0;0;0.003;0.003];
R2 = [3;4; -0.1;0.1;0.1;-0.1; 0;0;0.1;0.1];
R3 = [3;4; -0.1;0.1;0.1;-0.1; -0.125;-0.125;0.1;0.1];

gd = [R1, R2, R3];
sf = '(R3*R2+R3)+R1';
ns = char('R1','R2','R3');
ns = ns';
d = decsg(gd, sf, ns);
figure
pdegplot(d,'EdgeLabels','on','SubdomainLabels','on')
axis equal
figure
[p,e,t] = initmesh(d);
pdemesh(p,e,t,'ElementLabels','on')
% trimesh(t(1:3,:)', p(1,:), p(2,:))



% rect1 = [3
%     4
%     -1
%     1
%     1
%     -1
%     0
%     0
%     -0.5
%     -0.5];
% C1 = [1
%     1
%     -0.25
%     0.25];
% C2 = [1
%     -1
%     -0.25
%     0.25];
