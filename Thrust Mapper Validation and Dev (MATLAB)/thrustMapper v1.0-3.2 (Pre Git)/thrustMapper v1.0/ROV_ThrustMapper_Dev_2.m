%% ROV Triton: Thrust Mapper Understanding and Development
%  Development Driver 2: Plotting Erroneous Moments
%  2021.05.09

close all; clear; clc;

%tm = thrustMapper([1,1,0]);
tm = thrustMapper([0,0,0]);

th = (0:5:360)'; %[deg]
phi = (-90:5:90)'; %[deg]
n_th = length(th);
n_phi = length(phi);
[TH,PHI] = meshgrid(th,phi);

F_des = zeros(n_phi,n_th,3);
Mres_mapForce = zeros(n_phi,n_th,3);
Mres_mapBoth = zeros(n_phi,n_th,3);

for i = 1:n_phi
    for j = 1:n_th
        F_des(i,j,:) = [cosd(TH(i,j))*cosd(PHI(i,j)); sind(TH(i,j))*cosd(PHI(i,j)); sind(PHI(i,j))];
        thrustList = tm.mapForce(shiftdim(F_des(i,j,:)));
        Mres_mapForce(i,j,:) = tm.map_T2M*thrustList;
        
        thrustList = tm.mapBoth(shiftdim(F_des(i,j,:)),[0;0;0]);
        Mres_mapBoth(i,j,:) = tm.map_T2M*thrustList;
    end
end

figure;
subplot(1,3,1);
hold on; grid on; box on;
plot3(F_des(:,:,1),F_des(:,:,2),F_des(:,:,3),'r*');
hold off;
daspect([1 1 1]);
xlabel('x');
ylabel('y');
zlabel('z');
title('Desired Force Vectors');

subplot(1,3,2);
plot3(Mres_mapForce(:,:,1),Mres_mapForce(:,:,2),Mres_mapForce(:,:,3),'g.');
grid on; box on;
daspect([1 1 1]);
xlabel('x');
ylabel('y');
zlabel('z');
title('Erroneous Moments from mapForce()');

subplot(1,3,3);
plot3(Mres_mapBoth(:,:,1),Mres_mapBoth(:,:,2),Mres_mapBoth(:,:,3),'g.');
grid on; box on;
%daspect([1 1 1]);
xlabel('x');
ylabel('y');
zlabel('z');
title('Erroneous Moments from mapBoth()');
