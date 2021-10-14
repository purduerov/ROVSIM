%% ROV Triton: Thrust Mapper Understanding and Development
%  thrustMapper v2.0
%  Generating the Coarse/Fine Mode Force Envelope Plot for the Tech Doc
%  2021.06.25

close all; clear; clc; tic;

SHADING_INTERP = true;

COM = [0,0,0]; %[in]
tm = thrustMapper(COM);

% Discretize Domain
th = (linspace(0,360,181))'; %(0:5:360)'; %[deg]
phi = (linspace(-90,90,90))'; %[deg](-90:5:90)
n_th = length(th);
n_phi = length(phi);
[TH,PHI] = meshgrid(th,phi);

% Initialize Matrices
F_mag = 500; %[N] BIG
F_des = zeros(n_phi,n_th,3);
F_lim = zeros(n_phi,n_th,3);
F_lim_mag = zeros(n_phi,n_th);

F_lim_totalThrusterUtil = zeros(n_phi,n_th); %SUM of the abs. of all the thrusters at a given time
M_lim_totalThrusterUtil = zeros(n_phi,n_th);

% Iterate Through all the sample directions and run through thrust mapper
for i = 1:n_phi
    for j = 1:n_th
        F_des(i,j,:) = F_mag*[cosd(TH(i,j))*cosd(PHI(i,j)); sind(TH(i,j))*cosd(PHI(i,j)); sind(PHI(i,j))];
        thrustList2 = tm.limitedMap(shiftdim(F_des(i,j,:)),[0;0;0]);
        F_lim(i,j,:) = tm.getForce(thrustList2);
        F_lim_mag(i,j) = norm(shiftdim(F_lim(i,j,:)));
        F_lim_totalThrusterUtil(i,j) = sum(abs(thrustList2));
    end
end

%% Manual Envelope Points
% loop1 = [0, 0, -113.8;
%          -106.93, 0, -59.58;
%          -106.93, 0, 91.36;
%          0, 0, 145.58;
%          106.93, 0, 91.36;
%          106.93, 0, -59.58;
%          0, 0, -113.8];
chain1 = [0, 0, -113.8;
         -106.93, 0, -59.58;
         -106.93, 0, 91.36;
         0, 0, 145.58;
         106.93, 0, 91.36;
         106.93, 0, -59.58];
chain1b = [106.93, 0, -59.58; 0, 0, -113.8];
  
% loop2 = [0, 0, -113.8;
%          0, 38.92, -98.38;
%          0, 38.92, 130.16;
%          0, 0, 145.58;
%          0, -38.92, 130.16;
%          0, -38.92, -98.38;
%          0, 0, -113.8];
chain2 = [0, 0, -113.8;
         0, 38.92, -98.38;
         0, 38.92, 130.16;
         0, 0, 145.58;
         0, -38.92, 130.16];
chain2b = [0, -38.92, 130.16;
           0, -38.92, -98.38;
           0, 0, -113.8];

% loop3 = [106.93, 0, 91.36;
%          0, 38.92, 130.16;
%          -106.93, 0, 91.36;
%          0, -38.92, 130.16;
%          106.93, 0, 91.36];
chain3 = [106.93, 0, 91.36;
         0, 38.92, 130.16;
         -106.93, 0, 91.36;
         0, -38.92, 130.16];
chain3b = [0, -38.92, 130.16; 106.93, 0, 91.36];

% loop4 = [106.93, 0, -59.58;
%          0, 38.92, -98.38;
%          -106.93, 0, -59.58;
%          0, -38.92, -98.38;
%          106.93, 0, -59.58];
chain4 = [106.93, 0, -59.58;
         0, 38.92, -98.38;
         -106.93, 0, -59.58];
chain4b = [-106.93, 0, -59.58;
           0, -38.92, -98.38;
           106.93, 0, -59.58];
    

%% ENVELOPE PLOTS
fineMag = min(min(F_lim_mag));

figure;
hold on; 
ax = gca();
ax.GridAlpha = 0.5;
grid on; box on;
coarse = surf(F_lim_mag.*cosd(TH).*cosd(PHI), F_lim_mag.*sind(TH).*cosd(PHI), F_lim_mag.*sind(PHI), F_lim_totalThrusterUtil, 'DisplayName','Coarse');
%set(coarse,'FaceColor',[1 0 0], 'FaceAlpha',0.5,'EdgeAlpha', 0);
set(coarse,'FaceAlpha',0.5,'EdgeAlpha', 0);
% plot3(loop1(:,1), loop1(:,2), loop1(:,3), 'k--');
% plot3(loop2(:,1), loop2(:,2), loop2(:,3), 'k--');
% plot3(loop3(:,1), loop3(:,2), loop3(:,3), 'k--');
% plot3(loop4(:,1), loop4(:,2), loop4(:,3), 'k--');
plot3(chain1(:,1), chain1(:,2), chain1(:,3), 'k');
plot3(chain1b(:,1), chain1b(:,2), chain1b(:,3), 'k--');
plot3(chain2(:,1), chain2(:,2), chain2(:,3), 'k');
plot3(chain2b(:,1), chain2b(:,2), chain2b(:,3), 'k--');
plot3(chain3(:,1), chain3(:,2), chain3(:,3), 'k');
plot3(chain3b(:,1), chain3b(:,2), chain3b(:,3), 'k--');
plot3(chain4(:,1), chain4(:,2), chain4(:,3), 'k');
plot3(chain4b(:,1), chain4b(:,2), chain4b(:,3), 'k--');
fine = surf(fineMag*cosd(TH).*cosd(PHI), fineMag*sind(TH).*cosd(PHI), fineMag*sind(PHI),'DisplayName','Fine');
%set(fine,'FaceColor',[0 0 1], 'FaceAlpha',1,'EdgeAlpha', 0);
set(fine,'FaceColor',[0 0 1], 'FaceAlpha',1,'EdgeAlpha', 0);
hold off;
daspect([1 1 1]);
view([-140,22]);
xlabel('X [N]');
ylabel('Y [N]');
zlabel('Z [N]');
%legend('show');
colormap jet;
title('Thrust Envelope');
xlim(150*[-1,1]);

fprintf('Runtime: %.4f s\n',toc);
