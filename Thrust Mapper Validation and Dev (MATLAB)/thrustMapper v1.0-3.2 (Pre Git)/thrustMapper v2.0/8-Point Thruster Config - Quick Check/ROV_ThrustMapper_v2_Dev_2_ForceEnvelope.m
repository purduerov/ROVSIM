%% ROV Triton: Thrust Mapper Understanding and Development
%  thrustMapper v2.0-BRANCH
%  Development Driver 2: Plotting Force Envelope
%  2021.05.09
%  2021.11.15

close all; clear; clc; tic;

SHADING_INTERP = true;

COM = [0,0,0]; %[in]
tm = thrustMapper(COM);

% Discretize Domain
res = 'fine';
switch(res)
    case 'low'
        th = (linspace(0,360,181))'; %(0:5:360)'; %[deg]
        phi = (linspace(-90,90,90))'; %[deg](-90:5:90)
    case 'medium'
        th = (0:1:360)'; %(0:5:360)'; %[deg]
        phi = (-90:1:90)'; %[deg](-90:5:90)
    case 'fine'
        th = (0:0.5:360)'; %[deg]
        phi = (-90:0.5:90)'; %[deg]
end
n_th = length(th);
n_phi = length(phi);
[TH,PHI] = meshgrid(th,phi);

% Initialize Matrices
F_mag = 500; %[N] BIG
F_des = zeros(n_phi,n_th,3);
F_unl = zeros(n_phi,n_th,3);
F_lim = zeros(n_phi,n_th,3);
F_unl_mag = zeros(n_phi,n_th);
F_lim_mag = zeros(n_phi,n_th);

M_mag = 500; %[Nm] BIG
M_des = zeros(n_phi,n_th,3);
M_unl = zeros(n_phi,n_th,3);
M_lim = zeros(n_phi,n_th,3);
M_unl_mag = zeros(n_phi,n_th);
M_lim_mag = zeros(n_phi,n_th);

F_lim_totalThrusterUtil = zeros(n_phi,n_th); %SUM of the abs. of all the thrusters at a given time
M_lim_totalThrusterUtil = zeros(n_phi,n_th);

% Iterate Through all the sample directions and run through thrust mapper
for i = 1:n_phi
    for j = 1:n_th
        F_des(i,j,:) = F_mag*[cosd(TH(i,j))*cosd(PHI(i,j)); sind(TH(i,j))*cosd(PHI(i,j)); sind(PHI(i,j))];
        thrustList1 = tm.dumbMap(shiftdim(F_des(i,j,:)),[0;0;0]);
        thrustList2 = tm.limitedMap(shiftdim(F_des(i,j,:)),[0;0;0]);
        F_unl(i,j,:) = tm.getForce(thrustList1);
        F_lim(i,j,:) = tm.getForce(thrustList2);
        F_unl_mag(i,j) = norm(shiftdim(F_unl(i,j,:)));
        F_lim_mag(i,j) = norm(shiftdim(F_lim(i,j,:)));
        F_lim_totalThrusterUtil(i,j) = sum(abs(thrustList2));
        
        M_des(i,j,:) = M_mag*[cosd(TH(i,j))*cosd(PHI(i,j)); sind(TH(i,j))*cosd(PHI(i,j)); sind(PHI(i,j))];
        thrustList3 = tm.dumbMap([0;0;0],shiftdim(M_des(i,j,:)));
        thrustList4 = tm.limitedMap([0;0;0],shiftdim(M_des(i,j,:)));
        M_unl(i,j,:) = tm.getMoment(thrustList3);
        M_lim(i,j,:) = tm.getMoment(thrustList4);
        M_unl_mag(i,j) = norm(shiftdim(M_unl(i,j,:)));
        M_lim_mag(i,j) = norm(shiftdim(M_lim(i,j,:)));
        M_lim_totalThrusterUtil(i,j) = sum(abs(thrustList4));
    end
end

%% ENVELOPE PLOTS

figure;
subplot(1,2,1);
hold on; grid on; box on;
surf(F_lim_mag.*cosd(TH).*cosd(PHI), F_lim_mag.*sind(TH).*cosd(PHI), F_lim_mag.*sind(PHI), F_lim_totalThrusterUtil);
hold off;
daspect([1 1 1]);
%view([-90,90]);
view([-140,30]);
xlabel('x [N]');
ylabel('y [N]');
zlabel('z [N]');
title({'Force Envelope',sprintf('COM: [%d, %d, %d] in',COM(1),COM(2),COM(3)),...
    sprintf('Inscribed Sphere: %.2f N',min(min(F_lim_mag)))});
h = colorbar;
ylabel(h,'Total Thruster Utilization [N]','FontSize',11); %Label colorbar
colormap jet;
if(SHADING_INTERP) shading('interp'); end

subplot(1,2,2);
hold on; grid on; box on;
surf(M_lim_mag.*cosd(TH).*cosd(PHI), M_lim_mag.*sind(TH).*cosd(PHI), M_lim_mag.*sind(PHI), M_lim_totalThrusterUtil);
hold off;
daspect([1 1 1]);
%view([-90,90]);
view([-140,30]);
xlabel('Roll (About x) [Nm]');
ylabel('Pitch (About y) [Nm]');
zlabel('Yaw (About z) [Nm]');
title({'Moment Envelope',sprintf('COM: [%d, %d, %d] in',COM(1),COM(2),COM(3)),...
    sprintf('Inscribed Sphere: %.2f Nm',min(min(M_lim_mag)))});
h = colorbar;
ylabel(h,'Total Thruster Utilization [N]','FontSize',11); %Label colorbar
colormap jet;
if(SHADING_INTERP) shading('interp'); end

%% EXPORT STLs
% surf2stl('ForceEnvelope.stl',F_lim_mag.*cosd(TH).*cosd(PHI), F_lim_mag.*sind(TH).*cosd(PHI), F_lim_mag.*sind(PHI),'ascii');
% surf2stl('MomentEnvelope.stl',M_lim_mag.*cosd(TH).*cosd(PHI), M_lim_mag.*sind(TH).*cosd(PHI), M_lim_mag.*sind(PHI),'ascii');
fprintf('Runtime: %.4f s\n',toc);
