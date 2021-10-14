%% ROV Triton: Thrust Mapper Understanding and Development
%  thrustMapper v3.0
%  Development Driver 2b: Trialing Force Envelope
%  2021.05.26

%  I'm still developing the nullMap() function

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
x = 300;
F_mag = x; %[N] BIG
F_des = zeros(n_phi,n_th,3);
F_lim = zeros(n_phi,n_th,3);
F_null = zeros(n_phi,n_th,3);
F_lim_mag = zeros(n_phi,n_th);
F_null_mag = zeros(n_phi,n_th);

M_mag = x; %[Nm] BIG
M_des = zeros(n_phi,n_th,3);
M_lim = zeros(n_phi,n_th,3);
M_null = zeros(n_phi,n_th,3);
M_lim_mag = zeros(n_phi,n_th);
M_null_mag = zeros(n_phi,n_th);

F_lim_totalThrusterUtil = zeros(n_phi,n_th); %SUM of the abs. of all the thrusters at a given time
M_lim_totalThrusterUtil = zeros(n_phi,n_th);
F_null_totalThrusterUtil = zeros(n_phi,n_th);
M_null_totalThrusterUtil = zeros(n_phi,n_th);
F_null_numExceed = zeros(n_phi,n_th);
M_null_numExceed = zeros(n_phi,n_th);

% Iterate Through all the sample directions and run through thrust mapper
for i = 1:n_phi
    for j = 1:n_th
        F_des(i,j,:) = F_mag*[cosd(TH(i,j))*cosd(PHI(i,j)); sind(TH(i,j))*cosd(PHI(i,j)); sind(PHI(i,j))];
        thrustList1 = tm.limitedMap(shiftdim(F_des(i,j,:)),[0;0;0]);
        [thrustList2, numExceed] = tm.nullMap(shiftdim(F_des(i,j,:)),[0;0;0]);
        F_lim(i,j,:) = tm.getForce(thrustList1);
        F_null(i,j,:) = tm.getForce(thrustList2);
        F_lim_mag(i,j) = norm(shiftdim(F_lim(i,j,:)));
        F_null_mag(i,j) = norm(shiftdim(F_null(i,j,:)));
        F_lim_totalThrusterUtil(i,j) = sum(abs(thrustList1));
        F_null_totalThrusterUtil(i,j) = sum(abs(thrustList2));
        F_null_numExceed(i,j) = numExceed;
        
        M_des(i,j,:) = M_mag*[cosd(TH(i,j))*cosd(PHI(i,j)); sind(TH(i,j))*cosd(PHI(i,j)); sind(PHI(i,j))];
        thrustList3 = tm.limitedMap([0;0;0],shiftdim(M_des(i,j,:)));
        [thrustList4, numExceed] = tm.nullMap([0;0;0],shiftdim(M_des(i,j,:)));
        M_lim(i,j,:) = tm.getMoment(thrustList3);
        M_null(i,j,:) = tm.getMoment(thrustList4);
        M_lim_mag(i,j) = norm(shiftdim(M_lim(i,j,:)));
        M_null_mag(i,j) = norm(shiftdim(M_null(i,j,:)));
        M_lim_totalThrusterUtil(i,j) = sum(abs(thrustList3));
        M_null_totalThrusterUtil(i,j) = sum(abs(thrustList4));
        M_null_numExceed(i,j) = numExceed;
    end
end

%% ENVELOPE PLOTS

for i = 1:2
    switch(i)
        case 1 %limitedMap()
            F_mag = F_lim_mag;
            M_mag = M_lim_mag;
            F_color = F_lim_totalThrusterUtil;
            M_color = M_lim_totalThrusterUtil;
            mapString = 'limitedMap()';
        case 2 %nullMap()
            F_mag = F_null_mag;
            M_mag = M_null_mag;
            F_color = F_null_totalThrusterUtil;
            M_color = M_null_totalThrusterUtil;
%             F_color = F_null_numExceed;
%             M_color = M_null_numExceed;
            mapString = 'nullMap()';
    end
    
    figure;
    subplot(1,2,1);
    hold on; grid on; box on;
    surf(F_mag.*cosd(TH).*cosd(PHI), F_mag.*sind(TH).*cosd(PHI), F_mag.*sind(PHI), F_color);
    hold off;
    daspect([1 1 1]);
    view([-140,30]);
    xlabel('x [N]');
    ylabel('y [N]');
    zlabel('z [N]');
    title({mapString,'Force Envelope',sprintf('COM: [%d, %d, %d] in',COM(1),COM(2),COM(3)),...
        sprintf('Inscribed Sphere: %.2f N',min(min(F_mag)))});
    h = colorbar;
    ylabel(h,'Total Thruster Utilization [N]','FontSize',11); %Label colorbar
    colormap jet;
    if(SHADING_INTERP) shading('interp'); end

    subplot(1,2,2);
    hold on; grid on; box on;
    surf(M_mag.*cosd(TH).*cosd(PHI), M_mag.*sind(TH).*cosd(PHI), M_mag.*sind(PHI), M_color);
    hold off;
    daspect([1 1 1]);
    view([-140,30]);
    xlabel('Roll (About x) [Nm]');
    ylabel('Pitch (About y) [Nm]');
    zlabel('Yaw (About z) [Nm]');
    title({mapString,'Moment Envelope',sprintf('COM: [%d, %d, %d] in',COM(1),COM(2),COM(3)),...
        sprintf('Inscribed Sphere: %.2f Nm',min(min(M_mag)))});
    h = colorbar;
    ylabel(h,'Total Thruster Utilization [N]','FontSize',11); %Label colorbar
    colormap jet;
    if(SHADING_INTERP) shading('interp'); end
end

%% EXPORT STLs
% surf2stl('ForceEnvelope.stl',F_lim_mag.*cosd(TH).*cosd(PHI), F_lim_mag.*sind(TH).*cosd(PHI), F_lim_mag.*sind(PHI),'ascii');
% surf2stl('MomentEnvelope.stl',M_lim_mag.*cosd(TH).*cosd(PHI), M_lim_mag.*sind(TH).*cosd(PHI), M_lim_mag.*sind(PHI),'ascii');
fprintf('Runtime: %.4f s\n',toc);
