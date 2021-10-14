%% ROV Triton: Thrust Mapper Understanding and Development
%  Development Driver 3
%  2021.05.09

close all; clear; clc;

MAKE_GIF = 1;
GIF_FILENAME = 'ROV_ThrustMapper_v2_Dev_3_Animation_3.gif';

tm = thrustMapper([0,0,0]);

F_mag = 150; %[N]
th_n = deg2rad(0:1:360)'; %[rad]
n = length(th_n);

FORCE_SCALE = 0.2;
MOMENT_SCALE = 1;
TEXT_OFFSET = 2;

h = figure;
hold on; grid on; box on;
daspect([1 1 1]);
view([-170,10]);
xlabel(sprintf('x [in or %.2f N or %.2f Nm]',1/FORCE_SCALE,1/MOMENT_SCALE));
ylabel(sprintf('y [in or %.2f N or %.2f Nm]',1/FORCE_SCALE,1/MOMENT_SCALE));
zlabel(sprintf('z [in or %.2f N or %.2f Nm]',1/FORCE_SCALE,1/MOMENT_SCALE));
xlim(25*[-1,1]);
ylim(25*[-1,1]);
zlim(15*[-1,1]);
pause(5);
for j = 1:n
    F_des = F_mag*[cos(th_n(j));sin(th_n(j));0];
    thrustList = tm.limitedMap(F_des,[0;0;0]);
    netForce = tm.getForce(thrustList);
    netMoment = tm.getMoment(thrustList);
    
    cla;
    maxForceX = FORCE_SCALE*106.9;
    maxForceY = FORCE_SCALE*37.88;
    plot([-maxForceX,0,maxForceX],[0,maxForceY,0], 'LineWidth',1.5, 'Color', 0.7*[1,1,1]);
    plot([-maxForceX,0,maxForceX],[0,-maxForceY,0], 'LineWidth',1.5, 'Color', 0.7*[1,1,1]);
    % Thruster Vectors
    for i = 1:8
        x = tm.thrusterCoords_VCSYS(i,1)/0.0254;
        y = tm.thrusterCoords_VCSYS(i,2)/0.0254;
        z = tm.thrusterCoords_VCSYS(i,3)/0.0254;
        % Plot limits
        limF = FORCE_SCALE*tm.thrusterDirections(i,:) * tm.THRUST_LIM_FORWARD;
        limB = FORCE_SCALE*tm.thrusterDirections(i,:) * tm.THRUST_LIM_BACKWARD;
        quiver3(x,y,z, limF(1),limF(2),limF(3), 'LineStyle','--','Color',0.6*[1,1,1],'LineWidth',1.5,'AutoScale','off');
        quiver3(x,y,z, limB(1),limB(2),limB(3), 'LineStyle','--','Color',0.6*[1,1,1],'LineWidth',1.5,'AutoScale','off');
        % Plot Current Thrust
        vec = FORCE_SCALE*tm.thrusterDirections(i,:) .* thrustList(i);
        quiver3(x,y,z, vec(1),vec(2),vec(3), 'b','LineWidth',1.5,'AutoScale','off');

        plot3(x,y,z, 'ko','MarkerFaceColor','k');
        %text(x+sign(x)*TEXT_OFFSET, y+sign(y)*TEXT_OFFSET, z,sprintf('%d',i));
    end
    % COM, Net Force, and Net Moment
    comx = tm.COM_Coords(1);
    comy = tm.COM_Coords(2);
    comz = tm.COM_Coords(3);
    quiver3(comx,comy,comz, FORCE_SCALE*netForce(1), FORCE_SCALE*netForce(2), FORCE_SCALE*netForce(3),...
        'r','LineWidth',1.5,'DisplayName','Net Force','AutoScale','off');
    quiver3(comx,comy,comz, MOMENT_SCALE*netMoment(1), MOMENT_SCALE*netMoment(2), MOMENT_SCALE*netMoment(3),...
        'm','LineWidth',1.5,'DisplayName','Net Moment','AutoScale','off');
    plot3(comx,comy,comz,'ko');
    plot3(comx,comy,comz,'k+');
    
    if(MAKE_GIF)
        frame = getframe(h);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if(j==1)
            imwrite(imind,cm, GIF_FILENAME, 'gif', 'Loopcount',inf);
        else
            imwrite(imind,cm, GIF_FILENAME, 'gif', 'WriteMode','append', 'DelayTime', 0.02);
        end
    end
    pause = 0.02;
end