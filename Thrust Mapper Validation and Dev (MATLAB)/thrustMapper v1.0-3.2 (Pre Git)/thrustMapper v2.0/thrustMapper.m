%% Purdue IEEE ROV: ROV Triton (X12) Thrust-Mapper
%  Understanding, Validation, and Development
%  Tyler Stagge

%  The purpose of this class and the accompanying driver scripts is to
%  better understand the thrustmapper as a mathematical and practical
%  problem, as well as to validate the math being implemented in the actual
%  python thrustmapper class. It was written independently but concurrent
%  with Scott Hotchkiss' effort to do the necessary fixes/improvements/
%  validation that the Software Team -- whose responsibility this was --
%  failed to do until after we had already wasted 90% of our available pool
%  hours with a horrible, unintuitive thrust mapper. Forgive me for being
%  blunt.

% VERSION HISTORY =========================================================
% [v1.0] [2021.05.09]
%   -A basic, workable version with some plotting functions to help with
%    visualization and debugging.
%   -The mapForce() and mapMoment() functions were written just to confirm
%    to our satisfaction that you cannot, in fact, map force and moments
%    separately. If you do, you'll rotate when you want to translate and
%    translate when you want to rotate.
% [v2.0] [2021.05.10]
%   -Added thruster limiting; the Force and Moment Envelopes can now be
%    generated.

classdef thrustMapper %[v2.0]
    properties
%         X13
%         fowardDist = 7.75;
%         sideDist = 6.75;
%         verticalDist = 3;
%         
%         thrusterCoords_VCSYS = @(x,y,z) [[4.4375, 5.6791, 0];
%                                 [-4.4375, 5.6791, 0];
%                                 [-4.4375, -5.6791, 0];
%                                 [4.4375, -5.6791, 0];
%                                 [7.5, 7.3125, -2.25];
%                                 [-7.5, 7.3125, -2.25];
%                                 [-7.5, -7.3125, -2.25];
%                                 [7.5, -7.3125, -2.25]]*0.0254; %[m] In the Vehicle CSYS
%         thrusterCoords = []; %[in] in the COM CSYS
%         alpha = deg2rad(20); %[rad] Angle of horizontal thrusters
%         beta = 0;
%         thrusterDirections = @(a,b) [[0,0,1];
%                               [0,0,1];
%                               [0,0,1];
%                               [0,0,1];
%                               [cos(a), -sin(a), 0];
%                               [-cos(a), -sin(a), 0];
%                               [-cos(a), sin(a), 0];
%                               [cos(a), sin(a), 0]]; %This is gross, I know

%         X14
%         fowardDist = 8; %7.75
%         sideDist = 6.75; %6.75
%         verticalDist = 8; %6
%         thrusterCoords_VCSYS = @(x,y,z) [[x, y, z];
%                                 [-x, y, z];
%                                 [-x, -y, z];
%                                 [x, -y, z];
%                                 [x, y, -z];
%                                 [-x, y, -z];
%                                 [-x, -y, -z];
%                                 [x, -y, -z]]*0.0254; %[m] In the Vehicle CSYS
%         thrusterCoords = []; %[in] in the COM CSYS
%         alpha = deg2rad(30); %[rad] Yaw of thrusters
%         beta = deg2rad(20); %[rad] Pitch of thrusters
%         thrusterDirections = @(a,b) [[cos(a)*cos(b), -sin(a)*cos(b), -sin(b)];
%                               [-cos(a)*cos(b), -sin(a)*cos(b), -sin(b)];
%                               [-cos(a)*cos(b), sin(a)*cos(b), -sin(b)];
%                               [cos(a)*cos(b), sin(a)*cos(b), -sin(b)];
%                               [cos(a)*cos(b), -sin(a)*cos(b), sin(b)];
%                               [-cos(a)*cos(b), -sin(a)*cos(b), sin(b)];
%                               [-cos(a)*cos(b), sin(a)*cos(b), sin(b)];
%                               [cos(a)*cos(b), sin(a)*cos(b), sin(b)]]; %This is gross, I know

      X15
        fowardDist = 8; %9
        sideDist = 7.25; %6
        verticalDist = 4; %4 up, 4.25 down
        thrusterCoords_VCSYS = @(x,y,z) [[x, y, z];
                                [-x, y, z];
                                [-x, -y, z];
                                [x, -y, z];
                                [x, y, -z];
                                [-x, y, -z];
                                [-x, -y, -z];
                                [x, -y, -z]]*0.0254; %[m] In the Vehicle CSYS
        thrusterCoords = []; %[in] in the COM CSYS
        alpha = deg2rad(30); %[rad] Yaw of thrusters
        beta = deg2rad(25); %[rad] Pitch of thrusters
        thrusterDirections = @(a,b) [[cos(a)*cos(b), -sin(a)*cos(b), -sin(b)];
                              [-cos(a)*cos(b), -sin(a)*cos(b), -sin(b)];
                              [-cos(a)*cos(b), sin(a)*cos(b), -sin(b)];
                              [cos(a)*cos(b), sin(a)*cos(b), -sin(b)];
                              [cos(a)*cos(b), -sin(a)*cos(b), sin(b)];
                              [-cos(a)*cos(b), -sin(a)*cos(b), sin(b)];
                              [-cos(a)*cos(b), sin(a)*cos(b), sin(b)];
                              [cos(a)*cos(b), sin(a)*cos(b), sin(b)]]; %This is gross, I know

        COM_Coords = [0,0,0]; %[m] %Initialization line (real data fed into contstructor)

        map_T2M = []; %ThrusterList-to-Moment Map
        map_T2F = []; %ThrusterList-to-Force Map
        MAP_T2V = []; %Combined ThrusterList-to-Force/Moment(Vector) Map
        map_M2T = []; %Moment-to-ThrusterList Map
        map_F2T = []; %Force-to-ThrusterList Map
        MAP_V2T = []; %Combined Force/Moment(Vector)-to-ThrusterList Map
        
        % T200 Thrust Limits @ 12 V, per Blue Robitics
        THRUST_LIM_FORWARD = 3.71*(9.81); %[N]
        THRUST_LIM_BACKWARD = -2.9*(9.81); %[N]
    end
    
    methods
        %% CONSTRUCTOR METHOD
        %   COM = [in] Center of Mass Coordinates
        function self = thrustMapper(COM)
            self.COM_Coords = COM*0.0254;
            self.thrusterDirections = self.thrusterDirections(self.alpha, self.beta);
            
            % Coordiate Transformation
            self.thrusterCoords = self.thrusterCoords_VCSYS(self.fowardDist, self.sideDist, self.verticalDist) - self.COM_Coords;
            
            % Moment Matrix
            self.map_T2M = zeros(8,3);
            for i = 1:8
                self.map_T2M(i,:) = cross(self.thrusterCoords(i,:), self.thrusterDirections(i,:));
            end
            self.map_T2M = self.map_T2M'; %Transpose to 3x8
            
            % Force Matrix
            self.map_T2F = self.thrusterDirections';
            
            % Combined Matrix
            self.MAP_T2V = [self.map_T2F; self.map_T2M];
            
            % Pseudo-Inverses
            self.map_F2T = pinv(self.map_T2F);
            self.map_M2T = pinv(self.map_T2M);
            self.MAP_V2T = pinv(self.MAP_T2V);
        end
        
        %% MAPPING FUNCTIONS
        
        %------------------------------------------------------------------
        % No thrust limiting or scaling
        function thrustList = dumbMap(self,F_des,M_des)
            justForce = self.MAP_V2T*[F_des;0;0;0];
            justMoment = self.MAP_V2T*[0;0;0;M_des];
            thrustList = justForce + justMoment;
        end
        
        %------------------------------------------------------------------
        % Applies thruster limits
        function thrustList = limitedMap(self,F_des,M_des)
            thrustList = self.dumbMap(F_des,M_des);
            minThrust = min(thrustList);
            if(minThrust < self.THRUST_LIM_BACKWARD)
                thrustList = thrustList*(self.THRUST_LIM_BACKWARD/minThrust);
            end
            maxThrust = max(thrustList);
            if(maxThrust > self.THRUST_LIM_FORWARD)
                thrustList = thrustList*(self.THRUST_LIM_FORWARD/maxThrust);
            end
        end
        
        %% RECOVERY FUNCTIONS
        function force = getForce(self, thrustList)
            force = self.map_T2F*thrustList;
        end
        function moment = getMoment(self, thrustList)
            moment = self.map_T2M*thrustList;
        end
        
        %% PLOTTING FUNCTIONS
        function plotSetup(self)
            FORCE_SCALE = 4;
            TEXT_OFFSET = 1;
            figure;
            hold on; grid on; box on;
            % Vehicle CSYS
            quiver3(0,0,0,FORCE_SCALE,0,0,'r','LineWidth',1.5,'HandleVisibility','off'); % x-axis
            quiver3(0,0,0,0,FORCE_SCALE,0,'g','LineWidth',1.5,'HandleVisibility','off'); % y-axis
            quiver3(0,0,0,0,0,FORCE_SCALE,'b','LineWidth',1.5,'HandleVisibility','off'); % z-axis
            % Thruster Vectors
            for i = 1:8
                x = self.thrusterCoords_VCSYS(i,1)/0.0254;
                y = self.thrusterCoords_VCSYS(i,2)/0.0254;
                z = self.thrusterCoords_VCSYS(i,3)/0.0254;
                u = FORCE_SCALE*self.thrusterDirections(i,1);
                v = FORCE_SCALE*self.thrusterDirections(i,2);
                w = FORCE_SCALE*self.thrusterDirections(i,3);
                plot3(x,y,z, 'ko','MarkerFaceColor','k','HandleVisibility','off');
                text(x+sign(x)*TEXT_OFFSET, y+sign(y)*TEXT_OFFSET, z,sprintf('%d',i));
                quiver3(x,y,z, u,v,w, 'b','LineWidth',1.5,'HandleVisibility','off');
            end
            % COM
            plot3(self.COM_Coords(1),self.COM_Coords(2),self.COM_Coords(3),'ko','HandleVisibility','off');
            plot3(self.COM_Coords(1),self.COM_Coords(2),self.COM_Coords(3),'k+','HandleVisibility','off');
            hold off;
            daspect([1 1 1]);
            maxLim = max([xlim,ylim]);
            xlim(maxLim*[-1,1]);
            ylim(maxLim*[-1,1]);
            view([-90,90]);
            xlabel('x [in]');
            ylabel('y [in]');
            zlabel('z [in]');
            title('Thruster Setup');
        end
        
        %------------------------------------------------------------------
        function plotThrusters(self, thrustList, titleString)
            if(nargin<3) titleString=''; end
            FORCE_SCALE = 0.5;
            MOMENT_SCALE = 1;
            TEXT_OFFSET = 1;

            netForce = self.map_T2F*thrustList;
            netMoment = self.map_T2M*thrustList;
            
            %figure;
            hold on; grid on; box on;
            % Thruster Vectors
            for i = 1:8
                x = self.thrusterCoords_VCSYS(i,1)/0.0254;
                y = self.thrusterCoords_VCSYS(i,2)/0.0254;
                z = self.thrusterCoords_VCSYS(i,3)/0.0254;
                % Plot limits
                limF = FORCE_SCALE*self.thrusterDirections(i,:) * self.THRUST_LIM_FORWARD;
                limB = FORCE_SCALE*self.thrusterDirections(i,:) * self.THRUST_LIM_BACKWARD;
                quiver3(x,y,z, limF(1),limF(2),limF(3), 'LineStyle','--','Color',0.6*[1,1,1],'LineWidth',1.5,'HandleVisibility','off','AutoScale','off');
                quiver3(x,y,z, limB(1),limB(2),limB(3), 'LineStyle','--','Color',0.6*[1,1,1],'LineWidth',1.5,'HandleVisibility','off','AutoScale','off');
                % Plot Current Thrust
                vec = FORCE_SCALE*self.thrusterDirections(i,:) .* thrustList(i);
                quiver3(x,y,z, vec(1),vec(2),vec(3), 'b','LineWidth',1.5,'HandleVisibility','off','AutoScale','off');
                
                plot3(x,y,z, 'ko','MarkerFaceColor','k','HandleVisibility','off');
                text(x+sign(x)*TEXT_OFFSET, y+sign(y)*TEXT_OFFSET, z,sprintf('%d',i));
            end
            % COM, Net Force, and Net Moment
            comx = self.COM_Coords(1);
            comy = self.COM_Coords(2);
            comz = self.COM_Coords(3);
            quiver3(comx,comy,comz, FORCE_SCALE*netForce(1), FORCE_SCALE*netForce(2), FORCE_SCALE*netForce(3),...
                'r','LineWidth',1.5,'DisplayName','Net Force','AutoScale','off');
            quiver3(comx,comy,comz, MOMENT_SCALE*netMoment(1), MOMENT_SCALE*netMoment(2), MOMENT_SCALE*netMoment(3),...
                'm','LineWidth',1.5,'DisplayName','Net Moment','AutoScale','off');
            plot3(comx,comy,comz,'ko','HandleVisibility','off');
            plot3(comx,comy,comz,'k+','HandleVisibility','off');
            
            hold off;
            legend('show');
            daspect([1 1 1]);
            maxLim = max([xlim,ylim]);
            %xlim(maxLim*[-1,1]);
            %ylim(maxLim*[-1,1]);
            view([-90,90]);
            xlabel(sprintf('x [in or %.2f N or %.2f Nm]',1/FORCE_SCALE,1/MOMENT_SCALE));
            ylabel(sprintf('y [in or %.2f N or %.2f Nm]',1/FORCE_SCALE,1/MOMENT_SCALE));
            zlabel(sprintf('z [in or %.2f N or %.2f Nm]',1/FORCE_SCALE,1/MOMENT_SCALE));
            title({titleString,sprintf('Net Force: [%.1f, %.1f, %.1f] N',netForce(1),netForce(2),netForce(3)),...
                sprintf('Net Moment: [%.1f, %.1f, %.1f] Nm',netMoment(1),netMoment(2),netMoment(3))});
        end
    end
end