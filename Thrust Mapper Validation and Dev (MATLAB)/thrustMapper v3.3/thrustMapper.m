%% Purdue IEEE ROV: ROV Triton (X12) Thrust-Mapper, v3.3
%  Understanding, Validation, and Development
%  Tyler Stagge

% INTRODUCTION ============================================================
%  The purpose of this class and the accompanying driver scripts is to
%  better understand the thrustmapper as a mathematical and practical
%  problem, as well as to validate the math being implemented in the actual
%  Python thrustmapper class. It was written independently but concurrent
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
% [v3.0] [2021.05.23]
%   -Adding the ability to 'traverse the nullspace'. The 6x8 MAP_T2V has a
%    rank of 6 and a nullity of 2 (i.e. the nullspace has a dimension of
%    2). Using MATLAB's null() function, you can get the two basis vectors
%    that describe this nullspace. You can, as such, add any linear
%    combination of these vectors to any thrustList (output of the
%    pseudo-inverse map MAP_V2T) without changing the net force or moment.
% [v3.1] [2021.05.27]
%   -Starting over on nullMap(); its no longer recursive. Its now iterative
%    and calls nullTraverse() to do the actual traversal.
%   -nullTraverse() isnt very good
%   -nullTraverse2() collapses the cases by just checking the dimension and
%    rank of the B matrix for every case
%   -Future improvements:
%       > Rewrite nullMap() to use a binary search to more accurately find
%         the boundary of the thrust envelope
% [v3.2] [2021.05.28]
%   -Deleted the original nullTraverse() and renamed nullTraverse2()
%   -Deleted nullMap_DEFUNCT() that was written for v3.0
%   -Rewrote nullMap(). It now uses a binary search to hone in on the true
%    thrust envelope
%   -Eliminated the unused maps
%   -ADDED TO GITHUB (the ROVSIM library)
% [v3.3] [2021.05.28]
%   -Trying to make the rankB=1 case of nullTraverse() more mathematically
%    complete (trying to catch any potential cases where option1 and
%    option2 aren't viable, but some other (a1,a2) will be viable).

classdef thrustMapper %[v3.3]
    properties
        thrusterCoords_VCSYS = [[4.4375, 5.6791, 0];
                                [-4.4375, 5.6791, 0];
                                [-4.4375, -5.6791, 0];
                                [4.4375, -5.6791, 0];
                                [7.5, 7.3125, -2.25];
                                [-7.5, 7.3125, -2.25];
                                [-7.5, -7.3125, -2.25];
                                [7.5, -7.3125, -2.25]]*0.0254; %[m] In the Vehicle CSYS
        thrusterCoords = []; %[in] in the COM CSYS
        alpha = deg2rad(20); %[rad] Angle of horizontal thrusters
        thrusterDirections = @(a) [[0,0,1];
                              [0,0,1];
                              [0,0,1];
                              [0,0,1];
                              [cos(a), -sin(a), 0];
                              [-cos(a), -sin(a), 0];
                              [-cos(a), sin(a), 0];
                              [cos(a), sin(a), 0]]; %This is gross, I know
        COM_Coords = [0,0,0]; %[m]
        map_T2M = []; %[3x8] ThrusterList-to-Moment Map
        map_T2F = []; %[3x8] ThrusterList-to-Force Map
        MAP_T2V = []; %[6x8] Combined ThrusterList-to-Force/Moment(Vector) Map
        MAP_V2T = []; %[8x6] Combined Force/Moment(Vector)-to-ThrusterList Map
        nullBasis = []; %Set of basis vectors describing the nullspace of MAP_T2V
                        % -See version history notes for explanation of the significance.
        nb1 = []; % nullBasis = [nb1, nb2]
        nb2 = [];
        
        % T200 Thrust Limits @ 12 V, per Blue Robitics
        THRUST_LIM_FORWARD = 3.71*(9.81); %[N]
        THRUST_LIM_BACKWARD = -2.9*(9.81); %[N]
        
        % OTHER CONSTANTS
        TOL_RANK = 1e-12; %The tolerance I give to MATLAB's rank() function
        TOL_BINARY_SEARCH = 1e-5; %Termination tolerance for nullMap()
        ITERATION_LIMIT = 20;     %Termination criterion for nullMap()
    end
    
    methods
        %% CONSTRUCTOR METHOD =============================================
        %   COM = [in] Center of Mass Coordinates
        function self = thrustMapper(COM)
            self.COM_Coords = COM*0.0254;
            self.thrusterDirections = self.thrusterDirections(self.alpha);
            
            % Coordiate Transformation
            self.thrusterCoords = self.thrusterCoords_VCSYS - self.COM_Coords;
            
            % Moment Matrix
            self.map_T2M = zeros(8,3);
            for i = 1:8
                self.map_T2M(i,:) = cross(self.thrusterCoords(i,:), self.thrusterDirections(i,:));
            end
            self.map_T2M = self.map_T2M'; %Transpose to 3x8
            
            % Force Matrix
            self.map_T2F = self.thrusterDirections';
            
            % Combined Matrix
            self.MAP_T2V = [self.map_T2F; self.map_T2M]; %[6x8], rank=6, nullity=2
            self.nullBasis = null(self.MAP_T2V); %Two [8x1] vectors
            self.nb1 = self.nullBasis(:,1);
            self.nb2 = self.nullBasis(:,2);
            
            % Pseudo-Inverse
            self.MAP_V2T = pinv(self.MAP_T2V);
        end
        
        %% THRUST MAPPING FUNCTIONS =======================================
        
        %------------------------------------------------------------------
        % No thrust limiting or scaling
        function thrustList = dumbMap(self,F_des,M_des)
            justForce = self.MAP_V2T*[F_des;0;0;0];  %There is no actual reason to do these separately
            justMoment = self.MAP_V2T*[0;0;0;M_des];
            thrustList = justForce + justMoment;
        end
        
        %------------------------------------------------------------------
        % Applies thruster limits
        function thrustList = limitedMap(self,F_des,M_des)
            thrustList = self.MAP_V2T*[shiftdim(F_des); shiftdim(M_des)]; %Map
            minThrust = min(thrustList);
            if(minThrust < self.THRUST_LIM_BACKWARD)
                thrustList = thrustList*(self.THRUST_LIM_BACKWARD/minThrust);
            end
            maxThrust = max(thrustList);
            if(maxThrust > self.THRUST_LIM_FORWARD)
                thrustList = thrustList*(self.THRUST_LIM_FORWARD/maxThrust);
            end
        end
        
        %------------------------------------------------------------------
        % Uses nullspace-traversal to extend the thrust envelope. This
        %  version uses a binary search to hone in on the true boundary
        function thrustList = nullMap(self,F_des,M_des)
            V_des = [shiftdim(F_des); shiftdim(M_des)];
            [isDesiredThrustPossible, thrustList] = self.isThrustPossible(V_des);
            if(~isDesiredThrustPossible)
                terminateSearch = false;
                numIterations = 1;
                upper = V_des;      upperPossible = false;
                lower = zeros(6,1); lowerPossible = true;
                while(~terminateSearch)
                    middle = 0.5*(upper+lower);
                    [middlePossible,thrustList] = self.isThrustPossible(middle);
                    if(middlePossible && ~upperPossible) % search next between middle and upper
                        lower = middle; lowerPossible = middlePossible;
                    elseif(~middlePossible && lowerPossible) % search next between middle and lower
                        upper = middle; upperPossible = middlePossible;
                    else
                        fprintf('WARNING: This branch of the binary search shouldn''t be possible');
                    end
                    numIterations = numIterations + 1;
                    
                    % Check termination criterion
                    if(numIterations >= self.ITERATION_LIMIT)
                        terminateSearch = true;
                        %fprintf('%d iterations of binary search exceeded.\n',self.ITERATION_LIMIT);
                    end
                    if(abs(norm(lower)-norm(upper))/norm(lower) < self.TOL_BINARY_SEARCH)
                        terminateSearch = true;
                        %fprintf('Termination tolerance met after %d iterations.\n',numIterations);
                    end
                end
            end
        end
        
        %% SUPPORT FUNCTIONS ==============================================
        
        %------------------------------------------------------------------
        % Returns a 8x1 vector whose elements are only {-1,0,+1}
        %   +1 means the forward limit for that thruster is exceeded
        %   -1 means the backward limit for that thruster is exceeded
        function exceedsLimit = whichThrustersExceedLimits(self,thrustList)
            exceedsLimit = zeros(8,1);
            for i = 1:8
                if(thrustList(i) < self.THRUST_LIM_BACKWARD)
                    exceedsLimit(i) = -1;
                elseif(thrustList(i) > self.THRUST_LIM_FORWARD)
                    exceedsLimit(i) = +1;
                end
            end
        end
        
        %------------------------------------------------------------------
        % Checks if a requested thrust is at all possible without exceeding
        % thruster limits -- uses nullTraverse() as part of that check
        function [isPossible, thrustList] = isThrustPossible(self,V)
            thrustList = self.MAP_V2T*V; % V = [F;M]
            exceedsLimitVector = self.whichThrustersExceedLimits(thrustList);
            numExceedsLimit = sum(abs(exceedsLimitVector)); %1-norm
            if(numExceedsLimit>0) %If the pseudo-inverse doesnt give a valid option, try traversing the nullspace
                thrustList = self.nullTraverse(thrustList,exceedsLimitVector);
                numExceedsLimit = sum(abs(self.whichThrustersExceedLimits(thrustList))); %1-norm
                if(numExceedsLimit>0)
                    isPossible = false;
                else
                    isPossible = true;
                end
            else
                isPossible = true;
            end
        end
        
        %------------------------------------------------------------------
        % If a given thrustList exceeds thruster limits, this functions
        % checks to see if adding a linear combo of the nullspace basis
        % vectors can achieve a thrustList with the same net force and
        % moment but that doesn't exceed thruster limits (which is not 
        % always possible)
        %  -This version collapsed a bunch of the cases (and does things
        %   way better) by comparing the rank of the B matrix with the rank
        %   of the augmented matrix to figure out if the system is
        %   overconstrained or underconstrained
        function thrustList = nullTraverse(self,thrustList,exceedsLimitVector)
            % [B]{a} = {deltaT} = {T_lim - T_exceed}
            %   where {a} is [a1;a2], the coefficients by which the 
            %   nullspace basis vectors {nb1,nb2} get multiplied
            T_exceed = []; % sub-set of thrustList (just the ones that exceed)
            T_lim = [];    % the corresponding positive or negative limit
            B = [];        % matrix of the relevant components of the nb1,nb2 vectors
            for i = 1:8
                if(exceedsLimitVector(i) ~= 0)
                    T_exceed = [T_exceed; thrustList(i)];
                    B = [B; self.nb1(i), self.nb2(i)];
                    limit = (exceedsLimitVector(i)==1)*self.THRUST_LIM_FORWARD...
                            +(exceedsLimitVector(i)==-1)*self.THRUST_LIM_BACKWARD;
                    T_lim = [T_lim; limit];
                end
            end
            deltaT = T_lim - T_exceed;
            augmentedMatrix = [B, deltaT];
            rankB = rank(B,self.TOL_RANK);
            rankAug = rank(augmentedMatrix,self.TOL_RANK);
            
            % If rank(B) == rank(augmentedMatrix), then the system is NOT
            %  overconstrained
            if(rankB==rankAug)
                if(rankB==1) %Either only one eqn. or one eqn. with redundant copies -- can just choose one
                    % To make it easy, I'm artificially limiting myself to
                    %  adding coeff*nb1 OR coeff*nb2 (i.e. 2 possible
                    %  solutions instead of infinite solutions). Thus, 
                    %  there may be  cases where the requested thrust IS in
                    %  fact possible w/o exceeding limits, but this 
                    %  algorithm won't find it.
                    a1_intercept = (deltaT(1)/B(1,1));
                    a2_intercept = (deltaT(1)/B(1,2));
                    if(true) %DEBUG PLOT
                        a2 = @(a1) -(a2_intercept/a1_intercept)*a1 + a2_intercept;
                        a1_n = abs(a1_intercept)*[-5,5];
                        figure;
                        plot(a1_n, a2(a1_n), 'r', 'LineWidth',1.5);
                        grid on;
                        xlabel('a_1');
                        ylabel('a_2');
                    end
                    option1 = thrustList + a1_intercept*self.nb1;
                    option2 = thrustList + a2_intercept*self.nb2;
                    isOption1Viable = sum(abs(self.whichThrustersExceedLimits(option1))) < self.TOL_RANK;
                    isOption2Viable = sum(abs(self.whichThrustersExceedLimits(option2))) < self.TOL_RANK;
                    switch(isOption1Viable + isOption2Viable)
                        case 0 % If neither are viable, do nothing and the
                               %  requested F and M will be scaled down on
                               %  the next iteration
                            %fprintf('This case was reached.\n');
                        case 1 % If only one is viable, use that one
                            thrustList = (isOption1Viable*option1) + (isOption2Viable*option2);
                        case 2 % If both are viable, take the one that uses less total thrust
                            thrustUtil1 = sum(abs(option1)); %i.e. norm(option1,1)
                            thrustUtil2 = sum(abs(option2));
                            thrustList = (thrustUtil1 <= thrustUtil2)*option1 + (thrustUtil1 > thrustUtil2)*option2;
                    end
                elseif(rankB==2) %Either two eqn. or two eqn. plus redundant copies -- can just choose two
                    a = B(1:2,:)\deltaT(1:2);
                    thrustList = thrustList + a(1)*self.nb1 + a(2)*self.nb2;
                else
                    fprintf('WARNING: rankB=rankAug > 2\n');
                end
            
            % If rank(B) < rank(augmentedMatrix), then the system IS
            %  overconstrained
            elseif(rankB < rankAug)
                % Do nothing and F,M will be scaled down on the next
                %  iteration
            else
                fprintf('WARNING: rankB > rankAug. How???\n');
            end
        end
        
        %% RECOVERY FUNCTIONS =============================================
        function force = getForce(self, thrustList)
            force = self.map_T2F*thrustList;
        end
        function moment = getMoment(self, thrustList)
            moment = self.map_T2M*thrustList;
        end
        
        %% PLOTTING FUNCTIONS =============================================
        
        %------------------------------------------------------------------
        % Plots the thruster unit vectors (w/ labels) and the COM in 3D
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
        % Plots the thrusts of each thruster and the net force and moment
        % on the vehicle.
        %   - This has only gotten more awkward
        %   - {titleString, F_des, M_des} are optional arguments used in
        %     generating the title of the plot.
        function plotThrusters(self, thrustList, titleString, F_des, M_des)
            if(nargin<3) titleString=''; end
            if(nargin<5) F_des = false; M_des = false; end
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
            %maxLim = max([xlim,ylim]);
            %xlim(maxLim*[-1,1]);
            %ylim(maxLim*[-1,1]);
            view([-90,90]);
            xlabel(sprintf('x [1 in or %.2f N or %.2f Nm]',1/FORCE_SCALE,1/MOMENT_SCALE));
            ylabel(sprintf('y [1 in or %.2f N or %.2f Nm]',1/FORCE_SCALE,1/MOMENT_SCALE));
            zlabel(sprintf('z [1 in or %.2f N or %.2f Nm]',1/FORCE_SCALE,1/MOMENT_SCALE));
            if(sum(F_des)==0) %Not printing F_des or M_des
                title({titleString,sprintf('Net Force: [%.1f, %.1f, %.1f] N',netForce(1),netForce(2),netForce(3)),...
                    sprintf('Net Moment: [%.1f, %.1f, %.1f] Nm',netMoment(1),netMoment(2),netMoment(3))});
            else
                title({titleString,sprintf('Desired Force: [%.1f, %.1f, %.1f] N',F_des(1),F_des(2),F_des(3)),...
                    sprintf('Desired Moment: [%.1f, %.1f, %.1f] Nm',M_des(1),M_des(2),M_des(3)),...
                    sprintf('Net Force: [%.1f, %.1f, %.1f] N',netForce(1),netForce(2),netForce(3)),...
                    sprintf('Net Moment: [%.1f, %.1f, %.1f] Nm',netMoment(1),netMoment(2),netMoment(3))});
            end
        end
    end
end