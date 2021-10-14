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


classdef thrustMapper %[v3.1]
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
        map_M2T = []; %[8x3] Moment-to-ThrusterList Map | NO LONGER USED (was only to confirm that MAP_V2T was needed)
        map_F2T = []; %[8x3] Force-to-ThrusterList Map  | NO LONGER USED (was only to confirm that MAP_V2T was needed)
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
            
            % Pseudo-Inverses
            %self.map_F2T = pinv(self.map_T2F); %NO LONGER USED
            %self.map_M2T = pinv(self.map_T2M); %NO LONGER USED
            self.MAP_V2T = pinv(self.MAP_T2V);
        end
        
        %% MAPPING FUNCTIONS ==============================================
        
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
        % Uses nullspace-traversal to extend the thrust envelope (v3.1 version)
        %   >In order to do this smarter, I'd need to do it
        %    semi-recursively. Instead of just scaling F_des and M_des down
        %    by 1% until it works, I should do a binary search to more
        %    accurately find the boundary of the thrust envelope
        function thrustList = nullMap(self,F_des,M_des)
            terminate = false;
            numIterations = 0;
            while(~terminate)
                F = F_des*(1 - 0.01*numIterations); %Decrease linearly by 1% on every iteration
                M = M_des*(1 - 0.01*numIterations);
                thrustList = self.MAP_V2T*[shiftdim(F); shiftdim(M)]; %Map
                
                % Check if limits are exceeded
                exceedsLimitVector = self.whichThrustersExceedLimits(thrustList);
                numExceedsLimit = sum(abs(exceedsLimitVector)); %1-norm
                if(numExceedsLimit>0)
                    thrustList = nullTraverse2(self,thrustList,exceedsLimitVector);
                    numExceedsLimit = sum(abs(self.whichThrustersExceedLimits(thrustList))); %1-norm
                    if(numExceedsLimit==0) terminate = true; end
                else
                    terminate = true;
                end
                numIterations = numIterations + 1;
                assert(numIterations<100, 'Error: F_des and M_des have been scaled all the way to zero');
            end
        end
        
        %------------------------------------------------------------------
        % Private function called by nullMap() -- this version isn't very
        %  good
        function thrustList = nullTraverse(self,thrustList,exceedsLimitVector)
            numExceedsLimit = sum(abs(exceedsLimitVector)); %1-norm
            
            % If 1 thruster exceeds limits
            %   >Compute coefficient for just nb1 and just nb2; choose one
            if(numExceedsLimit==1)
                index = abs((1:8)*exceedsLimitVector); %Picks out the one nonzero element
                limit = (exceedsLimitVector(index)==1)*self.THRUST_LIM_FORWARD...
                        +(exceedsLimitVector(index)==-1)*self.THRUST_LIM_BACKWARD;
                deltaT = limit - thrustList(index);
                option1 = thrustList + (deltaT/self.nb1(index))*self.nb1; %SHIT -- I had this wrong in v3
                option2 = thrustList + (deltaT/self.nb2(index))*self.nb2;
                thrustUtil1 = sum(abs(option1)); %i.e. norm(option1,1)
                thrustUtil2 = sum(abs(option2));
                isOption1Viable = sum(abs(self.whichThrustersExceedLimits(option1)))~=0;
                isOption2Viable = sum(abs(self.whichThrustersExceedLimits(option2)))~=0;
                if(isOption1Viable && (thrustUtil1 < thrustUtil2))
                    thrustList = option1;
                elseif(isOption2Viable)
                    thrustList = option2;
                else
                    % Do nothing and the requested F and M will be scaled
                    %  down on the next iteration
                end
            
            % If 2 thrusters exceed limits
            %   >Set up and solve the fully constrained linear system
            %    for the coefficients on nb1 and nb2
            elseif(numExceedsLimit==2)
                T_exceed = []; % sub-set of thrustList (just the ones that exceed)
                T_lim = [];    % the corresponding positive or negative limit
                B = [];        % matrix -- the relevant components of the nb1,nb2 vectors
                for i = 1:8
                    if(exceedsLimitVector(i) ~= 0)
                        T_exceed = [T_exceed; thrustList(i)];
                        B = [B; self.nb1(i), self.nb2(i)];
                        limit = (exceedsLimitVector(i)==1)*self.THRUST_LIM_FORWARD...
                                +(exceedsLimitVector(i)==-1)*self.THRUST_LIM_BACKWARD;
                        T_lim = [T_lim; limit];
                    end
                end 
                if(abs(det(B)) > 1e-10) %Invertible
                    a = B\(T_lim-T_exceed); %Coefficients for nb1, nb2
                    thrustList = thrustList + a(1)*self.nb1 + a(2)*self.nb2;
                    % ^If this causes a different thruster to exceed a limit,
                    %  then that will be caught by nullMap() and the requested
                    %  F and M will be scaled down on the next iterations
                else %Singular matrix -- need to treat it like one thruster is exceeding
                    % I think it'll only get here if the two rows of B are
                    %  the same (and if the same limits apply?)
                    isWhatIExpectB = (sum(abs(B(1,:)-B(2,:)))<1e-10) + (sum(abs(B(1,:)+B(2,:)))<1e-10);
                    isWhatIExpectT = abs(T_lim(1)-T_lim(2))<1e-10;
                    isWhatIExpect = isWhatIExpectB * isWhatIExpectT;
                    if(~isWhatIExpectB) 
                        fprintf('Not what I expect!\n');
                    end
                    % FIXME: I still need to do the stuff here
                end
            
            % If 3 or more thrusters exceed limits
            elseif(numExceedsLimit>2)
                % Do nothing and the requested F and M will be scaled down
                %  on the next iteration
            else
                assert(true, 'Error: This case should be unreachable.');
            end
        end
        
        %------------------------------------------------------------------
        % This version collapsed a bunch of the cases (and does things
        %  way better) by comparing the rank of the B matrix with the rank
        %  of the augmented matrix to figure out if the system is
        %  overconstrained or underconstrained
        function thrustList = nullTraverse2(self,thrustList,exceedsLimitVector)
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
            
            % If rank(B) == rank(augmentedMatrix), then the system is not
            %  overconstrained
            if(rankB==rankAug)
                if(rankB==1) %Either only one eqn. or one eqn. with redundant copies -- can just choose one
                    option1 = thrustList + (deltaT(1)/B(1,1))*self.nb1;
                    option2 = thrustList + (deltaT(1)/B(1,2))*self.nb2;
                    thrustUtil1 = sum(abs(option1)); %i.e. norm(option1,1)
                    thrustUtil2 = sum(abs(option2));
                    isOption1Viable = sum(abs(self.whichThrustersExceedLimits(option1))) < self.TOL_RANK;
                    isOption2Viable = sum(abs(self.whichThrustersExceedLimits(option2))) < self.TOL_RANK;
                    if(isOption1Viable && (thrustUtil1 < thrustUtil2))
                        thrustList = option1;
                    elseif(isOption2Viable)
                        thrustList = option2;
                    else
                        % Do nothing and the requested F and M will be scaled
                        %  down on the next iteration
                    end
                elseif(rankB==2) %Either two eqn. or two eqn. plus redundany copies -- can just choose two
                    a = B(1:2,:)\deltaT(1:2);
                    thrustList = thrustList + a(1)*self.nb1 + a(2)*self.nb2;
                else
                    fprintf('WARNING: rankB=rankAug > 2\n');
                end
            
            % If rank(B) < rank(augmentedMatrix), then the system is
            %  actually overconstrained
            elseif(rankB < rankAug)
                % Do nothing and F,M will be scaled down on the next
                %  iteration
            else
                fprintf('WARNING: rankB > rankAug. How???\n');
            end
            
            %fprintf('rank(B)=%d\trank(H)=%d\n',rank(B),rank(augmentedMatrix));
        end
        
        %------------------------------------------------------------------
        % Uses nullspace-traversal to extend the thrust envelope
        %   -Doesn't currently work correctly
        function [thrustList, numExceedsLimit] = nullMap_DEFUNCT(self,F_des,M_des, counter)
            fprintf('Down\n');
            if(nargin<4) counter = 1; end %First iteration of the recursion
            if(counter >= 200) assert(true, 'ERROR: Recursion limit reached'); end
            RECURSION_FACTOR = 0.99;
            thrustList = self.MAP_V2T*[shiftdim(F_des); shiftdim(M_des)]; %Map
            
            % Check if limits are exceeded
            exceedsLimitVector = self.whichThrustersExceedLimits(thrustList);
            numExceedsLimit = sum(abs(exceedsLimitVector)); %1-norm
            %fprintf('%d Thrusters exceed their limits.\n',numExceedsLimit);
            
            % If 1 thruster exceeds limits
            %   >Compute coefficient for just nb1 and just nb2; choose one
            if(numExceedsLimit == 1)
                index = abs((1:8)*exceedsLimitVector); %Picks out the one nonzero element
                limit = (exceedsLimitVector(index)==1)*self.THRUST_LIM_FORWARD...
                        +(exceedsLimitVector(index)==-1)*self.THRUST_LIM_BACKWARD;
                deltaT = limit - thrustList(index);
                option1 = thrustList + (deltaT/self.nb1(index))*self.nb1(index);
                option2 = thrustList + (deltaT/self.nb2(index))*self.nb2(index);
                thrustUtil1 = sum(abs(option1)); %i.e. norm(option1,1)
                thrustUtil2 = sum(abs(option2));
                isOption1Viable = sum(abs(self.whichThrustersExceedLimits(option1)))~=0;
                isOption2Viable = sum(abs(self.whichThrustersExceedLimits(option2)))~=0;
                if(isOption1Viable && (thrustUtil1 < thrustUtil2))
                    thrustList = option1;
                elseif(isOption2Viable)
                    thrustList = option2;
                else %Both options have at least one thruster exceeding a limit
                    % Recursively decrement the requested force/moment
                    [thrustList, numExceedsLimit] = nullMap(self, RECURSION_FACTOR*F_des, RECURSION_FACTOR*M_des, counter+1);
                end
            end
            
            % If 2 thrusters exceed limits
            %   >Set up and solve the fully constrained linear system
            %    for the coefficients on nb1 and nb2
            if(numExceedsLimit==2)
                % {T_exceed - T_lim} = [B]{a}
                T_exceed = []; % sub-set of thrustList (just the ones that exceed)
                T_lim = [];    % the corresponding positive or negative limit
                B = [];        % matrix -- the relevant components of the nb1,nb2 vectors
                for i = 1:8
                    if(exceedsLimitVector(i) ~= 0)
                        T_exceed = [T_exceed; thrustList(i)];
                        B = [B; self.nb1(i), self.nb2(i)];
                        limit = (exceedsLimitVector(i)==1)*self.THRUST_LIM_FORWARD...
                                +(exceedsLimitVector(i)==-1)*self.THRUST_LIM_BACKWARD;
                        T_lim = [T_lim; limit];
                    end
                end
                if(false) %det(B)<1e-10) %If the matrix is singular or close to it
                    fprintf('Singular matrix\n');
                    [thrustList, numExceedsLimit] = nullMap(self, RECURSION_FACTOR*F_des, RECURSION_FACTOR*M_des, counter+1);
                else
                    a = B\(T_lim-T_exceed); %Coefficients for nb1, nb2
                    thrustList = thrustList + a(1)*self.nb1 + a(2)*self.nb2;
                    % Check if new thrustList exceeds any limits; if so,
                    %  recursively decrement the requested force/moment
                    if(sum(abs(self.whichThrustersExceedLimits(thrustList)))~=0)
                        [thrustList, numExceedsLimit] = nullMap(self, RECURSION_FACTOR*F_des, RECURSION_FACTOR*M_des, counter+1);
                    end
                end
            end
            
            % If 3+ thrusters exceed limits
            %   >Scale down until only 2 exceed limits
            if(numExceedsLimit>2)
                [thrustList, numExceedsLimit] = nullMap(self, RECURSION_FACTOR*F_des, RECURSION_FACTOR*M_des, counter+1);
            end
            
            fprintf('Up\n');
            assert(sum(abs(self.whichThrustersExceedLimits(thrustList)))~=0, 'Something is wrotten in the state of Denmark');
        end
        
        %% RECOVERY FUNCTIONS =============================================
        function force = getForce(self, thrustList)
            force = self.map_T2F*thrustList;
        end
        function moment = getMoment(self, thrustList)
            moment = self.map_T2M*thrustList;
        end
        
        %% OTHER HELPER FUNCTIONS
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
        
        %% PLOTTING FUNCTIONS =============================================
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
        % This has only gotten more awkward
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
            maxLim = max([xlim,ylim]);
            %xlim(maxLim*[-1,1]);
            %ylim(maxLim*[-1,1]);
            view([-90,90]);
            xlabel(sprintf('x [1 in or %.2f N or %.2f Nm]',1/FORCE_SCALE,1/MOMENT_SCALE));
            ylabel(sprintf('y [1 in or %.2f N or %.2f Nm]',1/FORCE_SCALE,1/MOMENT_SCALE));
            zlabel(sprintf('z [1 in or %.2f N or %.2f Nm]',1/FORCE_SCALE,1/MOMENT_SCALE));
            if(sum(F_des)==0)
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